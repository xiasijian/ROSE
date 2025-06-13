# Import necessary standard library modules
import os  # For operating system interactions (file paths, etc.)
import re  # For regular expressions
import sys  # For system-specific parameters and functions
import subprocess  # For spawning new processes
import datetime  # For date/time operations
from collections import defaultdict  # Specialized dictionary with default values

#==================================================================
#==========================I/O FUNCTIONS===========================
#==================================================================

def unParseTable(table, output, sep):
    """Write a table (2D list) to a file with specified separator"""
    with open(output, 'w') as fh_out:  # Open output file in write mode
        if len(sep) == 0:  # If no separator specified
            for i in table:  # Write each element directly
                fh_out.write(str(i) + '\n')  # Convert to string and add newline
        else:  # If separator specified
            for line in table:  # Process each line
                line = [str(x) for x in line]  # Convert all elements to strings
                fh_out.write(sep.join(line) + '\n')  # Join with separator and write

def parseTable(fn, sep, header=False, excel=False):
    """Parse a delimited text file into a 2D list (table)"""
    with open(fn) as fh:  # Open input file
        if header:  # If header flag is True
            next(fh)  # Skip the header line
        
        table = []  # Initialize empty table
        for line in fh:  # Read each line
            line = line.rstrip().split(sep)  # Remove trailing whitespace and split
            table.append(line)  # Add to table
    
    return table  # Return the parsed table

def formatBed(bed, output=''):
    """Format BED file to standard 6-column format"""
    newBed = []  # Initialize output list
    
    if isinstance(bed, str):  # If input is a file path string
        bed = parseTable(bed, '\t')  # Parse it as tab-delimited

    indexTicker = 1  # Counter for unique IDs
    for line in bed:  # Process each BED line
        newLine = line[0:4]  # Take first 4 columns (chr, start, end, name)
        try:
            strand = line[5]  # Try to get strand from column 6
        except IndexError:  # If column 6 doesn't exist
            strand = '.'  # Default to unknown strand
        newLine += [indexTicker, strand]  # Add ID and strand
        indexTicker += 1  # Increment ID counter
        newBed.append(newLine)  # Add formatted line
    
    if output:  # If output file specified
        unParseTable(newBed, output, '\t')  # Write to file
    else:
        return newBed  # Otherwise return formatted BED

def bedToGFF(bed, output=''):
    """Convert BED format to GFF format"""
    if isinstance(bed, str):  # If input is file path
        bed = parseTable(bed, '\t')  # Parse it
    
    bed = formatBed(bed)  # Ensure proper BED format
    gff = []  # Initialize GFF output
    
    for line in bed:  # Convert each BED line to GFF
        # Create GFF line with proper fields:
        # 0:chr, 1:name, 2:empty, 3:start, 4:end, 5:ID, 6:strand, 7:empty, 8:name
        gffLine = [line[0], line[3], '', line[1], line[2], line[4], line[5], '', line[3]]
        gff.append(gffLine)

    if output:  # If output file specified
        unParseTable(gff, output, '\t')  # Write to file
    else:
        return gff  # Otherwise return GFF data

def gffToBed(gff, output=''):
    """Convert GFF format to BED format"""
    bed = []  # Initialize BED output
    indexTicker = 1  # Counter for unique IDs
    
    for line in gff:  # Process each GFF line
        # Create BED line with proper fields:
        # 0:chr, 1:start, 2:end, 3:name, 4:ID, 5:strand
        newLine = [line[0], line[3], line[4], line[1], indexTicker, line[6]]
        indexTicker += 1  # Increment ID
        bed.append(newLine)  # Add to output
    
    if output:  # If output file specified
        unParseTable(bed, output, '\t')  # Write to file
    else:
        return bed  # Otherwise return BED data

def formatFolder(folderName, create=False):
    """Ensure folder path ends with slash and optionally create it"""
    if not folderName.endswith('/'):  # If missing trailing slash
        folderName += '/'  # Add it
    
    if os.path.exists(folderName):  # If folder exists
        return folderName  # Return path
    else:  # If folder doesn't exist
        print(f'folder {folderName} does not exist')  # Notify user
        if create:  # If creation requested
            os.makedirs(folderName, exist_ok=True)  # Create folder
            return folderName  # Return path
        else:
            return False  # Indicate failure

#==================================================================
#===================ANNOTATION FUNCTIONS===========================
#==================================================================

def makeStartDict(annotFile, geneList=[]):
    """Create dictionary of gene start/end positions from annotation file"""
    if isinstance(geneList, str):  # If geneList is a file path
        geneList = parseTable(geneList, '\t')  # Parse it
        geneList = [line[0] for line in geneList]  # Extract gene names
            
    if 'REFSEQ' in annotFile.upper():  # If REFSEQ format
        refseqTable, refseqDict = importRefseq(annotFile)  # Parse refseq file
        if not geneList:  # If no gene list provided
            geneList = list(refseqDict.keys())  # Use all genes
            
        startDict = {}  # Initialize output dict
        for gene in geneList:  # Process each gene
            if gene not in refseqDict:  # Skip if gene not found
                continue
                
            # Store gene information in dictionary
            startDict[gene] = {
                'sense': refseqTable[refseqDict[gene][0]][3],  # Strand
                'chr': refseqTable[refseqDict[gene][0]][2],  # Chromosome
                'start': getTSSs([gene], refseqTable, refseqDict),  # TSS position
                'name': refseqTable[refseqDict[gene][0]][12]  # Gene name
            }
            
            # Set end position based on strand
            if startDict[gene]['sense'] == '+':
                startDict[gene]['end'] = [int(refseqTable[refseqDict[gene][0]][5])]
            else:
                startDict[gene]['end'] = [int(refseqTable[refseqDict[gene][0]][4])]
                
    return startDict  # Return gene position dictionary

def getTSSs(geneList, refseqTable, refseqDict):
    """Get transcription start sites for given genes"""
    if not geneList:  # If no genes specified
        refseq = refseqTable  # Use all refseq entries
    else:  # Otherwise
        refseq = refseqFromKey(geneList, refseqDict, refseqTable)  # Get specific genes
        
    TSS = []  # Initialize TSS list
    for line in refseq:  # For each gene entry
        if line[3] == '+':  # If on positive strand
            TSS.append(line[4])  # Start is TSS
        if line[3] == '-':  # If on negative strand
            TSS.append(line[5])  # End is TSS
            
    return [int(t) for t in TSS]  # Return TSS positions as integers

def refseqFromKey(refseqKeyList, refseqDict, refseqTable):
    """Get refseq entries matching given keys"""
    typeRefseq = []  # Initialize output list
    for name in refseqKeyList:  # For each requested gene
        if name in refseqDict:  # If gene exists in refseq
            # Add corresponding refseq entry to output
            typeRefseq.append(refseqTable[refseqDict[name][0]])
    return typeRefseq  # Return matching entries

def importRefseq(refseqFile, returnMultiples=False):
    """Parse refseq annotation file into table and dictionary"""
    refseqTable = parseTable(refseqFile, '\t')  # Parse as tab-delimited
    refseqDict = {}  # Initialize gene dictionary
    
    # Build dictionary of gene names to line numbers
    for i, line in enumerate(refseqTable[1:], 1):  # Skip header, count from 1
        if line[1] in refseqDict:  # If gene already seen
            refseqDict[line[1]].append(i)  # Add line number to list
        else:  # New gene
            refseqDict[line[1]] = [i]  # Create new entry

    # Find genes with multiple entries
    multiples = [i for i in refseqDict if len(refseqDict[i]) > 1]
    
    if returnMultiples:  # If requested
        return refseqTable, refseqDict, multiples  # Return with duplicates
    else:
        return refseqTable, refseqDict  # Return just table and dict

#==================================================================
#========================LOCUS INSTANCE============================
#==================================================================

class Locus:
    """Class representing a genomic region with chromosome, coordinates, and strand"""
    
    # Class-level dictionaries for memory optimization
    __chrDict = {}  # Shared chromosome string storage
    __senseDict = {'+': '+', '-': '-', '.': '.'}  # Strand options
    
    def __init__(self, chr, start, end, sense, ID=''):
        """Initialize locus with coordinates and metadata"""
        coords = sorted([int(start), int(end)])  # Ensure start <= end
        if chr not in self.__chrDict:  # If new chromosome
            self.__chrDict[chr] = chr  # Store chromosome string
        self._chr = self.__chrDict[chr]  # Reference shared string
        self._sense = self.__senseDict[sense]  # Validate strand
        self._start = coords[0]  # Start coordinate
        self._end = coords[1]  # End coordinate
        self._ID = ID  # Optional identifier
        
    # Property access methods
    def ID(self): return self._ID  # Get identifier
    def chr(self): return self._chr  # Get chromosome
    def start(self): return self._start  # Get start position
    def end(self): return self._end  # Get end position
    def len(self): return self._end - self._start + 1  # Calculate length
    
    def getAntisenseLocus(self):
        """Return new locus on opposite strand"""
        if self._sense == '.':  # If unstranded
            return self  # Return unchanged
        else:  # Otherwise
            switch = {'+': '-', '-': '+'}  # Strand switching dict
            return Locus(self._chr, self._start, self._end, switch[self._sense], self._ID)
            
    def coords(self): return [self._start, self._end]  # Return coordinate pair
    def sense(self): return self._sense  # Return strand
    
    def overlaps(self, otherLocus):
        """Check if two loci overlap on same strand"""
        if self.chr() != otherLocus.chr():  # Different chromosomes
            return False
        elif not (self._sense == '.' or  # Either unstranded or
                 otherLocus.sense() == '.' or  # other unstranded or
                 self.sense() == otherLocus.sense()):  # same strand
            return False
        elif self.start() > otherLocus.end() or otherLocus.start() > self.end():  # No overlap
            return False
        else:  # Otherwise overlaps
            return True
            
    def contains(self, otherLocus):
        """Check if this locus completely contains another"""
        if self.chr() != otherLocus.chr():  # Different chromosomes
            return False
        elif not (self._sense == '.' or  # Strand check
                 otherLocus.sense() == '.' or 
                 self.sense() == otherLocus.sense()):
            return False
        elif self.start() > otherLocus.start() or otherLocus.end() > self.end():  # Containment
            return False
        else:
            return True
            
    def overlapsAntisense(self, otherLocus):
        """Check overlap on opposite strands"""
        return self.getAntisenseLocus().overlaps(otherLocus)
        
    def containsAntisense(self, otherLocus):
        """Check containment on opposite strands"""
        return self.getAntisenseLocus().contains(otherLocus)
        
    def __hash__(self):
        """Hash function for use in dictionaries/sets"""
        return self._start + self._end  # Simple hash based on coordinates
        
    def __eq__(self, other):
        """Equality comparison for loci"""
        if self.__class__ != other.__class__:  # Different classes
            return False
        if self.chr() != other.chr():  # Different chromosomes
            return False
        if self.start() != other.start():  # Different starts
            return False
        if self.end() != other.end():  # Different ends
            return False
        if self.sense() != other.sense():  # Different strands
            return False
        return True  # Otherwise equal
        
    def __ne__(self, other):
        """Inequality comparison"""
        return not(self.__eq__(other))
        
    def __str__(self):
        """String representation of locus"""
        return f"{self.chr()}({self.sense()}):{self._start}-{self._end}"
        
    def checkRep(self):
        """Check representation invariants (currently unused)"""
        pass

class LocusCollection:
    """Efficient collection of loci with spatial indexing"""
    
    def __init__(self, loci, windowSize):
        """Initialize collection with loci and window size for spatial indexing"""
        self.__chrToCoordToLoci = {}  # Nested dict for spatial lookup
        self.__loci = {}  # Dictionary of all loci
        self.__winSize = windowSize  # Size of spatial bins
        for lcs in loci:  # Add all initial loci
            self.__addLocus(lcs)

    def __addLocus(self, lcs):
        """Internal method to add a locus to the collection"""
        if lcs not in self.__loci:  # If new locus
            self.__loci[lcs] = None  # Add to master dictionary
            
            # Determine chromosome/strand combinations to index
            if lcs.sense() == '.':  # If unstranded
                chrKeyList = [lcs.chr() + '+', lcs.chr() + '-']  # Index both strands
            else:  # If stranded
                chrKeyList = [lcs.chr() + lcs.sense()]  # Index only this strand
                
            for chrKey in chrKeyList:  # For each chromosome/strand combo
                if chrKey not in self.__chrToCoordToLoci:  # If new chromosome/strand
                    self.__chrToCoordToLoci[chrKey] = {}  # Initialize dict
                    
                # Add to all spatial bins this locus overlaps
                for n in self.__getKeyRange(lcs):
                    if n not in self.__chrToCoordToLoci[chrKey]:
                        self.__chrToCoordToLoci[chrKey][n] = []  # Initialize bin
                    self.__chrToCoordToLoci[chrKey][n].append(lcs)  # Add to bin

    def __getKeyRange(self, locus):
        """Calculate which spatial bins a locus falls into"""
        start = locus.start() // self.__winSize  # Start bin
        end = locus.end() // self.__winSize + 1  # End bin (inclusive)
        return range(start, end)  # All bins covered

    def __len__(self):
        """Number of loci in collection"""
        return len(self.__loci)
        
    def append(self, new):
        """Add a single locus to collection"""
        self.__addLocus(new)
        
    def extend(self, newList):
        """Add multiple loci to collection"""
        for lcs in newList:
            self.__addLocus(lcs)
            
    def hasLocus(self, locus):
        """Check if collection contains a locus"""
        return locus in self.__loci
        
    def remove(self, old):
        """Remove a locus from the collection"""
        if old not in self.__loci:  # If locus not present
            raise ValueError("requested locus isn't in collection")
            
        del self.__loci[old]  # Remove from master dict
        
        # Determine chromosome/strand combinations to update
        if old.sense() == '.':
            senseList = ['+', '-']  # Both strands if unstranded
        else:
            senseList = [old.sense()]  # Only this strand
            
        # Remove from all spatial bins
        for k in self.__getKeyRange(old):
            for sense in senseList:
                chrKey = old.chr() + sense
                if k in self.__chrToCoordToLoci[chrKey]:
                    try:
                        self.__chrToCoordToLoci[chrKey][k].remove(old)  # Remove from bin
                    except ValueError:  # If not present (shouldn't happen)
                        pass

    def getWindowSize(self):
        """Get the spatial bin window size"""
        return self.__winSize
        
    def getLoci(self):
        """Get all loci in collection"""
        return list(self.__loci.keys())
        
    def getChrList(self):
        """Get list of all chromosomes in collection"""
        tempKeys = {}
        for k in self.__chrToCoordToLoci.keys():
            tempKeys[k[:-1]] = None  # Remove strand suffix
        return list(tempKeys.keys())
            
    def __subsetHelper(self, locus, sense):
        """Internal method to find loci matching strand criteria"""
        sense = sense.lower()  # Normalize sense parameter
        if sense not in ['sense', 'antisense', 'both']:  # Validate
            raise ValueError(f"sense command invalid: '{sense}'.")
            
        matches = {}  # Dictionary for results (preserves uniqueness)
        senses = ['+', '-']  # Possible strands
        
        # Define filter function based on sense parameter
        if locus.sense() == '.' or sense == 'both':
            def lamb(s): return True  # No strand filtering
        elif sense == 'sense':
            def lamb(s): return s == locus.sense()  # Same strand
        elif sense == 'antisense':
            def lamb(s): return s != locus.sense()  # Opposite strand
        else:
            raise ValueError(f"sense value was inappropriate: '{sense}'.")
            
        # Find all loci in overlapping bins that match strand criteria
        for s in filter(lamb, senses):  # For each relevant strand
            chrKey = locus.chr() + s  # Construct chromosome/strand key
            if chrKey in self.__chrToCoordToLoci:  # If chromosome exists
                for n in self.__getKeyRange(locus):  # For each relevant bin
                    if n in self.__chrToCoordToLoci[chrKey]:
                        for lcs in self.__chrToCoordToLoci[chrKey][n]:
                            matches[lcs] = None  # Add to matches
                            
        return list(matches.keys())  # Return unique loci
        
    def getOverlap(self, locus, sense='sense'):
        """Get loci overlapping the query locus"""
        matches = self.__subsetHelper(locus, sense)  # Get candidate loci
        realMatches = {}  # Dictionary for actual overlaps
        
        # Check for real overlaps (not just bin co-membership)
        if sense in ['sense', 'both']:
            for i in [lcs for lcs in matches if lcs.overlaps(locus)]:
                realMatches[i] = None
                
        if sense in ['antisense', 'both']:
            for i in [lcs for lcs in matches if lcs.overlapsAntisense(locus)]:
                realMatches[i] = None
                
        return list(realMatches.keys())  # Return overlapping loci

    def getContained(self, locus, sense='sense'):
        """Get loci completely contained within query locus"""
        matches = self.__subsetHelper(locus, sense)
        realMatches = {}
        
        if sense in ['sense', 'both']:
            for i in [lcs for lcs in matches if locus.contains(lcs)]:
                realMatches[i] = None
                
        if sense in ['antisense', 'both']:
            for i in [lcs for lcs in matches if locus.containsAntisense(lcs)]:
                realMatches[i] = None
                
        return list(realMatches.keys())

    def getContainers(self, locus, sense='sense'):
        """Get loci that completely contain query locus"""
        matches = self.__subsetHelper(locus, sense)
        realMatches = {}
        
        if sense in ['sense', 'both']:
            for i in [lcs for lcs in matches if lcs.contains(locus)]:
                realMatches[i] = None
                
        if sense in ['antisense', 'both']:
            for i in [lcs for lcs in matches if lcs.containsAntisense(locus)]:
                realMatches[i] = None
                
        return list(realMatches.keys())

    def stitchCollection(self, stitchWindow=1, sense='both'):
        """Merge overlapping loci into contiguous regions"""
        locusList = self.getLoci()  # Get all loci
        oldCollection = LocusCollection(locusList, 500)  # Temporary collection
        stitchedCollection = LocusCollection([], 500)  # Output collection

        for locus in locusList:  # Process each locus
            if oldCollection.hasLocus(locus):  # If not already processed
                oldCollection.remove(locus)  # Mark as processed
                # Find overlapping loci within stitch window
                overlappingLoci = oldCollection.getOverlap(
                    Locus(locus.chr(), locus.start()-stitchWindow, 
                         locus.end()+stitchWindow, locus.sense(), locus.ID()),
                    sense
                )
                
                stitchTicker = 1  # Count merged loci
                while overlappingLoci:  # While there are overlaps
                    stitchTicker += len(overlappingLoci)  # Increment count
                    overlapCoords = locus.coords()  # Current coordinates
                    
                    # Merge all overlapping loci
                    for overlappingLocus in overlappingLoci:
                        overlapCoords += overlappingLocus.coords()  # Add their coords
                        oldCollection.remove(overlappingLocus)  # Mark as processed
                        
                    # Create new merged locus
                    if sense == 'both':
                        locus = Locus(
                            locus.chr(),
                            min(overlapCoords),  # New start
                            max(overlapCoords),  # New end
                            '.',  # Unstranded
                            locus.ID()
                        )
                    else:
                        locus = Locus(
                            locus.chr(),
                            min(overlapCoords),
                            max(overlapCoords),
                            locus.sense(),
                            locus.ID()
                        )
                        
                    # Find new overlaps with merged locus
                    overlappingLoci = oldCollection.getOverlap(
                        Locus(locus.chr(), locus.start()-stitchWindow, 
                             locus.end()+stitchWindow, locus.sense()),
                        sense
                    )
                
                # Set ID indicating number of merged loci
                locus._ID = f"{stitchTicker}_{locus.ID()}_lociStitched"
                stitchedCollection.append(locus)  # Add to output
                
        return stitchedCollection  # Return merged loci

#==================================================================
#========================LOCUS FUNCTIONS===========================
#==================================================================

def locusCollectionToGFF(locusCollection):
    """Convert LocusCollection to GFF format"""
    lociList = locusCollection.getLoci()  # Get all loci
    gff = []  # Initialize output
    
    for locus in lociList:  # Convert each locus to GFF line
        newLine = [
            locus.chr(),  # Chromosome
            locus.ID(),  # Feature name
            '',  # Empty field
            locus.coords()[0],  # Start
            locus.coords()[1],  # End
            '',  # Empty field
            locus.sense(),  # Strand
            '',  # Empty field
            locus.ID()  # Attributes
        ]
        gff.append(newLine)  # Add to output
        
    return gff  # Return GFF data

def gffToLocusCollection(gff, window=500):
    """Convert GFF file to LocusCollection"""
    lociList = []  # Initialize loci list
    
    if isinstance(gff, str):  # If input is file path
        gff = parseTable(gff, '\t')  # Parse it

    nameList = []  # For checking unique names
    for line in gff:  # Process each GFF line
        if not line or line[0].startswith('#'):  # Skip empty/comments
            continue
            
        if len(line) < 7:  # If line has too few fields
            print('SKIPPING THIS LINE')
            print(line)
            continue

        # Determine feature name - priority to column 2, then 9, then construct
        name = line[1] if line[1] else (line[8] if len(line) > 8 and line[8] else 
              f"{line[0]}:{line[6]}:{line[3]}-{line[4]}")
              
        nameList.append(name)  # Add to name list
        # Create locus from GFF fields
        lociList.append(Locus(line[0], line[3], line[4], line[6], name))

    # Verify all names are unique
    if len(nameList) != len(uniquify(nameList)):
        print('ERROR: FOR GFFS, ALL REGIONS MUST HAVE A UNIQUE IDENTIFIER IN COLUMN 2')
        sys.exit()
        
    return LocusCollection(lociList, window)  # Return collection

def makeTranscriptCollection(annotFile, upSearch, downSearch, window=500, geneList=[]):
    """Create LocusCollection of transcripts with search windows"""
    if 'REFSEQ' in annotFile.upper():  # If refseq format
        refseqTable, refseqDict = importRefseq(annotFile)  # Parse file
        locusList = []  # Initialize output
        
        if not geneList:  # If no gene list provided
            geneList = list(refseqDict.keys())  # Use all genes
            
        for line in refseqTable[1:]:  # Skip header
            if line[1] in geneList:  # If gene in target list
                if line[3] == '-':  # Negative strand
                    # Create locus with downstream extension
                    locus = Locus(line[2], int(line[4])-downSearch, 
                                int(line[5])+upSearch, line[3], line[1])
                else:  # Positive strand
                    # Create locus with upstream extension
                    locus = Locus(line[2], int(line[4])-upSearch, 
                                int(line[5])+downSearch, line[3], line[1])
                locusList.append(locus)  # Add to list

    return LocusCollection(locusList, window)  # Return collection

def makeTSSLocus(gene, startDict, upstream, downstream):
    """Create locus around transcription start site"""
    start = startDict[gene]['start'][0]  # Get TSS position
    
    if startDict[gene]['sense'] == '-':  # Negative strand
        return Locus(startDict[gene]['chr'], start-downstream, 
                    start+upstream, '-', gene)
    else:  # Positive strand
        return Locus(startDict[gene]['chr'], start-upstream, 
                    start+downstream, '+', gene)

def makeSearchLocus(locus, upSearch, downSearch):
    """Extend locus by specified amounts on each side"""
    if locus.sense() == '-':  # Negative strand
        return Locus(locus.chr(), locus.start()-downSearch, 
                    locus.end()+upSearch, locus.sense(), locus.ID())
    else:  # Positive strand
        return Locus(locus.chr(), locus.start()-upSearch, 
                    locus.end()+downSearch, locus.sense(), locus.ID())

#==================================================================
#==========================BAM CLASS===============================
#==================================================================

def checkChrStatus(bamFile):
    """Check if BAM file uses 'chr' prefix in chromosome names"""
    command = f'samtools view {bamFile} | head -n 1'  # Command to get first read
    stats = subprocess.Popen(command, stdin=subprocess.PIPE, 
                            stderr=subprocess.PIPE, 
                            stdout=subprocess.PIPE, shell=True, text=True)
    statLines = stats.stdout.readlines()  # Get output
    stats.stdout.close()  # Close pipe
    
    chrPattern = re.compile('chr')  # Pattern to check
    for line in statLines:  # Check each line
        sline = line.split("\t")  # Split SAM fields
        if re.search(chrPattern, sline[2]):  # If chromosome has 'chr'
            return 1  # Has chr prefix
        else:
            return 0  # No chr prefix

def convertBitwiseFlag(flag):
    """Convert SAM flag to strand ('+' or '-')"""
    if int(flag) & 16:  # If reverse strand bit set
        return "-"  # Negative strand
    else:
        return "+"  # Positive strand

class Bam:
    """Class for working with BAM alignment files"""
    
    def __init__(self, bamFile):
        """Initialize with BAM file path"""
        self._bam = bamFile  # Store file path

    def getTotalReads(self, readType='mapped'):
        """Get total or mapped reads from BAM"""
        command = f'samtools flagstat {self._bam}'  # Command to get stats
        stats = subprocess.Popen(command, stdin=subprocess.PIPE, 
                                stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE, shell=True, text=True)
        statLines = stats.stdout.readlines()  # Get output
        stats.stdout.close()  # Close pipe
        
        if readType == 'mapped':  # If mapped reads requested
            for line in statLines:  # Find mapped count line
                if 'mapped (' in line:
                    return int(line.split(' ')[0])  # Return count
        elif readType == 'total':  # If total reads requested
            return int(statLines[0].split(' ')[0])  # First line is total

    def getRawReads(self, locus, sense, unique=False, includeJxnReads=False, printCommand=False):
        """Get raw reads overlapping a locus"""
        locusLine = f"{locus.chr()}:{locus.start()}-{locus.end()}"  # Region string
        command = f"samtools view {self._bam} {locusLine}"  # samtools command
        
        if printCommand:  # If debugging
            print(command)  # Show command
            
        getReads = subprocess.Popen(command, stdin=subprocess.PIPE, 
                                   stderr=subprocess.PIPE,
                                   stdout=subprocess.PIPE, shell=True, text=True)
        reads = getReads.communicate()[0].split('\n')[:-1]  # Get and split output
        reads = [read.split('\t') for read in reads]  # Split each read into fields
        
        if not includeJxnReads:  # If excluding junction reads
            reads = [read for read in reads if 'N' not in read[5]]  # Filter by CIGAR

        if sense == '-':  # If antisense requested
            strand = '+' if locus.sense() == '-' else '-'  # Opposite of locus
        else:
            strand = locus.sense()  # Same as locus
            
        keptReads = []  # Initialize output
        seqDict = defaultdict(int)  # For tracking duplicates
        
        for read in reads:  # Process each read
            readStrand = convertBitwiseFlag(read[1])  # Get strand from flag
            
            if sense in ['both', '.'] or readStrand == strand:  # If strand matches
                if unique and seqDict[read[9]] == 0:  # If unique requested and new
                    keptReads.append(read)  # Add to output
                elif not unique:  # If not checking uniqueness
                    keptReads.append(read)  # Add to output
                    
            seqDict[read[9]] += 1  # Count sequence occurrences
            
        return keptReads  # Return filtered reads

    def readsToLoci(self, reads, IDtag=''):
        """Convert SAM format reads to Locus objects"""
        if not IDtag or IDtag not in ['sequence', 'seqID', 'none']:  # Validate
            print('please specify one of the three options: sequence, seqID, none')
            return []
            
        loci = []  # Initialize output
        numPattern = re.compile(r'\d*')  # For parsing CIGAR
        
        for read in reads:  # Process each read
            chrom = read[2]  # Chromosome
            strand = convertBitwiseFlag(read[1])  # Strand
            
            if IDtag == 'sequence':  # If using sequence as ID
                ID = read[9]  # Sequence from SAM
            elif IDtag == 'seqID':  # If using read name as ID
                ID = read[0]  # Read name
            else:  # If no ID
                ID = ''
                
            length = len(read[9])  # Sequence length
            start = int(read[3])  # Start position
            
            if 'N' in read[5]:  # If junction read (has skip in CIGAR)
                parts = [int(x) for x in filter(None, re.findall(numPattern, read[5]))][0:3]
                if len(parts) == 3:  # If proper junction
                    first, gap, second = parts  # Get segments
                    
                    if IDtag == 'sequence':  # If using sequence as ID
                        # Split into two loci with sequence portions
                        loci.append(Locus(chrom, start, start+first, strand, ID[0:first]))
                        loci.append(Locus(chrom, start+first+gap, 
                                        start+first+gap+second, strand, ID[first:]))
                    else:  # If not using sequence as ID
                        # Split into two loci with same ID
                        loci.append(Locus(chrom, start, start+first, strand, ID))
                        loci.append(Locus(chrom, start+first+gap, 
                                        start+first+gap+second, strand, ID))
            else:  # If continuous read
                loci.append(Locus(chrom, start, start+length, strand, ID))  # Single locus
                
        return loci  # Return loci

    def getReadsLocus(self, locus, sense='both', unique=True, IDtag='', includeJxnReads=False):
        """Get loci for reads overlapping a region"""
        reads = self.getRawReads(locus, sense, unique, includeJxnReads)  # Get reads
        return self.readsToLoci(reads, IDtag)  # Convert to loci

    def getReadSequences(self, locus, sense='both', unique=True, includeJxnReads=False):
        """Get sequences of reads overlapping a region"""
        reads = self.getRawReads(locus, sense, unique, includeJxnReads)  # Get reads
        return [read[9] for read in reads]  # Extract sequences
    
    def getReadStarts(self, locus, sense='both', unique=False, includeJxnReads=False):
        """Get start positions of reads overlapping a region"""
        reads = self.getRawReads(locus, sense, unique, includeJxnReads)  # Get reads
        return [int(read[3]) for read in reads]  # Extract starts
        
    def getReadCount(self, locus, sense='both', unique=True, includeJxnReads=False):
        """Count reads overlapping a region"""
        reads = self.getRawReads(locus, sense, unique, includeJxnReads)  # Get reads
        return len(reads)  # Return count

#==================================================================
#========================MISC FUNCTIONS============================
#==================================================================

def uniquify(seq, idfun=None):
    """Remove duplicates from sequence while preserving order"""
    if idfun is None:  # If no key function
        def idfun(x): return x  # Use element itself
        
    seen = {}  # Dictionary for tracking seen items
    result = []  # Output list
    
    for item in seq:  # Process each item
        marker = idfun(item)  # Get key
        if marker in seen:  # If already seen
            continue  # Skip
        seen[marker] = 1  # Mark as seen
        result.append(item)  # Add to output
        
    return result  # Return unique items

def order(x, NoneIsLast=True, decreasing=False):
    """Return indices that would sort a list"""
    n = len(x)  # Get length
    ix = list(range(n))  # Initial indices
    
    if None not in x:  # If no None values
        ix.sort(key=lambda j: x[j], reverse=decreasing)  # Sort by values
    else:  # If None values present
        def key(i):  # Key function handling None
            elem = x[i]
            if decreasing == NoneIsLast:  # If None should be last
                return not(elem is None), elem  # Sort None first
            else:  # If None should be first
                return elem is None, elem  # Sort None last
                
        ix.sort(key=key, reverse=decreasing)  # Sort with key
    
    if NoneIsLast is None:  # If excluding None
        n = len([i for i in ix if x[i] is not None])  # Count non-None
        return ix[:n]  # Return only those indices
        
    return ix  # Return all indices
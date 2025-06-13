# bamToGFF.py
# script to grab reads from a bam that align to a .gff file
import sys
import re
import ROSE_utils_new
from collections import defaultdict
import os

#=====================================================================
#====================MAPPING BAM READS TO GFF REGIONS=================
#=====================================================================

def mapBamToGFF(bamFile, gff, sense='both', extension=200, floor=0, rpm=False, matrix=None):
    '''maps reads from a bam to a gff'''
    floor = int(floor)
    
    # USING BAM CLASS
    bam = ROSE_utils_new.Bam(bamFile)

    # new GFF to write to
    newGFF = []
    
    # millionMappedReads
    if rpm:    
        MMR = round(float(bam.getTotalReads('mapped'))/1000000, 4)
    else:
        MMR = 1

    print(f'using a MMR value of {MMR}')
    
    if ROSE_utils_new.checkChrStatus(bamFile) == 1:
        print("has chr")
        hasChrFlag = 1
    else:
        print("does not have chr")
        hasChrFlag = 0
      
    if isinstance(gff, str):
        gff = ROSE_utils_new.parseTable(gff, '\t')
        
    # setting up a matrix table
    if matrix:
        bamName = os.path.basename(bamFile)
        newGFF.append(['GENE_ID', 'locusLine'] + [f'bin_{n}_{bamName}' for n in range(1, int(matrix)+1)])        

    # getting and processing reads for gff lines
    ticker = 0
    print('Number lines processed')
    for line in gff:
        line = line[0:9]
        if ticker % 100 == 0:
            print(ticker)
        ticker += 1
        
        if not hasChrFlag:
            line[0] = re.sub(r"chr", r"", line[0])
            
        gffLocus = ROSE_utils_new.Locus(line[0], int(line[3]), int(line[4]), line[6], line[1])
        searchLocus = ROSE_utils_new.makeSearchLocus(gffLocus, int(extension), int(extension))
        
        reads = bam.getReadsLocus(searchLocus, 'both', False, 'none')
        
        # now extend the reads and make a list of extended reads
        extendedReads = []
        for locus in reads:
            if locus.sense() in ['+', '.']:
                newLocus = ROSE_utils_new.Locus(locus.chr(), locus.start(), locus.end()+extension, locus.sense(), locus.ID())
            else:
                newLocus = ROSE_utils_new.Locus(locus.chr(), locus.start()-extension, locus.end(), locus.sense(), locus.ID())
            extendedReads.append(newLocus)
            
        if gffLocus.sense() in ['+', '.']:
            senseReads = [x for x in extendedReads if x.sense() in ['+', '.']]
            antiReads = [x for x in extendedReads if x.sense() == '-']
        else:
            senseReads = [x for x in extendedReads if x.sense() in ['-', '.']]
            antiReads = [x for x in extendedReads if x.sense() == '+']

        senseHash = defaultdict(int)
        antiHash = defaultdict(int)

        # filling in the readHashes             
        if sense in ['+', 'both', '.']:
            for read in senseReads:
                for x in range(read.start(), read.end()+1):
                    senseHash[x] += 1
                    
        if sense in ['-', 'both', '.']:
            for read in antiReads:
                for x in range(read.start(), read.end()+1):
                    antiHash[x] += 1

        # now apply flooring and filtering for coordinates
        keys = list(set(list(senseHash.keys()) + list(antiHash.keys())))
        if floor > 0:
            keys = [x for x in keys if (senseHash[x] + antiHash[x]) > floor]
            
        # coordinate filtering
        keys = [x for x in keys if gffLocus.start() < x < gffLocus.end()]

        # setting up the output table
        clusterLine = [gffLocus.ID(), str(gffLocus)]

        if matrix:
            # getting the binsize
            binSize = (gffLocus.len()-1)/int(matrix)
            nBins = int(matrix)
            if binSize == 0:
                clusterLine += ['NA'] * int(matrix)
                newGFF.append(clusterLine)
                continue
                
            n = 0
            if gffLocus.sense() in ['+', '.', 'both']:
                i = gffLocus.start()
                while n < nBins:
                    n += 1
                    binKeys = [x for x in keys if i < x < i+binSize]
                    binDen = float(sum([senseHash[x] + antiHash[x] for x in binKeys]))/binSize
                    clusterLine.append(round(binDen/MMR, 4))
                    i += binSize
            else:
                i = gffLocus.end()
                while n < nBins:
                    n += 1
                    binKeys = [x for x in keys if i-binSize < x < i]
                    binDen = float(sum([senseHash[x] + antiHash[x] for x in binKeys]))/binSize
                    clusterLine.append(round(binDen/MMR, 4))
                    i -= binSize
                    
        newGFF.append(clusterLine)
        
    return newGFF
        
#=====================================================================
#============================MAIN METHOD==============================
#=====================================================================

def main():
    from optparse import OptionParser
    usage = "usage: %prog [options] -b [SORTED BAMFILE] -i [INPUTFILE] -o [OUTPUTFILE]"
    parser = OptionParser(usage=usage)
    
    # required flags
    parser.add_option("-b", "--bam", dest="bam", nargs=1, default=None,
                    help="Enter .bam file to be processed.")
    parser.add_option("-i", "--input", dest="input", nargs=1, default=None,
                    help="Enter .gff or ENRICHED REGION file to be processed.")
    
    # output flag
    parser.add_option("-o", "--output", dest="output", nargs=1, default=None,
                    help="Enter the output filename.")
    
    # additional options
    parser.add_option("-s", "--sense", dest="sense", nargs=1, default='both',
                    help="Map to '+','-' or 'both' strands. Default maps to both.")
    parser.add_option("-f", "--floor", dest="floor", nargs=1, default=0,
                    help="Sets a read floor threshold necessary to count towards density")    
    parser.add_option("-e", "--extension", dest="extension", nargs=1, default=200,
                    help="Extends reads by n bp. Default value is 200bp")
    parser.add_option("-r", "--rpm", dest="rpm", action='store_true', default=False,
                    help="Normalizes density to reads per million (rpm)")
    parser.add_option("-m", "--matrix", dest="matrix", nargs=1, default=None,
                    help="Outputs a variable bin sized matrix. User must specify number of bins.")

    (options, args) = parser.parse_args()

    print(options)
    print(args)

    if options.bam:
        bamFile = options.bam
        fullPath = os.path.abspath(bamFile)
        bamName = os.path.splitext(os.path.basename(bamFile))[0]
        pathFolder = os.path.dirname(fullPath)
        fileList = os.listdir(pathFolder)
        hasBai = any(fileName.startswith(bamName) and fileName.endswith('.bai') for fileName in fileList)

        if not hasBai:
            print('ERROR: no associated .bai file found with bam. Must use a sorted bam with accompanying index file')
            parser.print_help()
            sys.exit()
   
    if options.sense and options.sense not in ['+', '-', '.', 'both']:
        print('ERROR: sense flag must be followed by +,-,.,both')
        parser.print_help()
        sys.exit()

    if options.matrix:
        try:
            int(options.matrix)
        except ValueError:
            print('ERROR: User must specify an integer bin number for matrix (try 50)')
            parser.print_help()
            sys.exit()
            
    if options.input and options.bam:
        inputFile = options.input
        gffFile = inputFile
        bamFile = options.bam
        
        if options.output is None:
            output = os.path.join(os.getcwd(), os.path.basename(inputFile) + '.mapped')
        else:
            output = options.output
            
        if options.matrix:
            print('mapping to GFF and making a matrix with fixed bin number')
            newGFF = mapBamToGFF(bamFile, gffFile, options.sense, int(options.extension), 
                               options.floor, options.rpm, options.matrix)
        else:
            newGFF = mapBamToGFF(bamFile, gffFile, options.sense, int(options.extension), 
                               options.floor, options.rpm, None)
        
        ROSE_utils_new.unParseTable(newGFF, output, '\t')
    else:
        parser.print_help()
                
if __name__ == "__main__":
    main()
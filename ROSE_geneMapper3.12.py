#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ROSE_geneMapper.py - Maps genes to enhancer regions from ROSE_Main output

Functionality:
1. Takes enhancer region table from ROSE_Main and maps genes to each enhancer
2. Creates two output files:
   - Enhancer-to-gene table (each row represents an enhancer)
   - Gene-to-enhancer table (each row represents a gene)
3. By default processes super-enhancers only

Usage:
python ROSE_geneMapper.py -i <input_file> -g <genome_build> [options]
"""

import sys
import os
from collections import defaultdict
from typing import List, Dict, Tuple, DefaultDict

# Import ROSE utility functions
import ROSE_utils


def mapEnhancerToGene(annotFile: str, 
                     enhancerFile: str, 
                     transcribedFile: str = '', 
                     uniqueGenes: bool = True, 
                     searchWindow: int = 50000, 
                     noFormatTable: bool = False) -> Tuple[List[List[str]], List[List[str]]]:
    """
    Maps genes to enhancer regions based on genomic coordinates
    
    Args:
        annotFile: Path to genome annotation file
        enhancerFile: Path to enhancer region file from ROSE
        transcribedFile: Optional file with transcribed genes to consider
        uniqueGenes: If True, returns only unique gene names
        searchWindow: Window size (bp) to search for nearby genes
        noFormatTable: Preserve original table formatting if True
        
    Returns:
        Tuple of two tables:
        1. Enhancer-to-gene mapping table
        2. Gene-to-enhancer mapping table
    """
    
    # Create dictionary of gene start positions
    startDict = ROSE_utils.makeStartDict(annotFile)
    enhancerTable = ROSE_utils.parseTable(enhancerFile, '\t')

    # Get list of transcribed genes to consider
    if transcribedFile:
        transcribedTable = ROSE_utils.parseTable(transcribedFile, '\t')
        transcribedGenes = [line[1] for line in transcribedTable]
    else:
        transcribedGenes = list(startDict.keys())

    print('Creating transcript collection')
    transcribedCollection = ROSE_utils.makeTranscriptCollection(
        annotFile, 0, 0, 500, transcribedGenes)

    print('Creating TSS collection')
    tssLoci = []
    for geneID in transcribedGenes:
        tssLoci.append(ROSE_utils.makeTSSLocus(geneID, startDict, 0, 0))

    # Create collection of transcription start sites
    tssCollection = ROSE_utils.LocusCollection(tssLoci, 50)

    # Initialize dictionaries to store gene-enhancer relationships
    geneDict = {'overlapping': defaultdict(list), 'proximal': defaultdict(list)}
    rankDict = defaultdict(list)  # Stores enhancer ranks
    superDict = defaultdict(list)  # Stores super-enhancer status
    overallGeneList = []  # All genes found in analysis

    # Initialize output tables
    if noFormatTable:
        enhancerToGeneTable = [enhancerTable[0] + ['OVERLAP_GENES', 'PROXIMAL_GENES', 'CLOSEST_GENE']]
    else:
        enhancerToGeneTable = [enhancerTable[0][0:9] + ['OVERLAP_GENES', 'PROXIMAL_GENES', 'CLOSEST_GENE'] + enhancerTable[5][-2:]]

    geneToEnhancerTable = [['GENE_NAME', 'REFSEQ_ID', 'PROXIMAL_ENHANCERS', 'ENHANCER_RANKS', 'IS_SUPER']]

    # Process each enhancer region
    for line in enhancerTable:
        if line[0].startswith(('#', 'R')):
            continue

        # Create enhancer locus object
        enhancerString = f"{line[1]}:{line[2]}-{line[3]}"
        enhancerLocus = ROSE_utils.Locus(line[1], line[2], line[3], '.', line[0])

        # Find overlapping genes (genes whose transcripts overlap enhancer)
        overlappingLoci = transcribedCollection.getOverlap(enhancerLocus, 'both')
        overlappingGenes = [locus.ID() for locus in overlappingLoci]

        # Find proximal genes (TSS within search window of enhancer)
        proximalLoci = tssCollection.getOverlap(
            ROSE_utils.makeSearchLocus(enhancerLocus, searchWindow, searchWindow), 'both')
        proximalGenes = [locus.ID() for locus in proximalLoci]

        # Find distal genes (TSS within 1Mb of enhancer)
        distalLoci = tssCollection.getOverlap(
            ROSE_utils.makeSearchLocus(enhancerLocus, 1000000, 1000000), 'both')
        distalGenes = [locus.ID() for locus in distalLoci]

        # Remove duplicates between gene categories
        overlappingGenes = ROSE_utils.uniquify(overlappingGenes)
        proximalGenes = ROSE_utils.uniquify(proximalGenes)
        distalGenes = ROSE_utils.uniquify(distalGenes)
        
        # Remove genes that appear in multiple categories
        allEnhancerGenes = overlappingGenes + proximalGenes + distalGenes
        proximalGenes = [g for g in proximalGenes if g not in overlappingGenes]
        distalGenes = [g for g in distalGenes if g not in overlappingGenes + proximalGenes]

        # Find closest gene to enhancer center
        closestGene = ''
        if allEnhancerGenes:
            enhancerCenter = (int(line[2]) + int(line[3])) // 2
            distList = [abs(enhancerCenter - startDict[geneID]['start'][0]) 
                       for geneID in allEnhancerGenes]
            closestGene = startDict[allEnhancerGenes[distList.index(min(distList))]]['name']

        # Add row to enhancer-to-gene table
        if noFormatTable:
            newLine = list(line)
            newLine.append(','.join(ROSE_utils.uniquify([startDict[x]['name'] for x in overlappingGenes])))
            newLine.append(','.join(ROSE_utils.uniquify([startDict[x]['name'] for x in proximalGenes])))
            newLine.append(closestGene)
        else:
            newLine = line[0:9]
            newLine.append(','.join(ROSE_utils.uniquify([startDict[x]['name'] for x in overlappingGenes])))
            newLine.append(','.join(ROSE_utils.uniquify([start_dict[x]['name'] for x in proximalGenes])))
            newLine.append(closestGene)
            newLine += line[-2:]
        
        enhancerToGeneTable.append(newLine)

        # Record gene-enhancer relationships
        overallGeneList.extend(overlappingGenes)
        for refID in overlappingGenes:
            geneDict['overlapping'][refID].append(enhancerString)
            rankDict[refID].append(int(line[-2]))
            superDict[refID].append(int(line[-1]))
            
        overallGeneList.extend(proximalGenes)
        for refID in proximalGenes:
            geneDict['proximal'][refID].append(enhancerString)
            rankDict[refID].append(int(line[-2]))
            superDict[refID].append(int(line[-1]))

    # Create gene-to-enhancer table
    overallGeneList = ROSE_utils.uniquify(overallGeneList)
    rankOrder = ROSE_utils.order([min(rankDict[x]) for x in overallGeneList])
    
    usedNames = []
    for i in rankOrder:
        refID = overallGeneList[i]
        geneName = startDict[refID]['name']
        
        if geneName in usedNames and uniqueGenes:
            continue
            
        usedNames.append(geneName)
        proxEnhancers = geneDict['overlapping'][refID] + geneDict['proximal'][refID]
        superStatus = max(superDict[refID])
        enhancerRanks = ','.join(map(str, rankDict[refID]))
        
        geneToEnhancerTable.append([
            geneName, refID, ','.join(proxEnhancers), enhancerRanks, superStatus
        ])

    # Sort enhancer table by rank if not preserving format
    if not noFormatTable:
        enhancerOrder = ROSE_utils.order([int(line[-2]) for line in enhancerToGeneTable[1:]])
        sortedTable = [enhancerToGeneTable[0]] + [enhancerToGeneTable[i+1] for i in enhancerOrder]
        return sortedTable, geneToEnhancerTable
    
    return enhancerToGeneTable, geneToEnhancerTable


def main():
    """Main function to parse arguments and run the analysis"""
    from optparse import OptionParser
    
    # Set up command line argument parser
    parser = OptionParser(usage="usage: %prog [options] -g GENOME -i INPUT_ENHANCER_FILE")
    
    # Required arguments
    parser.add_option("-i", "--input", dest="input", 
                     help="ROSE ranked enhancer or super-enhancer file")
    parser.add_option("-g", "--genome", dest="genome",
                     help="Genome build (MM10,MM9,MM8,HG18,HG19,HG38)")
    
    # Optional arguments
    parser.add_option("-l", "--list", dest="geneList",
                     help="Gene list file to filter through")
    parser.add_option("-o", "--out", dest="out",
                     help="Output folder (default: same as input file)")
    parser.add_option("-w", "--window", dest="window", type="int", default=50000,
                     help="Search distance for genes (bp, default: 50000)")
    parser.add_option("-f", "--format", dest="formatTable", action="store_true",
                     help="Maintain original formatting of input table")
    
    # Parse arguments
    options, _ = parser.parse_args()
    
    if not options.input or not options.genome:
        parser.print_help()
        sys.exit(1)

    # Set up paths and parameters
    enhancerFile = options.input
    window = options.window
    outFolder = (ROSE_utils.formatFolder(options.out, True) if options.out 
                else os.path.join(os.path.dirname(enhancerFile), ''))
    
    # Get genome annotation file
    genome = options.genome.upper()
    print(f'Using {genome} as the genome')
    
    cwd = os.getcwd()
    genomeDict = {
        'HG18': f'{cwd}/annotation/hg18_refseq.ucsc',
        'MM9': f'{cwd}/annotation/mm9_refseq.ucsc',
        'HG19': f'{cwd}/annotation/hg19_refseq.ucsc',
        'HG38': f'{cwd}/annotation/hg38_refseq.ucsc',
        'MM8': f'{cwd}/annotation/mm8_refseq.ucsc',
        'MM10': f'{cwd}/annotation/mm10_refseq.ucsc',
    }
    
    annotFile = genomeDict.get(genome)
    if not annotFile:
        print(f"Error: Unsupported genome {genome}")
        sys.exit(1)

    # Run analysis
    enhancerToGeneTable, geneToEnhancerTable = mapEnhancerToGene(
        annotFile, enhancerFile, options.geneList, True, window, options.formatTable)
    
    # Write output files
    baseName = os.path.splitext(os.path.basename(enhancerFile))[0]
    suffix = f"_{window//1000}KB" if window != 50000 else ""
    
    out1 = f"{outFolder}{baseName}_ENHANCER_TO_GENE{suffix}.txt"
    ROSE_utils.unParseTable(enhancerToGeneTable, out1, '\t')
    
    out2 = f"{outFolder}{baseName}_GENE_TO_ENHANCER{suffix}.txt"
    ROSE_utils.unParseTable(geneToEnhancerTable, out2, '\t')


if __name__ == "__main__":
    main()

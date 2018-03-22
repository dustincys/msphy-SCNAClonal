#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: getGCMap.py
#          Desc: get GC mappability
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-07-05 08:24:13
#       History:
# =============================================================================
'''
import pysam

import os
from collections import Counter


def getGCMap(bedFileName, outFileFolder, reference_genome):
    """

    :bedFileName: TODO
    :outFileFolder: TODO
    :reference_genome: TODO
    :begwig: TODO
    :returns: TODO

    """
    bedFile = open(bedFileName)
    outFile = open(
        outFileFolder +
        os.path.basename(bedFileName).rstrip(".txt") +
        ".gc.txt",
     'w')

    ref_genome_fasta = pysam.Fastafile(reference_genome)

    outFile.write(
        "chrom\tstart\tend\tsampleReads\treferenceReads\tgc\n")

    for line in bedFile:
        line = line.strip()
        if line == "":
            continue
        listLine = line.split("\t")
        if line.startswith("chrom"):
            continue

        chrom = listLine[0]
        start = int(float(listLine[1]))
        end = int(float(listLine[2]))
        sampleReads = listLine[3]
        referenceReads = listLine[4]

        ref_base = ref_genome_fasta.fetch(chrom, start, end).upper()
        counter = Counter(ref_base)
        N_count = counter['N']
        gc_ratio = (counter['G'] + counter['C']
                    ) * 1.0 / (end - start - N_count + 1)

        outFile.write(
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                chrom,
                start,
                end,
                sampleReads,
                referenceReads,
                gc_ratio))

        print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                chrom,
                start,
                end,
                sampleReads,
                referenceReads,
                gc_ratio)
        pass

    ref_genome_fasta.close()
    outFile.close()
    bedFile.close()

if __name__ == "__main__":
    import optparse
    usage = "usage: <file>"
    parser = optparse.OptionParser(usage)

    options, args = parser.parse_args()

    if len(args) <= 1:
        parser.print_help()
        exit(1)

    getGCMap(args[0], args[1], args[2])

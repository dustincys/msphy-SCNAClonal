#!/usr/bin/env python
'''
# =============================================================================
#      FileName: filterRange.py
#          Desc: filter Range of log(D^T/D^N)
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-12-16 15:34:03
#       History:
# =============================================================================
'''
import argparse
import numpy as np

CHROM_LIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
              '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
              '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']


def main():
    parser = argparse.ArgumentParser(
        description='Converting BICseq segments file to BED file')

    parser.add_argument('inBICseq', help='''Input BICseq file of segments.''')
    parser.add_argument('outBED', help='''Output bed file of segments.''')
    parser.add_argument( '--upper_bound', default=1.0, type=float,
        help='''upper bound for log(D^T/D^N), default 1.0''')
    parser.add_argument( '--lower_bound', default=-1.0, type=float,
        help='''lower bound for log(D^T/D^N), default -1.0''')
    args = parser.parse_args()

    infile = open(args.inBICseq)
    outfile = open(args.outBED, 'w')

    for line in infile:
        if line[0:5] == 'chrom':
            continue

        chrom, start, end, sampleReads, referenceReads, gc\
            = line.strip('\n').split('\t')[0:6]

        temp_y = np.log(int(float(sampleReads)) + 1) - np.log(int(float(referenceReads)) + 1)

        if temp_y < float(args.lower_bound) or \
                temp_y > float(args.upper_bound) or\
                chrom not in CHROM_LIST:
            continue

        outfile.write('\t'.join([chrom, start, end,str(int(float(sampleReads))),
                                 str(int(float(referenceReads))), gc]) + '\n')

    infile.close()
    outfile.close()

if __name__ == '__main__':
    main()

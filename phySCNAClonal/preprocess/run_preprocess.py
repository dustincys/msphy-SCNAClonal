'''
# =============================================================================
#      FileName: run_preprocess.py
#          Desc: run_preprocess
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-11 10:57:17
#       History: Yi Li
# =============================================================================
'''

import sys
from BamToDataConverter import pSCNAClonal_Converter

import time


def run_preprocess_pSCNAClonal(args):
    '''
    args.gc_correction_method: manual, auto
    args.baseline_selection_method: manual, auto
    '''

    print "run preprocess pSCNAClonal"
    print "pkl_path"
    print args.pkl_path
    print "pkl_flag"
    print args.pkl_flag
    time_start = time.time()

    converter = pSCNAClonal_Converter(
        args.normal_bam,
        args.tumor_bam,
        args.reference_genome,
        args.input_filename_base,
        args.segments_bed,
        args.BICseq_bed_corrected,
        args.pkl_path,

        args.max_copynumber,
        args.subclone_num,
        args.baseline_thred_LOH,
        args.baseline_thred_APM,

        min_depth=args.min_depth,
        min_bqual=args.min_base_qual,
        min_mqual=args.min_map_qual,
        process_num=args.process_num
    )

    print "pilflag"
    print args.pkl_flag
    converter.convert(args.gc_correction_method, args.pkl_flag)

    time_end = time.time()

    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    sys.stdout.flush()

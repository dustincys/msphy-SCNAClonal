'''
# =============================================================================
#      FileName: run_preprocess.py
#          Desc: run_preprocess
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-03-11 22:40:41
#       History: Yi Li
# =============================================================================
'''

import sys
import time

from phySCNAClonal.preprocess.converter import BamConverter


def process(args):
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

    converter = BamConverter(
        args.nBamName,
        args.tBamNameL,
        args.refFaName,
        args.pathPreFix,
        args.bedNameL,
        args.bedCorrectedPath,
        args.pklPath,

        args.maxCopyNumber,
        args.subcloneNumL,
        args.baselineThredLOH,
        args.baselineThredAPM,

        minDepth=args.minDepth,
        minBqual=args.minBqual,
        minMqual=args.minMqual,
        processNum=args.processNum
    )

    # print "pilflag"
    # print args.pkl_flag

    converter.convert(readFromBed = True, method = args.gcCorrectionMethod, pklFlag = args.pklFlag)

    time_end = time.time()

    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    sys.stdout.flush()

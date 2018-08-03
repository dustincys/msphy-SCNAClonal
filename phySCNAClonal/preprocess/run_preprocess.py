#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
    print "run preprocess phy-SCNAClonal"
    print "pklPath"
    print args.pklPath
    print "pklFlag"
    print args.pklFlag
    print "minDepth"
    print args.minDepth
    print "processNum"
    print args.processNum
    print "bedCorrectedPath"
    print args.bedCorrectedPath

    time_start = time.time()

    converter = BamConverter(
        args.nBamName,
        args.tBamNameL,
        args.bedNameL,
        args.refFaName,
        args.pathPreFix,
        args.subcloneNumL,
        args.coverageL,
        args.maxCopyNumber,
        args.baselineThredLOH,
        args.baselineThredAPM,
        minDepth=int(args.minDepth),
        minBqual=float(args.minBqual),
        minMqual=float(args.minMqual),
        processNum=int(args.processNum),
        bedCorrectedPath=args.bedCorrectedPath,
        pklPath=args.pklPath)

    # print "pilflag"
    # print args.pkl_flag

    converter.convert(readFromBed = False, method = args.gcCorrectionMethod, pklFlag = True)

    time_end = time.time()

    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    sys.stdout.flush()

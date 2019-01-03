#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: run.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-03-13 14:38:10
#       History:
# =============================================================================
'''

import argparse

from phySCNAClonal.preprocess.run_preprocess import process as run_preprocess
from phySCNAClonal.model.evolve import process as run_modelprocess
from phySCNAClonal.postprocess.write_results import process as run_postprocess


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(
    description=
        'Run phySCNAClonal to infer subclonal composition from SCNA stripes',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

subparsers = parser.add_subparsers()

################
#  preprocess  #
################

parserPreprocess = subparsers.add_parser('preprocess',
                                        help='''Out put preprocess format''')

parserPreprocess.add_argument('--nBamName',
                             help='''BAM file for normal sample.''')

parserPreprocess.add_argument('--tBamNameL', nargs='+',
                             help='''BAM files for tumor samples sorted in \
                               chronological order.''')

parserPreprocess.add_argument('--bedNameL', nargs='+',
                             help='''BED files for segments of each sample in\
                               chronological order.''')

parserPreprocess.add_argument('--refFaName',
                             help='''FASTA file for reference genome.''')

parserPreprocess.add_argument( '--pathPreFix',
                             help='''Base name of the preprocessed input
                             file to be created.''')

parserPreprocess.add_argument('--subcloneNumL', nargs='+', type=int,
                          help='''Set the subclone numbers''')

parserPreprocess.add_argument('--coverageL', nargs='+', type=int,
                          help='''Set the coverage numbers''')

parserPreprocess.add_argument('--maxCopyNumber', default=6, type=int,
                          help='''Set the maximum copy number''')

parserPreprocess.add_argument('--baselineThredLOH', default=0.3, type=float,
                          help='''baseline Thred  of LOH''')

parserPreprocess.add_argument('--baselineThredAPM', default=0.01, type=float,
                          help='''baseline Thred of APM''')

parserPreprocess.add_argument( '--minDepth', default=20, type=int,
                             help='''Minimum reads depth required for both
                             normal and tumor samples.  Default is 20.''')

parserPreprocess.add_argument( '--minBqual', default=10, type=int,
                             help='''Minimum base quality required.
                             Default is 10.''')

parserPreprocess.add_argument( '--minMqual', default=10, type=int,
                             help='''Minimum mapping quality required.
                             Default is 10.''')

parserPreprocess.add_argument( '--processNum', default=1, type=int,
                             help='''Number of processes to launch for
                             preprocessing. Default is 1.''')

parserPreprocess.add_argument('--bedCorrectedPath',
                          help='''The name of corrected BICseq result file''')

parserPreprocess.add_argument('--pklPath',
                          help='''Load the pkl path''')

parserPreprocess.add_argument('--answerFilePath',
                          help='''Load the answer file path''')

parserPreprocess.add_argument('--gcCorrectionMethod', default="auto",
                             help='''The gc correction method, one of auto and
                             visual''')

parserPreprocess.add_argument('--readFromBed', default=False, type=str2bool,
                               help='''get read from Bed (True), from bam file if set it False ''')

parserPreprocess.add_argument('--mergeSeg', default=False, type=str2bool,
                               help='''to merge segment or not to''')

parserPreprocess.add_argument('--pklFlag', default=False, type=str2bool,
                               help='''The pkl flag''')

parserPreprocess.add_argument('--isFixedC', default=False, type=str2bool,
                               help='''Fix Copy number''')

parserPreprocess.set_defaults(func=run_preprocess)

###################
#  model process  #
###################

parserModel = subparsers.add_parser('model',
                                    help='''Output model parameters format''')

parserModel.add_argument( '-b', '--write-backups-every',
                         dest='writeBackupsEvery', default=100, type=int,
                         help='Number of iterations to go between writing\
                         backups of program state')

parserModel.add_argument( '-S', '--write-state-every', dest='writeStateEvery',
                         default=10, type=int, help= 'Number of iterations\
                         between writing program state to disk. Higher values\
                         reduce IO burden at the cost of losing progress made\
                         if program is interrupted.')

parserModel.add_argument('-k', '--top-k-trees', dest='topKTrees',
                         default='topKTrees', help='Output file to save top-k\
                         trees in text format')

parserModel.add_argument('-f', '--clonal-freqs', dest='clonalFreqs',
                         default='clonalFrequencies', help='Output file to save\
                         clonal frequencies')

parserModel.add_argument('-B', '--burnin-samples',
                         dest='burninSampleNum', default=1000, type=int,
                         help='Number of burnin samples')

parserModel.add_argument('-s', '--mcmc-samples', dest='mcmcSampleNum',
                         default=2500, type=int, help='Number of MCMC samples')

parserModel.add_argument('-i', '--mh-iterations', dest='mhIterations',
                         default=5000, type=int, help='Number of\
                         Metropolis-Hastings iterations')

parserModel.add_argument('-r', '--random-seed', dest='randomSeed',
                         type=int, help='Random seed for initializing MCMC\
                         sampler')

parserModel.add_argument('-t', '--tmp-dir', dest='tmpDir', help='Path to\
                         directory for temporary files')


parserModel.add_argument('-p', '--params', dest='paramsFile', help='JSON\
                         file listing run parameters, generated by the parser')

parserModel.add_argument('--inputDataFile', dest='inputDataFile', help= 'File\
                         listing data(SCNA data, either semgent or stripe). For\
                         proper format, see README.md.')

parserModel.add_argument('--inputDataTextFile', dest='inputDataTextFile', help= 'Text\
                         file listing data(SCNA stripes). For proper format,\
                         see README.md.')

parserModel.add_argument('--isMerged', default=False, type=str2bool,
                               help='''is merged data file''')

parserModel.add_argument('--isCrossing', default=False, type=str2bool,
                               help='''using crossing file''')

parserModel.add_argument('--crossingFile', dest='crossingFile', help= 'The crossing file.')

parserModel.add_argument('--isSingleCell', default=False, type=str2bool,
                               help='''using single cell file''')

parserModel.add_argument('--singleCellFile', dest='singleCellFile', help= 'The single cell file.')

parserModel.add_argument('--maxCopyNumber', dest='maxCopyNumber', default=6,
                         type=int, help= 'Max copy number')

parserModel.add_argument('--noTag', default=False, type=str2bool,
                         help='''to remove all tag or not''')

parserModel.set_defaults(func=run_modelprocess)


#################
#  postprocess  #
#################


parserPostprocess = subparsers.add_parser('postprocess',
                                help='''Output postprocess parameters format''')

parserPostprocess.add_argument('--treeFile', help='File containing sampled trees')

parserPostprocess.add_argument('--SCNAPoolFile', help='File containing SCNA pool')

parserPostprocess.add_argument('--answerFilePath', default="", help='Answer file path')

parserPostprocess.add_argument('--outputFolder', default="", help='Output folder')

# parserPostprocess.add_argument('mutlistOutput', help='Output file for\
                               # JSON-formatted list of mutations')
#
# parserPostprocess.add_argument('mutassOutput', help= 'Output file for\
#                                JSON-formatted list of SSMs and CNVs assigned to\
#                                each subclone')

parserPostprocess.set_defaults(func=run_postprocess)

args = parser.parse_args()
args.func(args)

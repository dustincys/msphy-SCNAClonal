'''
# =============================================================================
#      FileName: data.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-04-24 15:30:17
#       History: YI lI
# =============================================================================
'''
import sys
import numpy as np
import scipy.cluster.hierarchy as hcluster
from scipy.stats.mstats import gmean
from collections import Counter

import util2 as u


# from GCBASELINE import constants
# from GCBASELINE.preprocess.utils import *
from utils import get_chrom_format, get_chrom_lens_idxs, BEDnParser,\
    chrom_name_to_idx, chrom_idx_to_name, get_segment_name, BEDParser,\
    get_LOH_frac, get_APM_frac_MAXMIN, get_LOH_status, get_APM_status


from utils import get_loga, get_cn_allele_config, get_mu_E_joint,\
    log_binomial_likelihood, mad_based_outlier


import constants


class Segment:

    def __init__(self):
        self.name = ""
        self.chromIdx = -1
        self.chromName = ""
        self.start = -1
        self.end = -1
        self.nReadNum = -1
        self.tReadNum = -1
        self.gc = -1

        self.LOHFrac = -1
        self.LOHStatus = 'NONE'
        self.APMFrac = -1
        self.APMStatus = 'NONE'

        self.pairedCounts = None
        self.BAFCounts = None

        self.baselineLabel = 'FALSE'
        self.tag = 'SCNA'
        self.stripeIdx = -1
        self.alleleType = 'NONE'

        self.copyNumber = -1

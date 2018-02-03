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
        self.chrom_idx = -1
        self.chrom_name = ""
        self.start = -1
        self.end = -1
        self.normal_reads_num = -1
        self.tumor_reads_num = -1
        # self.sites_num = 0
        self.LOH_frac = -1
        self.LOH_status = 'NONE'
        self.APM_frac = -1
        self.APM_status = 'NONE'

        # 此处应该设置一下类别
        # "baseline", "1", "2", "3", "4"
        # default 1, all the segs belong to the SCNA happens on time 1
        # self.baseline_label = 'FALSE'
        self.tag = '1'

        self.log2_ratio = 0.0

        self.paired_counts = None
        self.BAF_counts = None
        self.copy_number = -1
        self.stripe_number = -1
        self.allele_type = 'NONE'
        self.subclone_prev = -1
        self.subclone_cluster = 'NONE'

        self.gc = -1
        self.log_ratio = 0

        # save the likelihood of the last phi, {phi: (likelihood, copy_number,
        # pi)}
        self.phi_last = None

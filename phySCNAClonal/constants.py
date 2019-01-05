'''
# =============================================================================
#      FileName: constants.py
#          Desc: constants/parameters for preprocess
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-03-03 08:52:28
#       History:
# =============================================================================
'''

import numpy as np


###################
#  RD parameters  #
###################

MINIMUM_POSITIVE=0.0000001

###########################
#  MCMC model parameters  #
###########################

ZOOM_P = 5
X_ZOOM_IN_FACTOR = 10

DOWN_GC_BOUNDARY_PERCENTILE = 30
UP_GC_BOUNDARY_PERCENTILE = 70
DOWN_LOGA_BOUNDARY_PERCENTILE = 30
UP_LOGA_BOUNDARY_PERCENTILE = 70
SLOPE_RANGE = 5


##############
#  COVERAGE  #
##############

COVERAGE = 300

####################
#  BAF parameters  #
####################

COVERAGE_PROPERATION_BAF = 0.8
BAF_THRESHOLD = 0.1
BAF_N_MIN = 0.35

# BAF_COUNTS_MIN = 10
# BAF_COUNTS_MAX = 95

BAF_BINS = np.array(range(0, 100 + 1))/100.0
LOH_FRAC_MAX = 0.35
SITES_NUM_MIN = 5
BINOM_TEST_P = 0.5
# default BINOM_TEST_THRED = 0.025
# for 30x read depth,
# BINOM_TEST_THRED = 0.025
BINOM_TEST_THRED = 0.002
BINOM_TEST_THRED_APM = 0.18

# Them
APM_N_MIN = 0.43 # This parameter is very important for baseline selection

EMPIRI_BAF = 0.5
EMPIRI_AAF = 1.0 - EMPIRI_BAF
MU_N = EMPIRI_BAF/(EMPIRI_BAF + EMPIRI_AAF)

########################################
#  Hierarchical clustering parameters  #
########################################

# HCC1954  n5t95 ~ n60t40
HC_THRESH = 0.005
# HCC1954  n80t20
#HC_THRESH = 0.001
# HCC1954  n95t5
#HC_THRESH = 0.0001


########################
#  Copy number config  #
########################

COPY_NUMBER_NORMAL = 2


######################
#  Chrom parameters  #
######################

CHROM_START = 0

CHROM_IDX_LIST = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
              20, 21, 22]


################################
#  Model structure parameters  #
################################


RD_WEIGHT = 0.5
RD_WEIGHT_TSSB = 0.5

VARPI = 0.8 ## gap parameters

#############################
#  Stripe range parameters  #
#############################

YDOWNL = [-0.8, -0.75, -5, -5]
YUPL = [0.6, 0.5, 5, 5]
STRIPENUML = [10, 15, 12, 12]
NOISESTRIPENUML = [0, 0, 1, 1]

#################################
#  Stripe decompose parameters  #
#################################

DECOMPOSE_NUMBER_THRESHOLD = 10

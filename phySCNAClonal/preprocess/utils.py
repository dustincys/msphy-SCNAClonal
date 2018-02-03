'''
Created on 2013-07-27

@author: Yi Li

pyloh.preprocess.utils

================================================================================

Modified on 2014-04-09

@author: Yi Li
'''

import sys
from multiprocessing import Pool

import numpy as np
from scipy.misc import comb
from scipy.special import gammaln
from scipy.stats import beta, binom

import constants

def show(*string):
    print string
    sys.stdout.flush()

def BEDnParser(bed_file_name):
    """TODO: Docstring for BEDnParser.
    :returns: TODO

    """
    inbed = open(bed_file_name)

    chroms = []
    starts = []
    ends = []
    tumorReadsn = []
    normalReads = []
    gcs = []

    for line in inbed:
        fields = line.split('\t')
        chrom_name = fields[0]
        chrom_idx = chrom_name_to_idx(chrom_name)

        if chrom_idx == -1:
            continue

        chrom_name, start, end, tumorrd, normalrd, gc = fields[0:6]

        chroms.append(chrom_name)
        starts.append(int(start))
        ends.append(int(end))
        tumorReadsn.append(int(tumorrd))
        normalReads.append(int(normalrd))
        gcs.append(float(gc))

    inbed.close()

    return (chroms, starts, ends, tumorReadsn, normalReads, gcs)


def BEDParser(bed_file_name):
    inbed = open(bed_file_name)

    chroms = []
    starts = []
    ends = []
    gcs = []

    for line in inbed:
        fields = line.split('\t')
        chrom_name = fields[0]
        chrom_idx = chrom_name_to_idx(chrom_name)

        if chrom_idx == -1:
            continue

        chrom_name, start, end, gc = fields[0:4]

        chroms.append(chrom_name)
        starts.append(int(start))
        ends.append(int(end))
        gcs.append(float(gc))

    inbed.close()

    return (chroms, starts, ends, gcs)


def chrom_idx_to_name(idx, format):
    if format == 'UCSC':
        chrom_name = 'chr' + str(idx)
    elif format == 'ENSEMBL':
        chrom_name = str(idx)
    else:
        print 'Error: %s not supported' % (format)
        sys.exit(1)

    return chrom_name


def chrom_name_to_idx(chrom_name):
    idx = -1

    try:
        idx = int(chrom_name.strip('chr'))
    except:
        pass

    return idx


def get_chrom_format(chroms):
    format = 'NONE'

    for chrom in chroms:
        if chrom[0:3] == 'chr':
            format = 'UCSC'
            break
        else:
            try:
                int(chrom)
                format = 'ENSEMBL'
                break
            except:
                pass

    if format == 'NONE':
        print 'Error: %s not supported' % (chrom)
        sys.exit(-1)
    else:
        return format


def get_chrom_lens_idxs(chrom_idx_list, sam_SQ):
    chrom_lens = []
    chrom_idxs = []
    for i in range(0, len(chrom_idx_list)):
        chrom_idx = chrom_idx_list[i]

        for j in range(0, len(sam_SQ)):
            if chrom_idx == chrom_name_to_idx(sam_SQ[j]['SN']):
                chrom_lens.append(int(sam_SQ[j]['LN']))
                chrom_idxs.append(chrom_idx)
                break

    return (chrom_lens, chrom_idxs)


def get_segment_name(chrom_name, start, end):

    return '_'.join([chrom_name, 'start', str(start), 'end', str(end)])


def normal_heterozygous_filter(counts):
    BAF_N_MIN = constants.BAF_N_MIN

    I = counts.shape[0]
    idx_keep = []

    for i in xrange(0, I):
        a_N = counts[i, 0] * 1.0
        b_N = counts[i, 1] * 1.0
        d_N = a_N + b_N
        BAF_N = b_N / d_N

        if BAF_N >= BAF_N_MIN and BAF_N <= 0.5:
            idx_keep.append(i)

    counts = counts[idx_keep]

    return counts


def get_BAF_counts(counts):
    BAF_bins = constants.BAF_BINS

    a_N = counts[:, 0] * 1.0
    b_N = counts[:, 1] * 1.0
    a_T = counts[:, 2] * 1.0
    b_T = counts[:, 3] * 1.0

    BAF_N = b_N / (a_N + b_N)
    BAF_T = b_T / (a_T + b_T)

    BAF_counts, _, _ = np.histogram2d(BAF_N, BAF_T, bins=(BAF_bins, BAF_bins))

    return BAF_counts


def get_APM_frac_MAXMIN_SNP(counts):
    """get the baf position that are average in the tumor bam

    :counts: TODO
    :returns: TODO

    """

    I = counts.shape[0]

    sites_num_min = 1

    APM_N_MIN = constants.APM_N_MIN

    if I < sites_num_min:
        APM_frac = -1
        return APM_frac

    a_T = counts[:, 2]
    b_T = counts[:, 3]
    d_T = a_T + b_T
    l_T = np.min(counts[:, 2:4], axis=1)
    p_T = l_T * 1.0 / d_T

    APM_num = np.where(np.logical_and(p_T > APM_N_MIN, p_T <= 0.5))[0].shape[0]
    APM_frac = APM_num * 1.0 / (I + 1.0)

    return APM_frac
    pass


def get_APM_frac_MAXMIN(counts):
    """get the baf position that are average in the tumor bam

    :counts: TODO
    :returns: TODO

    """

    I = counts.shape[0]

    sites_num_min = constants.SITES_NUM_MIN

    APM_N_MIN = constants.APM_N_MIN

    if I < sites_num_min:
        APM_frac = -1

        return APM_frac

    a_T = counts[:, 2]
    b_T = counts[:, 3]
    d_T = a_T + b_T
    l_T = np.min(counts[:, 2:4], axis=1)
    p_T = l_T * 1.0 / d_T

    APM_num = np.where(np.logical_and(p_T > APM_N_MIN, p_T <= 0.5))[0].shape[0]
    APM_frac = APM_num * 1.0 / I

    return APM_frac

    pass


def get_LOH_frac_SNP(counts):
    I = counts.shape[0]

    sites_num_min = 1
    p = constants.BINOM_TEST_P
    thred = constants.BINOM_TEST_THRED

    if I < sites_num_min:
        LOH_frac = -1
        return LOH_frac

    a_T = counts[:, 2]
    b_T = counts[:, 3]
    d_T = a_T + b_T
    l_T = np.min(counts[:, 2:4], axis=1)
    p_T = binom.cdf(l_T, d_T, p)

    LOH_num = np.where(p_T < thred)[0].shape[0]
    LOH_frac = (LOH_num + 1.0) * 1.0 / (I + 1.0)

    return LOH_frac


def get_LOH_frac(counts):
    I = counts.shape[0]

    sites_num_min = constants.SITES_NUM_MIN
    p = constants.BINOM_TEST_P
    thred = constants.BINOM_TEST_THRED

    if I < sites_num_min:
        LOH_frac = -1

        return LOH_frac

    a_T = counts[:, 2]
    b_T = counts[:, 3]
    d_T = a_T + b_T
    l_T = np.min(counts[:, 2:4], axis=1)
    p_T = binom.cdf(l_T, d_T, p)

    LOH_num = np.where(p_T < thred)[0].shape[0]
    LOH_frac = LOH_num * 1.0 / I

    return LOH_frac


def get_APM_frac(counts):
    """get the baf position that are average in the tumor bam

    :counts: TODO
    :returns: TODO

    """

    I = counts.shape[0]

    sites_num_min = constants.SITES_NUM_MIN
    p = constants.BINOM_TEST_P
    thred = constants.BINOM_TEST_THRED_APM

    if I < sites_num_min:
        APM_frac = -1

        return APM_frac

    a_T = counts[:, 2]
    b_T = counts[:, 3]
    d_T = a_T + b_T
    l_T = np.min(counts[:, 2:4], axis=1)
    p_T = binom.cdf(l_T, d_T, p)

    APM_num = np.where(p_T > thred)[0].shape[0]
    APM_frac = APM_num * 1.0 / I

    return APM_frac

    pass


def get_LOH_status(LOH_frac, baseline_thred):
    LOH_FRAC_MAX = constants.LOH_FRAC_MAX

    if LOH_frac < 0:
        LOH_status = 'NONE'
    elif LOH_frac < baseline_thred:
        LOH_status = 'FALSE'
    elif LOH_frac >= baseline_thred and LOH_frac < LOH_FRAC_MAX:
        LOH_status = 'UNCERTAIN'
    elif LOH_frac >= LOH_FRAC_MAX:
        LOH_status = 'TRUE'
    else:
        LOH_status = 'ERROR'

    return LOH_status


def get_APM_status(APM_frac, baseline_thred_APM):
    if APM_frac < 0:
        APM_status = "NONE"
    elif APM_frac > baseline_thred_APM:
        APM_status = "TRUE"
    else:
        APM_status = "FALSE"

    return APM_status


def remove_outliers(X):
    std_thred = 0.05

    idx_keep = []

    n = X.shape[0]

    for i in range(0, n):
        if np.abs(X[i] - X.mean()) <= X.std():
            idx_keep.append(i)

    if len(idx_keep) == 0 or len(idx_keep) == n:
        return X

    X = X[idx_keep]

    if X.std() < std_thred:
        return X
    else:
        return remove_outliers(X)


def calculate_BAF(tumorData, normalData, chrmsToUse, minSNP, gamma,
                  process_num):

    # function to select columns from a 2D list
    select_col = lambda array, colNum: map(lambda x: x[colNum], array)

    # vectors of tumor data
    tumorMutCount = select_col(tumorData, 3)
    tumorRefCount = select_col(tumorData, 2)

    # vectors of normal data
    normalMutCount = select_col(normalData, 3)
    normalRefCount = select_col(normalData, 2)

    # denominators for BAFs
    tumorDenom = map(sum, zip(tumorMutCount, tumorRefCount))
    normalDenom = map(sum, zip(normalMutCount, normalRefCount))

    tumorBAF = []
    normalBAF = []
    newTumorData = []
    newNormalData = []
    print "Determining heterozygosity."
    p = Pool(process_num)
    repGamma = [gamma for i in range(len(tumorData))]
    isHet = p.map(is_heterozygous,
                  zip(normalRefCount, normalMutCount, repGamma))
    print "Calculating BAFs."
    for i in range(len(tumorData)):
        chrm = tumorData[i][0]
        # filter out data that uses irrelevant chromosomes
        if chrm not in chrmsToUse:
            continue
        # filter out data where there aren't enough tumor reads
        if tumorMutCount[i] + tumorRefCount[i] < minSNP:
            continue
        # filter out data where there aren't enough normal reads
        if normalMutCount[i] + normalRefCount[i] < minSNP:
            continue

        currTumorNum = tumorMutCount[i]
        currTumorDenom = tumorDenom[i]

        currNormalNum = normalMutCount[i]
        currNormalDenom = normalDenom[i]

        # filter out points with denominators that are 0 to prevent zero
        # division error
        if currTumorDenom == 0 or currNormalDenom == 0:
            continue
        else:
            tumorBAFj = currTumorNum / currTumorDenom
            normalBAFj = currNormalNum / currNormalDenom
            # filter out data where normal BAFs do not fit in bounds correctly
            if isHet[i]:
                tumorBAF.append(tumorBAFj)
                normalBAF.append(normalBAFj)
                newTumorData.append(tumorData[i])
                newNormalData.append(normalData[i])
            else:
                continue

    return tumorBAF, normalBAF, newTumorData, newNormalData


def filter_normal_heterozygous(tumorData, normalData, gamma, process_num):
    # function to select columns from a 2D list
    select_col = lambda array, colNum: map(lambda x: x[colNum], array)

    # vectors of tumor data
    tumorMutCount = select_col(tumorData, 3)
    tumorRefCount = select_col(tumorData, 2)

    # vectors of normal data
    normalMutCount = select_col(normalData, 3)
    normalRefCount = select_col(normalData, 2)

    # denominators for BAFs
    tumorDenom = map(sum, zip(tumorMutCount, tumorRefCount))
    normalDenom = map(sum, zip(normalMutCount, normalRefCount))

    tumorBAF = []
    normalBAF = []
    newTumorData = []
    newNormalData = []
    print "Determining heterozygosity."
    repGamma = [gamma for i in range(len(tumorData))]
    wholeData = zip(normalRefCount, normalMutCount, repGamma)
    p = Pool(process_num)
    print wholeData[0]
    isHet = p.map(is_heterozygous, wholeData)

    tumorData_filtered = select_col(
        filter(lambda x: x[1], zip(tumorData, isHet)), 0)
    normalData_filtered = select_col(
        filter(lambda x: x[1], zip(normalData, isHet)), 0)

    return tumorData_filtered, normalData_filtered


def is_heterozygous(xxx_todo_changeme):
    """
    Determines if an allele should be considered heterozygous.

    Arguments:
            n_a (int): number of a alleles counted. Used as alpha parameter for the beta distribution
            n_b (int): number of b alleles counted. Used as beta parameter for the beta distribution
            gamma (float): parameter used for deciding heterozygosity; determined via a beta distribution
                                            with 1 - gamma confidence

    Returns:
            A boolean indicating whether or not the allele should be considered heterozygous.
    """
    (n_a, n_b, gamma) = xxx_todo_changeme
    if n_a == -1 or n_b == -1:
        return False
    if n_a == 0 or n_b == 0:
        return False

    p_lower = gamma / 2.0
    p_upper = 1 - p_lower

    [c_lower, c_upper] = beta.ppf([p_lower, p_upper], n_a + 1, n_b + 1)
    return c_lower <= 0.5 and c_upper >= 0.5


def read_snp_file(filename):
    """
    Converts an SNP file to a 2D array.

    Args:
            filename: location of SNP file
    Returns:
            data: data from SNP file in a 2D array, where each row is [chromosome, position, ref count, mut count]
    """

    print "Reading SNP file at " + filename

    data = []

    f = gzip.open(filename, "r") if ".gz" in filename else open(filename, "r")
    splitChar = "," if ".csv" in filename else "\t"

    chrmInd = 0
    posInd = 1

    for line in f:
        if line.strip() == "":
            continue
        if line.startswith("#"):
            continue

        if len(line.split(splitChar)) < 8:
            refInd = 2
            mutInd = 3
        else:
            refInd = 7
            mutInd = 8

        vals = line.split(splitChar)

        chrm = vals[chrmInd].lower()

        if chrm.startswith("chrm"):
            chrm = chrm[4:]
        if chrm.startswith("chr"):
            chrm = chrm[3:]

        if chrm == "x":
            chrm = 23
        elif chrm == "y":
            chrm = 24
        else:
            chrm = int(chrm)

        position = int(vals[posInd])
        refCount = float(vals[refInd])
        mutCount = float(vals[mutInd])
        data.append([chrm, position, refCount, mutCount])

    return data


def get_row_by_segment(tumorData, normalData, segment):
    """get the tumor and normal snp in this segment

    :tumorData: TODO
    :normalData: TODO
    :segment: TODO
    :returns: TODO

    """
    tumorData_temp = filter(
        lambda item: item[0] == int(segment.chrom_name) and (item[1] >= segment.start and item[1] <= segment.end),
        tumorData)
    normalData_temp = filter(
        lambda item: item[0] == int(segment.chrom_name) and (item[1] >= segment.start and item[1] <= segment.end),
        normalData)

    return tumorData_temp, normalData_temp


def get_paired_counts(tumorData, normalData):
    """

    :tumorData: TODO
    :normalData: TODO
    :returns: TODO

    """

    paired_counts_temp = []
    for i in range(len(normalData)):
        paired_counts_temp.append([
            int(normalData[i][2]),
            int(normalData[i][3]),
            int(tumorData[i][2]),
            int(tumorData[i][3]),
            int(normalData[i][0]),
            int(normalData[i][1])
        ])
    paired_counts_j = np.array(paired_counts_temp)

    return paired_counts_j


def get_loga(data):
    return np.log(data.tumor_reads_num + 1) - np.log(data.normal_reads_num + 1)


def log_poisson_likelihood(k, Lambda):
    return k * np.log(Lambda) - Lambda - gammaln(k + 1)


def get_cn_allele_config(max_copynumber):
    cn_allele_config = {}
    for cn in range(0, max_copynumber + 1):
        allele_config = {}
        for M_num in range(0, (cn + 2) / 2):
            P_num = cn - M_num
            if P_num == 0 and M_num == 0:
                mu_T = constants.EMPIRI_BAF / (
                    constants.EMPIRI_AAF + constants.EMPIRI_BAF)
                pi_T = 'NULL'
            elif P_num == M_num:
                mu_T = 0.5
                pi_T = 'P' * P_num + 'M' * M_num
            else:
                mu_T = (M_num * 1.0) / cn
                pi_T = 'P' * P_num + 'M' * M_num + '/' + 'P' * M_num + 'M' * P_num
            allele_config[pi_T] = mu_T
        cn_allele_config[cn] = allele_config
        # {2:{PP/MM:0, PM:0.5},     3:{...},...}
    return cn_allele_config


def get_mu_E_joint(mu_N, mu_G, c_N, c_H, phi):
    return ((1.0 - phi) * c_N * mu_N + phi * c_H * mu_G) / (
        (1.0 - phi) * c_N + phi * c_H)


def log_binomial_likelihood(k, n, mu):
    # example:
        # k: array([2, 3, 2])
        # n: array([2, 3, 2])
        # mu: 0.3
    # return array([2, 3, 2])
    return np.log(comb(n, k)) + k * np.log(mu) + (n - k) * np.log(1 - mu)


def mad_based_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score > thresh

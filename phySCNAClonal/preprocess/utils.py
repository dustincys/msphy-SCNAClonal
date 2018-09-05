'''
# =============================================================================
#      FileName: utils.py
#          Desc: Functions for preprocess
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-03-03 08:55:34
#       History: Yi Li
# =============================================================================
'''

import sys
from multiprocessing import Pool

import numpy as np
from scipy.misc import comb
from scipy.special import gammaln
from scipy.stats import beta, binom

import phySCNAClonal.constants as constants

class AnswerIndex(object):

    """Break points for segments. for searching index"""

    def __init__(self, answerFilePath):
        """init with answer file """
        self._answerFilePath = answerFilePath

        self._colNames = ""
        self._answerDict = {}
        # {1_12123_123123: "xxxx",}

        self._startBreakpoints = {}
        self._endBreakpoints = {}
        # {1:[12,23,1,23],}
        self.__index()

    def valueName(self):
        return self._colNames

    def __index(self):
        with open(self._answerFilePath) as answerFile:
            totalLine = 0
            for line in answerFile:
                totalLine = totalLine + 1
                line = line.strip()
                if line == "":
                    continue
                if 1 == totalLine:
                    self._colNames = line
                    continue
                listLine = line.split("\t")

                answer = "\t".join(listLine[3:])
                key = "_".join(listLine[0:3])
                self._answerDict[key] = answer

                chrom = listLine[0]
                start = int(listLine[1])
                end = int(listLine[2])

                if chrom not in self._startBreakpoints.keys():
                    self._startBreakpoints[chrom] = []
                if chrom not in self._endBreakpoints.keys():
                    self._endBreakpoints[chrom] = []

                self._startBreakpoints[chrom].append(start)
                self._endBreakpoints[chrom].append(end)

        for chrom in self._startBreakpoints.keys():
            self._startBreakpoints[chrom] = \
                np.array(self._startBreakpoints[chrom])
            self._endBreakpoints[chrom] = \
                np.array(self._endBreakpoints[chrom])
            self._startBreakpoints[chrom].sort()
            self._endBreakpoints[chrom].sort()

    def getValue(self, chrom, start, end):
        chrom = str(chrom)
        start = int(start)
        end = int(end)
        head = np.max(self._startBreakpoints[chrom][self._startBreakpoints[chrom] <=
                                             start])
        tail = np.min(self._endBreakpoints[chrom][self._endBreakpoints[chrom] >=
                                                  end])

        headIdx = np.where(head == self._startBreakpoints[chrom])[0][0]
        tailIdx = np.where(tail == self._endBreakpoints[chrom])[0][0]

        key = ""
        if 0 == tailIdx - headIdx:
            key = "{0}_{1}_{2}".format(chrom, str(head), str(tail))
        elif 1 == tailIdx - headIdx:
            leftDist = np.abs(head - start)
            rightDist = np.abs(tail - end)
            if leftDist < rightDist:
                key = "{0}_{1}_{2}".format(chrom, str(head),
                                           str(self._endBreakpoints[chrom][tailIdx-1]))
            else:
                key = "{0}_{1}_{2}".format(chrom, str(self._startBreakpoints[chrom][headIdx+1]),
                                           str(tail))
        else:
            key = "{0}_{1}_{2}".format(chrom, str(self._startBreakpoints[chrom][headIdx+1]),
                                           str(self._endBreakpoints[chrom][headIdx+2]))

        return self._answerDict[key]



def BEDnParser(bedName):
    """TODO: Docstring for BEDnParser.
    :returns: TODO

    """
    inBed = open(bedName)

    chromNameL = []
    startL = []
    endL = []
    tReadNumL = []
    nReadNumL = []
    gcL = []

    for line in inBed:
        fields = line.split('\t')
        chromName = fields[0]
        chromIdx = chrom_name_to_idx(chromName)

        if chromIdx == -1:
            continue

        chromName, start, end, tReadNum, nReadNum, gc = fields[0:6]

        chromNameL.append(chromName)
        startL.append(int(start))
        endL.append(int(end))
        tReadNumL.append(int(tReadNum))
        nReadNumL.append(int(nReadNum))
        gcL.append(float(gc))

    inBed.close()

    return (chromNameL, startL, endL, tReadNumL, nReadNumL, gcL)


def BEDParser(bedName):
    inBed = open(bedName)

    chromNameL = []
    startL = []
    endL = []
    gcL = []

    for line in inBed:
        fields = line.split('\t')
        chromName = fields[0]
        chromIdx = chrom_name_to_idx(chromName)

        if chromIdx == -1:
            continue

        chromName, start, end = fields[0:3]
        gc = fields[-1]

        chromNameL.append(chromName)
        startL.append(int(start))
        endL.append(int(end))
        gcL.append(float(gc))

    inBed.close()

    return (chromNameL, startL, endL, gcL)


def chrom_idx_to_name(idx, format):
    if format == 'UCSC':
        chromName = 'chr' + str(idx)
    elif format == 'ENSEMBL':
        chromName = str(idx)
    else:
        print 'Error: %s not supported' % (format)
        sys.exit(1)

    return chromName


def chrom_name_to_idx(chromName):
    idx = -1

    try:
        idx = int(chromName.strip('chr'))
    except:
        pass

    return idx


def get_chrom_format(chromNameL):
    format = 'NONE'

    for chrom in chromNameL:
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


def get_segment_name(chromName, start, end):

    return '_'.join([chromName, 'start', str(start), 'end', str(end)])


def normal_heterozygous_filter(counts):
    BAFNMin = constants.BAF_N_MIN

    I = counts.shape[0]
    idxKeep = []

    for i in xrange(0, I):
        aN = counts[i, 0] * 1.0
        bN = counts[i, 1] * 1.0
        dN = aN + bN
        BAFN = bN / dN

        if BAFN >= BAFNMin and BAFN <= 0.5:
            idxKeep.append(i)

    counts = counts[idxKeep]

    return counts


def get_BAF_counts(counts):
    BAFbins = constants.BAF_BINS

    aN = counts[:, 0] * 1.0
    bN = counts[:, 1] * 1.0
    aT = counts[:, 2] * 1.0
    bT = counts[:, 3] * 1.0

    BAFN = bN / (aN + bN)
    BAFT = bT / (aT + bT)

    BAFCounts, _, _ = np.histogram2d(BAFN, BAFT, bins=(BAFbins, BAFbins))

    return BAFCounts


def get_APM_frac_MAXMIN_SNP(counts):
    """get the baf position that are average in the tumor bam

    :counts: TODO
    :returns: TODO

    """

    I = counts.shape[0]

    sitesNumMin = 1

    APMNMin = constants.APM_N_MIN

    if I < sitesNumMin:
        APMFrac = -1
        return APMFrac

    aT = counts[:, 2]
    bT = counts[:, 3]
    dT = aT + bT
    lT = np.min(counts[:, 2:4], axis=1)
    pT = lT * 1.0 / dT

    APMNum = np.where(np.logical_and(pT > APMNMin, pT <= 0.5))[0].shape[0]
    APMFrac = APMNum * 1.0 / (I + 1.0)

    return APMFrac
    pass


def get_APM_frac_MAXMIN(counts):
    """get the baf position that are average in the tumor bam

    :counts: TODO
    :returns: TODO

    """

    I = counts.shape[0]

    sitesNumMin = constants.SITES_NUM_MIN

    APMNMin = constants.APM_N_MIN

    if I < sitesNumMin:
        APMFrac = -1

        return APMFrac

    aT = counts[:, 2]
    bT = counts[:, 3]
    dT = aT + bT
    lT = np.min(counts[:, 2:4], axis=1)
    pT = lT * 1.0 / dT

    APMNum = np.where(np.logical_and(pT > APMNMin, pT <= 0.5))[0].shape[0]
    APMFrac = APMNum * 1.0 / I

    return APMFrac

    pass


def get_LOH_frac_SNP(counts):
    I = counts.shape[0]

    sitesNumMin = 1
    p = constants.BINOM_TEST_P
    thred = constants.BINOM_TEST_THRED

    if I < sitesNumMin:
        LOHFrac = -1
        return LOHFrac

    aT = counts[:, 2]
    bT = counts[:, 3]
    dT = aT + bT
    lT = np.min(counts[:, 2:4], axis=1)
    pT = binom.cdf(lT, dT, p)

    LOHNum = np.where(pT < thred)[0].shape[0]
    LOHFrac = (LOHNum + 1.0) * 1.0 / (I + 1.0)

    return LOHFrac


def get_LOH_frac(counts):
    I = counts.shape[0]

    sitesNumMin = constants.SITES_NUM_MIN
    p = constants.BINOM_TEST_P
    thred = constants.BINOM_TEST_THRED

    if I < sitesNumMin:
        LOHFrac = -1

        return LOHFrac

    aT = counts[:, 2]
    bT = counts[:, 3]
    dT = aT + bT
    lT = np.min(counts[:, 2:4], axis=1)
    pT = binom.cdf(lT, dT, p)

    LOHNum = np.where(pT < thred)[0].shape[0]
    LOHFrac = LOHNum * 1.0 / I

    return LOHFrac


def get_APM_frac(counts):
    """get the baf position that are average in the tumor bam

    :counts: TODO
    :returns: TODO

    """

    I = counts.shape[0]

    sitesNumMin = constants.SITES_NUM_MIN
    p = constants.BINOM_TEST_P
    thred = constants.BINOM_TEST_THRED_APM

    if I < sitesNumMin:
        APMFrac = -1

        return APMFrac

    aT = counts[:, 2]
    bT = counts[:, 3]
    dT = aT + bT
    lT = np.min(counts[:, 2:4], axis=1)
    pT = binom.cdf(lT, dT, p)

    APMNum = np.where(pT > thred)[0].shape[0]
    APMFrac = APMNum * 1.0 / I

    return APMFrac

    pass


def get_LOH_status(LOHFrac, baselineThred):
    LOHFracMax = constants.LOH_FRAC_MAX

    if LOHFrac < 0:
        LOHStatus = 'NONE'
    elif LOHFrac < baselineThred:
        LOHStatus = 'FALSE'
    elif LOHFrac >= baselineThred and LOHFrac < LOHFracMax:
        LOHStatus = 'UNCERTAIN'
    elif LOHFrac >= LOHFracMax:
        LOHStatus = 'TRUE'
    else:
        LOHStatus = 'ERROR'

    return LOHStatus


def get_APM_status(APMFrac, baselineThredAPM):
    if APMFrac < 0:
        APMStatus = "NONE"
    elif APMFrac > baselineThredAPM:
        APMStatus = "TRUE"
    else:
        APMStatus = "FALSE"

    return APMStatus


def remove_outliers(X):
    stdThred = 0.05

    idxKeep = []

    n = X.shape[0]

    for i in range(0, n):
        if np.abs(X[i] - X.mean()) <= X.std():
            idxKeep.append(i)

    if len(idxKeep) == 0 or len(idxKeep) == n:
        return X

    X = X[idxKeep]

    if X.std() < stdThred:
        return X
    else:
        return remove_outliers(X)


def calculate_BAF(tumorData, normalData, chrmsToUse, minSNP, gamma,
                  processNum):

    # function to select columns from a 2D list
    selectCol = lambda array, colNum: map(lambda x: x[colNum], array)

    # vectors of tumor data
    tumorMutCount = selectCol(tumorData, 3)
    tumorRefCount = selectCol(tumorData, 2)

    # vectors of normal data
    normalMutCount = selectCol(normalData, 3)
    normalRefCount = selectCol(normalData, 2)

    # denominators for BAFs
    tumorDenom = map(sum, zip(tumorMutCount, tumorRefCount))
    normalDenom = map(sum, zip(normalMutCount, normalRefCount))

    tumorBAF = []
    normalBAF = []
    newTumorData = []
    newNormalData = []
    print "Determining heterozygosity."
    p = Pool(processNum)
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


def filter_normal_heterozygous(tumorData, normalData, gamma, processNum):
    # function to select columns from a 2D list
    selectCol = lambda array, colNum: map(lambda x: x[colNum], array)

    # vectors of tumor data
    tumorMutCount = selectCol(tumorData, 3)
    tumorRefCount = selectCol(tumorData, 2)

    # vectors of normal data
    normalMutCount = selectCol(normalData, 3)
    normalRefCount = selectCol(normalData, 2)

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
    p = Pool(processNum)
    print wholeData[0]
    isHet = p.map(is_heterozygous, wholeData)

    tumorData_filtered = selectCol(
        filter(lambda x: x[1], zip(tumorData, isHet)), 0)
    normalData_filtered = selectCol(
        filter(lambda x: x[1], zip(normalData, isHet)), 0)

    return tumorData_filtered, normalData_filtered


def is_heterozygous(xxxTodoChangeme):
    """
    Determines if an allele should be considered heterozygous.

    Arguments:
            nA (int): number of a alleles counted. Used as alpha parameter for the beta distribution
            nB (int): number of b alleles counted. Used as beta parameter for the beta distribution
            gamma (float): parameter used for deciding heterozygosity; determined via a beta distribution
                                            with 1 - gamma confidence

    Returns:
            A boolean indicating whether or not the allele should be considered heterozygous.
    """
    (nA, nB, gamma) = xxxTodoChangeme
    if nA == -1 or nB == -1:
        return False
    if nA == 0 or nB == 0:
        return False

    pLower = gamma / 2.0
    pUpper = 1 - pLower

    [cLower, cUpper] = beta.ppf([pLower, pUpper], nA + 1, nB + 1)
    return cLower <= 0.5 and cUpper >= 0.5


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
    tumorDataTemp = filter(
        lambda item: item[0] == int(segment.chromName) and (item[1] >= segment.start and item[1] <= segment.end),
        tumorData)
    normalDataTemp = filter(
        lambda item: item[0] == int(segment.chromName) and (item[1] >= segment.start and item[1] <= segment.end),
        normalData)

    return tumorDataTemp, normalDataTemp


def get_paired_counts(tumorData, normalData):
    """

    :tumorData: TODO
    :normalData: TODO
    :returns: TODO

    """

    pairedCountsTemp = []
    for i in range(len(normalData)):
        pairedCountsTemp.append([
            int(normalData[i][2]),
            int(normalData[i][3]),
            int(tumorData[i][2]),
            int(tumorData[i][3]),
            int(normalData[i][0]),
            int(normalData[i][1])
        ])
    pairedCountsJ = np.array(pairedCountsTemp)

    return pairedCountsJ


def get_loga(data):
    return np.log(data.tReadNum + 1) - np.log(data.nReadNum + 1)


def log_poisson_likelihood(k, Lambda):
    return k * np.log(Lambda) - Lambda - gammaln(k + 1)


def get_cn_allele_config(maxCopyNumber):
    cnAlleleConfig = {}
    for cn in range(0, maxCopyNumber + 1):
        alleleConfig = {}
        for MNum in range(0, (cn + 2) / 2):
            PNum = cn - MNum
            if PNum == 0 and MNum == 0:
                muT = constants.EMPIRI_BAF / (
                    constants.EMPIRI_AAF + constants.EMPIRI_BAF)
                piT = 'NULL'
            elif PNum == MNum:
                muT = 0.5
                piT = 'P' * PNum + 'M' * MNum
            else:
                muT = (MNum * 1.0) / cn
                piT = 'P' * PNum + 'M' * MNum + '/' + 'P' * MNum + 'M' * PNum
            alleleConfig[piT] = muT
        cnAlleleConfig[cn] = alleleConfig
        # {2:{PP/MM:0, PM:0.5},     3:{...},...}
    return cnAlleleConfig


def get_mu_E_joint(muN, muG, cN, cH, phi):
    return ((1.0 - phi) * cN * muN + phi * cH * muG) / (
        (1.0 - phi) * cN + phi * cH)


def log_binomial_likelihood(k, n, mu):
    # example:
        # k: array([2, 3, 2])
        # n: array([2, 3, 2])
        # mu: [0.3, 0.2]
    # return array([2, 3, 2])
    # print k.shape, n.shape, mu.shape
    nn = n * np.ones((mu.shape[0], n.shape[0]))
    kk = k * np.ones((mu.shape[0], k.shape[0]))
    mumu = mu[np.newaxis,:].T * np.ones((mu.shape[0], n.shape[0]))
    ll = np.log(comb(nn, kk)) + kk * np.log(mumu) + (nn - kk) * np.log(1 - mumu)
    return ll.transpose()


def mad_based_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    if med_abs_deviation == 0:
        modified_z_score = [-float("Inf")]
    else:
        modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score > thresh

def getBAFofSeg(seg):
    aT = seg.pairedCounts[:, 2]
    bT = seg.pairedCounts[:, 3]
    dT = aT + bT
    lT = np.min(seg.pairedCounts[:, 2:4], axis=1)
    pT = lT * 1.0 / dT
    outlier = mad_based_outlier(pT)
    return np.median(pT[np.invert(outlier)])

def main():
    ansIdx = AnswerIndex("/media/f/simdata/rd30/perfectSegmentation/rd30stage1.seg.gc.txt.answer")
    print ansIdx.valueName()


if __name__ == "__main__":
    main()

'''
# =============================================================================
#      FileName: converter.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-01-31 14:07:22
#       History:
# =============================================================================
'''

import sys
import pickle as pkl
from multiprocessing import Pool
import numpy as np
import pysam

from mcmc import MCMCLM

from phySCNAClonal.preprocess.data import SegmentPool
from phySCNAClonal.preprocess.stripe import Stripe, DataStripes
from phySCNAClonal.preprocess.iofun import PairedCountsIterator, PairedPileupIterator

from phySCNAClonal.preprocess.utils import get_BAF_counts, normal_heterozygous_filter

from plotGC import GCStripePlot


class BamConverter:

    def __init__(self, nBamName, tBamNameL, bedNameL, refFaName, pathPrefix="",
                 subcloneNumberL=[2], coverageL = [30], maxCopyNumber=6,
                 baselineThredLOH=0.3, baselineThredAPM=0.01, minDepth=20,
                 minBqual=10, minMqual=10, processNum=1, bedCorrectedPath="",
                 pklPath=""):
        self._nBamName = nBamName
        self._tBamNameL = tBamNameL
        self._bedNameL = bedNameL
        self._refFaName = refFaName

        self.__pathPrefix = pathPrefix

        self.__subcloneNumberL = subcloneNumberL
        self.__coverageL = coverageL
        self.__maxCopyNumber = maxCopyNumber
        self.__baselineThredLOH = baselineThredLOH
        self.__baselineThredAPM = baselineThredAPM

        self.__minDepth = minDepth
        self.__minBqual = minBqual
        self.__minMqual = minMqual

        self.__processNum = processNum
        self.__bedCorrectedPath=bedCorrectedPath
        self.__pklPath = pklPath

        self._segPoolL = []

    def convert(self, readFromBed=True, method="auto", pkl_flag=False):
        self._load_segs(readFromBed)
        self._correct_bias(method)
        blSegsL = self._get_baseline()
        self._mark_timestamp(blSegsL)
        stripePool = self._generate_stripe()
        self._dump(stripePool)

    def _dump(self, stripePool):
        fileName = self.__pathPrefix + self.__pklPath
        outFile = open(fileName, 'wb')
        pkl.dump(stripePool, outFile, protocol=2)
        outFile.close()

    def _generate_stripe(self):
        """
        generate stripe from segs
        """
        # 此处应该对每一个tag进行单独的聚类
        # 也有可能无法确定有多少个条带
        #
        # 另一种方案应该是对所有的条带进行聚类，然后从中挑选出目标条带，分解为新
        # 条带。

        yDown = constants.YDOWNL[self._segPoolL[-1].idx]
        yUp = constants.YUPL[self._segPoolL[-1].idx]
        # 此处应该近似为最大拷贝数与亚克隆数的乘积，作为手工输入也可以
        stripeNum = constants.STRIPENUML[self._segPoolL[-1].idx]
        noiseStripeNum = constants.NOISESTRIPENUML[self._segPoolL[-1].idx]

        tempSP = StripePool(self._segPoolL[-1], self._segPoolL[-1].baseline,
                            yDown, yUp, stripeNum, noiseStripeNum)
        tempSP.get(byTag = True)

    def _mark_timestamp(self, blSegsL):
        """
        mark segs in final sample
        """
        # 此处应用R来进行求解

        # 首先，求解每相邻数据的基线之差的集合
        #
        # 或直接列出所有基线

        # 然后，根据相邻数据的基线之差，映射到数据的非基线之上，确定归宿于哪一个
        # 基线之差
        #
        # 或找出落入基线之中的最大索引

        # 最后，所有的数据点中最先落入基线之差的为目标时间戳
        #
        # 根据该索引作为时间戳

        from rpy2.robjects.packages import importr
        from rpy2.robjects import IntVector, StrVector, globalenv
        import rpy2.robjects as robjects

        GR = importr('GenomicRanges')
        IR = importr('IRanges')

        GRL = GR.GRangesList()
        globalenv["GRL"] = GRL
        for blSegs, idx in zip(blSegsL, range(len(blSegsL))):
            chromNames = StrVector([seg.chromName for seg in blSegs])
            starts = IntVector([seg.start for seg in blSegs])
            ends = IntVector([seg.end for seg in blSegs])
            tempGR = GR.GRanges(seqnames = chromNames, ranges=IR.IRanges(starts, ends))
            globalenv["tempGR"] = tempGR
            GRL = robjects.r("GRL[[{0}]]=tempGR".format(idx))

        # 此处由于list中保存的是指向目标Seg的指针，所以更新nonBLSegs即可
        nonBlSegs = list(set(self._segPoolL[-1].segments) - set(blSegsL[-1]))
        chromNames = StrVector([seg.chromName for seg in nonBlSegs])
        starts = IntVector([seg.start for seg in nonBlSegs])
        ends = IntVector([seg.end for seg in nonBlSegs])
        nonBlGR = GR.GRanges(seqnames = chromNames, ranges=IR.IRanges(starts, ends))
        fo = IR.findOverlaps(GR2, GRL)
        globalenv["fo"] = fo
        robjects.reval("fom <- as.matrix(fo)")
        overlapIdx = np.array(list(robjects.r.fom)).reshape(tuple(reversed(robjects.r.fom.dim))) - 1
        # [[2, 2, 3, 3],
        # [1, 2, 1, 2]]
        #
        for index in set(overlapIdx[0,]):
            yIdx = np.where(overlapIdx[0,]==index)[0]
            ts = np.max(re[1,yIdx])
            # tag = 0, 1, 2, 3, 4, ..., BASELINE
            nonBlGR[index].tag = str(ts)

    def _load_segs(self, readFromBed=True):
        """
        load segments for each tumor sample
        """
        assert len(self._tBamNameL) == len(self._bedNameL)
        assert len(self._tBamNameL) == len(self.__subcloneNumberL)

        for tBamName, bedName, coverage, subcloneNumber in zip(self._tBamNameL,
            self._bedNameL, self.__coverageL, self.__subcloneNumberL):
            print >> sys.stdout, 'Loading segments from bam file:\n{0}\n'.format(tBamName)
            print >> sys.stdout, 'and bed file with gc:\n{0}\n'.format(bedName)
            tempSP = SegmentPool(self.__maxCopyNumber, coverage)
            if not readFromBed:
                nBam = pysam.Samfile(self._nBamName, 'rb')
                tBam = pysam.Samfile(tBamName, 'rb')
                tempSP.load_seg_bam(nBam, tBam, bedName)
                nBam.close()
                tBam.close()
            else:
                tempSP.load_seg_bed (self.segments_bed)
            self._segPoolL.append(tempSP)

    def _correct_bias(self, method="auto"):
        """
        correct bias of each tumor sample
        """
        assert len(self._segPoolL) == len(self.__subcloneNumberL)
        for segPool, subcloneNumber in zip(self._segPoolL,
                                               self.__subcloneNumberL):
            if "auto" == method:
                self._MCMC_GC_C(SegmentPool, subcloneNumber)
            elif "visual" == method:
                self._V_GC_C(SegmentPool, len(SegmentPool.segments))

    def _get_baseline(self):
        """
        get the baseline segments
        calculate baseline of each SegmentPool
        return: the baseline segments list of each SegmentPool
        """

        blSegsL = []

        for segPool, idx in zip(self._segPoolL, range(len(self._segPoolL))):
            tempBL = segPool.get_baseline(self,maxCopyNumber,
                                              self.subcloneNumberL[idx])
            blSegsL.append(tempBL)

        return blSegsL

    def _MCMC_GC_C(self, segPool, subcloneNumber):
        """
        The interception is irrelevant for correction, set as median
        MCMCLM only returns the m and c, then correct the segPool here
        """

        mcmclm = MCMCLM(segPool, 0, subcloneNumber, self.__maxCopyNumber)
        m, c = mcmclm.run()
        print "MCMC slope = {}".format(m)

        x = np.array(map(lambda seg: seg.gc, segPool.segments))
        y = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                         np.log(seg.normal_reads_num + 1),
                         segPool.segments))

        y_corrected = self._correct(x, y, m, c)

        for i in range(len(y_corrected)):
            segPool.segments[i].tumor_reads_num = np.exp(
                y_corrected[i] +
                np.log(segPool.segments[i].normal_reads_num + 1)
            ) - 1
            segPool.segments[i].log_ratio = np.log(
                (y_corrected[i] + 1.0) /
                (segPool.segments[i].normal_reads_num + 1.0)
            )

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(slope, intercept)

    def _correct(self, x, y, slope, intercept):
        K = np.percentile(y, 50)
        A = slope * x + intercept
        return y - A + K

    def visualize(self):
        gsp = GCStripePlot(self.segPool.segments, len(self.segPool.segments))
        print "total number: {}".format(self.segPool.segNum)
        gsp.plot()
        x, y, m, c = gsp.output()
        print "x, y, m, c"
        print x, y, m, c

    def _V_GC_C(self, segPool, sampleNumber=10000):
        gsp = GCStripePlot(segPool.segments, sampleNumber)
        print >> sys.stdout, "total number: {}".format(len(segPool.segments))
        gsp.plot()
        print >> sys.stdout, "x, y, m, c"
        print >> sys.stdout, gsp.output()

        x = np.array(map(lambda seg: seg.gc, segPool.segments))
        y = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                         np.log(seg.normal_reads_num + 1), segPool.segments))
        y_corrected = self._correct(x, y, m, c)

        for i in range(len(y_corrected)):
            segPool.segments[i].tumor_reads_num = np.exp(
                y_corrected[i] +
                np.log(segPool.segments[i].normal_reads_num + 1)
            ) - 1
            segPool.segments[i].log_ratio = np.log(
                (y_corrected[i] + 1.0) /
                (segPool.segments[i].normal_reads_num + 1.0)
            )

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(slope, intercept)

    def _baseline_selection(self):
        print "begin baseline selection.."
        self._get_LOH_frac()
        self._get_LOH_status()
        self._get_APM_frac()
        self._get_APM_status()
        self._compute_Lambda_S()

    def _get_APM_status(self):
        self.segPool.get_APM_status(self.__baselineThredAPM)

    def _get_LOH_status(self):
        self.segPool.get_LOH_status(self.__baselineThredLOH,
                                 flag_runpreprocess=True)

    def _compute_Lambda_S(self):
        print "begin compute lambda s .."
        self.segPool.compute_Lambda_S_LOH(self.__maxCopyNumber, self.__subcloneNumberL,
                                       flag_runpreprocess=True)

    def _get_counts(self, tBamName, segPool):
        """
        get allele counts of target bam file
        save the counts into segPool
        """

        segNum = len(self.segPool.segments)
        processNum = self.__processNum
        print "processNum = {}".format(processNum)

        if processNum > segNum:
            processNum = segNum

        pool = Pool(processes=processNum)

        argsL = []

        for j in range(0, segNum):
            segName = self.segPool.segments[j].name
            chromName = self.segPool.segments[j].chromName
            chromIdx = self.segPool.segments[j].chromIdx
            start = self.segPool.segments[j].start
            end = self.segPool.segments[j].end

            argsT = (
                segName,
                chromName,
                chromIdx,
                start,
                end,
                self._nBamName,
                tBamName,
                self._refFaName,
                self.__minDepth,
                self.__minBqual,
                self.__minMqual)

            argsL.append(argsT)

        countsTL = pool.map(process_by_segment, argsL)

        for j in range(0, segNum):
            pairedCounts_j, BAFCounts_j = countsTL[j]
            segPool.segments[j].pairedCounts = pairedCounts_j
            segPool.segments[j].BAFCounts = BAFCounts_j

    def _get_LOH_frac(self):
        self.segPool.get_LOH_frac()

    def _get_APM_frac(self):
        self.segPool.get_APM_frac()
# ===============================================================================
#  Function
# ===============================================================================


def process_by_segment(argsT):
    seg_name, chrom_name, chrom_idx, start, end, nBamName,\
        tBamName, refFaName, minDepth, minBqual,\
        minMqual = argsT

    print 'Preprocessing segment {0}...'.format(seg_name)
    sys.stdout.flush()

    nBam = pysam.Samfile(nBamName, 'rb')
    tBam = pysam.Samfile(tBamName, 'rb')
    ref_genome_fasta = pysam.Fastafile(refFaName)

    normal_pileup_iter = nBam.pileup(chrom_name, start, end)
    tumor_pileup_iter = tBam.pileup(chrom_name, start, end)

    paired_pileup_iter = PairedPileupIterator(
        normal_pileup_iter, tumor_pileup_iter, start, end)
    paired_counts_iter = PairedCountsIterator(
        paired_pileup_iter,
        ref_genome_fasta,
        chrom_name,
        chrom_idx,
        minDepth,
        minBqual,
        minMqual)

    paired_counts_j, BAF_counts_j = iterator_to_counts(paired_counts_iter)
    counts_tuple_j = (paired_counts_j, BAF_counts_j)

    nBam.close()
    tBam.close()
    ref_genome_fasta.close()

    return counts_tuple_j


def iterator_to_counts(paired_counts_iter):
    buff = 100000

    paired_counts_j = np.array([[], [], [], [], [], []], dtype=int).transpose()
    BAF_counts_j = np.zeros((100, 100))
    buff_counts = []
    i = 0

    for counts in paired_counts_iter:
        buff_counts.append(counts)
        i = i + 1

        if i < buff:
            continue

        buff_counts = np.array(buff_counts)

        if buff_counts.shape[0] != 0:
            BAF_counts_buff = get_BAF_counts(buff_counts)
            BAF_counts_j += BAF_counts_buff

        buff_counts_filtered = normal_heterozygous_filter(buff_counts)

        if buff_counts_filtered.shape[0] != 0:
            paired_counts_j = np.vstack((paired_counts_j, buff_counts_filtered))

        buff_counts = []
        i = 0

    buff_counts = np.array(buff_counts)

    if buff_counts.shape[0] != 0:
        BAF_counts_buff = get_BAF_counts(buff_counts)
        BAF_counts_j += BAF_counts_buff

    buff_counts_filtered = normal_heterozygous_filter(buff_counts)

    if buff_counts_filtered.shape[0] != 0:
        paired_counts_j = np.vstack((paired_counts_j, buff_counts_filtered))

    return (paired_counts_j, BAF_counts_j)

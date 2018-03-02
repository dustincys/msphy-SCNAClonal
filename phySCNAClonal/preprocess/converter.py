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

from phySCNAClonal.preprocess.data import Data
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

    def convert(self, method, pkl_flag=False):
        self._load_segs()
        self._correct_bias()
        blSegsL = self._get_baseline()
        self._mark_timestamp(blSegsL)
        self._generate_stripe()
        self._mark_stripe()
        self._dump()

    def _generate_stripe(self):
        """
        generate stripe from segs
        """
        # generate the stripe according to the tag




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

    def _load_segs(self):
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
            nBam = pysam.Samfile(self._nBamName, 'rb')
            tBam = pysam.Samfile(tBamName, 'rb')
            tempSP.load_segments(nBam, tBam, bedName)
            nBam.close()
            tBam.close()
            self._segPoolL.append(tempSP)

    def _correct_bias(self, method=""):
        """
        correct bias of each tumor sample
        """
        assert len(self._segPoolL) == len(self.__subcloneNumberL)
        for segmentPool, subcloneNumber in zip(self._segPoolL,
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

        for segmentPool, idx in zip(self._segPoolL, range(len(self._segPoolL))):
            tempBL = segmentPool.get_baseline(self,maxCopyNumber,
                                              self.subcloneNumberL[idx])
            blSegsL.append(tempBL)

        return blSegsL

    def _MCMC_GC_C(self, data, subcloneNumber):
        """
        The interception is irrelevant for correction, set as median
        MCMCLM only returns the m and c, then correct the data here
        """

        mcmclm = MCMCLM(data, 0, subcloneNumber, self.__maxCopyNumber)
        m, c = mcmclm.run()
        print "MCMC slope = {}".format(m)

        x = np.array(map(lambda seg: seg.gc, data.segments))
        y = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                         np.log(seg.normal_reads_num + 1),
                         data.segments))

        y_corrected = self._correct(x, y, m, c)

        for i in range(len(y_corrected)):
            data.segments[i].tumor_reads_num = np.exp(
                y_corrected[i] +
                np.log(data.segments[i].normal_reads_num + 1)
            ) - 1
            data.segments[i].log_ratio = np.log(
                (y_corrected[i] + 1.0) /
                (data.segments[i].normal_reads_num + 1.0)
            )

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(slope, intercept)

    def _correct(self, x, y, slope, intercept):
        K = np.percentile(y, 50)
        A = slope * x + intercept
        return y - A + K

    def visualize(self):
        gsp = GCStripePlot(self.data.segments, len(self.data.segments))
        print "total number: {}".format(self.data.seg_num)
        gsp.plot()
        x, y, m, c = gsp.output()
        print "x, y, m, c"
        print x, y, m, c

    def _V_GC_C(self, segmentPool, sampleNumber=10000):
        gsp = GCStripePlot(segmentPool.segments, sampleNumber)
        print >> sys.stdout, "total number: {}".format(len(segmentPool.segments))
        gsp.plot()
        print >> sys.stdout, "x, y, m, c"
        print >> sys.stdout, gsp.output()

        x = np.array(map(lambda seg: seg.gc, data.segments))
        y = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                         np.log(seg.normal_reads_num + 1), data.segments))
        y_corrected = self._correct(x, y, m, c)

        for i in range(len(y_corrected)):
            data.segments[i].tumor_reads_num = np.exp(
                y_corrected[i] +
                np.log(data.segments[i].normal_reads_num + 1)
            ) - 1
            data.segments[i].log_ratio = np.log(
                (y_corrected[i] + 1.0) /
                (data.segments[i].normal_reads_num + 1.0)
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
        self.data.get_APM_status(self.__baselineThredAPM)

    def _get_LOH_status(self):
        self.data.get_LOH_status(self.__baselineThredLOH,
                                 flag_runpreprocess=True)

    def _compute_Lambda_S(self):
        print "begin compute lambda s .."
        self.data.compute_Lambda_S_LOH(self.__maxCopyNumber, self.__subcloneNumberL,
                                       flag_runpreprocess=True)

    def _load_segmentsn(self):
        """
        :returns: TODO

        """
        nBam = pysam.Samfile(self._nBamName, 'rb')
        tBam = pysam.Samfile(self.tBam_filename, 'rb')

        print 'Loading normalized segments by {0}...'.format(self.segments_bed)
        sys.stdout.flush()
        self.data.load_segmentsn(nBam, tBam, self.segments_bed)

        nBam.close()
        tBam.close()


    def _load_segments_bed(self, segments_bed, data):
        print 'Loading segments with gc by {0}...'.format(self.segments_bed)
        sys.stdout.flush()
        data.load_segments_bed(segments_bed)

    def _get_counts(self, tBam_filename, data):
        seg_num = self.data.seg_num
        processNum = self.__processNum
        print "processNum = {}".format(processNum)

        if processNum > seg_num:
            processNum = seg_num

        pool = Pool(processes=processNum)

        args_list = []

        for j in range(0, seg_num):
            seg_name = self.data.segments[j].name
            chrom_name = self.data.segments[j].chrom_name
            chrom_idx = self.data.segments[j].chrom_idx
            start = self.data.segments[j].start
            end = self.data.segments[j].end

            args_tuple = (
                seg_name,
                chrom_name,
                chrom_idx,
                start,
                end,
                self._nBamName,
                tBam_filename,
                self._refFaName,
                self.__minDepth,
                self.__minBqual,
                self.__minMqual)

            args_list.append(args_tuple)

        counts_tuple_list = pool.map(process_by_segment, args_list)

        for j in range(0, seg_num):
            paired_counts_j, BAF_counts_j = counts_tuple_list[j]

            data.segments[j].paired_counts = paired_counts_j
            data.segments[j].BAF_counts = BAF_counts_j

    def _get_LOH_frac(self):
        self.data.get_LOH_frac()

    def _get_APM_frac(self):
        self.data.get_APM_frac()
# ===============================================================================
#  Function
# ===============================================================================


def process_by_segment(args_tuple):
    seg_name, chrom_name, chrom_idx, start, end, nBamName,\
        tBam_filename, refFaName, minDepth, minBqual,\
        minMqual = args_tuple

    print 'Preprocessing segment {0}...'.format(seg_name)
    sys.stdout.flush()

    nBam = pysam.Samfile(nBamName, 'rb')
    tBam = pysam.Samfile(tBam_filename, 'rb')
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

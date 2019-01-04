#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

import pickle as pkl
import sys
from multiprocessing import Pool

import numpy as np
import pysam

from copy import deepcopy

import phySCNAClonal.constants as constants
from phySCNAClonal.preprocess.data.pools.segmentPool import SegmentPool
from phySCNAClonal.preprocess.data.pools.stripePool import StripePool
from phySCNAClonal.preprocess.iofun import (PairedCountsIterator,
                                            PairedPileupIterator)
from phySCNAClonal.preprocess.mcmc import MCMCLM
from phySCNAClonal.preprocess.plotGC import GCStripePlot, GCStripePoolPlot
from phySCNAClonal.preprocess.utils import (get_BAF_counts,AnswerIndex,
                                            dump_seg_to_txt,
                                            dump_seg_to_txt_with_answer,
                                            dump_seg_to_txt_list,
                                            normal_heterozygous_filter,
                                            updateFixedCValue)

np.set_printoptions(threshold=np.inf)




class BamConverter:

    def __init__(self,
                 nBamName,
                 tBamNameL,
                 bedNameL,
                 refFaName,
                 pathPrefix="",
                 subcloneNumberL=[2],
                 coverageL=[30],
                 maxCopyNumber=6,
                 baselineThredLOH=0.3,
                 baselineThredAPM=0.01,
                 minDepth=20,
                 minBqual=10,
                 minMqual=10,
                 processNum=1,
                 bedCorrectedPath="",
                 pklPath="",
                 answerFilePath="",
                 isFixedC=False):
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
        self.__answerFilePath = answerFilePath
        self.__isFixedC = isFixedC

        self._segPoolL = []

    def convert(self, readFromBed=True, method="auto", mergeSeg=False,
                pklFlag=False):
        if not pklFlag:
            self._load_segs(readFromBed)
            # self._correct_bias(method)
            self._load_allele_counts()
            self._dump(self._segPoolL, ".temp.segPoolL")

            pklFile = open(self.__pklPath, 'wb')
            pkl.dump(self._segPoolL, pklFile, protocol=2)
            pklFile.close()
        else:
            pklFile = open(self.__pklPath, 'rb')
            self._segPoolL = pkl.load(pklFile )
            pklFile .close()


        # self._correct_bias(method)
        # self._dump(self._segPoolL, ".temp.gccorrected.segPoolL")
        blSegsL = self._get_baseline(mergeSeg)
        self._mark_timestamp(blSegsL)
        debug = True
        if debug:
            for blSegL in blSegsL:
                logaL = [np.log(seg.tReadNum + 1) - np.log(seg.nReadNum + 1) for
                         seg in blSegL]
                print "blSegL: min(logaL) = {}".format(min(logaL))
                print "blSegL: max(logaL) = {}".format(max(logaL))

            for segPool in self._segPoolL:
                logaL = [np.log(seg.tReadNum + 1) - np.log(seg.nReadNum + 1) for
                         seg in segPool.segments if seg.tag == "BASELINE"]
                print "tag: min(logaL) = {}".format(min(logaL))
                print "tag: max(logaL) = {}".format(max(logaL))

        if self.__isFixedC:
            # only update the segment pool it seems there is no need to update
            # the fixed C value for stripe pool
            self._updateFixedCValue()

        if mergeSeg:
            stripePool = self._generate_stripe()
            self._dump(stripePool, "stripePool.pkl")
            self._dump_txt(stripePool, "stripePool.txt")
        else:
            segmentPool = self._generate_segment()
            self._dump(segmentPool, "lastSegPoolNoBlSegs.pkl")
            self._dump_txt(segmentPool, "lastSegPoolNoBlSegs.txt")

        self._dump(self._segPoolL, "allSegPoolL.pkl")
        self._dump_seg_to_txt()

        if debug:
            for blSegL in blSegsL:
                logaL = [np.log(seg.tReadNum + 1) - np.log(seg.nReadNum + 1) for
                         seg in blSegL]
                print "min(logaL) = {}".format(min(logaL))
                print "max(logaL) = {}".format(max(logaL))
            for segPool in self._segPoolL:
                logaL = [np.log(seg.tReadNum + 1) - np.log(seg.nReadNum + 1) for
                         seg in segPool.segments if seg.tag == "BASELINE"]
                print "tag: min(logaL) = {}".format(min(logaL))
                print "tag: max(logaL) = {}".format(max(logaL))

    def _updateFixedCValue(self):
        updateFixedCValue(self._segPoolL[-1], self.__answerFilePath)

    def _dump_seg_to_txt(self):
        """
        output table for R, draw figures
        """
        if self.__answerFilePath != "" and not self.__answerFilePath is None:
            dump_seg_to_txt_with_answer(self._segPoolL[-1],
                                        len(self._segPoolL) - 1,
                                        self.__answerFilePath,
                                        self.__pathPrefix)
        else:
            # dump_seg_to_txt(self._segPoolL[-1], len(self._segPoolL) - 1,
                            # self.__pathPrefix)
            # dump_seg_to_txt(self._segPoolL[0], len(self._segPoolL) - 2,
                            # self.__pathPrefix)
            dump_seg_to_txt_list(self._segPoolL, self.__pathPrefix)

    def _dump_txt(self, pool, outFilePath):
        """
        out put pool into plain txt
        """
        pool.output_txt(self.__pathPrefix + "/" + outFilePath)

    def _load_allele_counts(self):
        for tBamName, segPool in zip(self._tBamNameL, self._segPoolL):
            self._get_counts(tBamName, segPool)

    def _dump(self, data, dpFileName):
        fileName = self.__pathPrefix + "/" + dpFileName
        outFile = open(fileName, 'wb')
        pkl.dump(data, outFile, protocol=2)
        outFile.close()

    def _generate_segment(self):
        # or no need to copy it?
        lastSegPool = deepcopy(self._segPoolL[-1])
        lastSegPool.segments = filter(lambda item:item.baselineLabel != "TRUE",
                                      lastSegPool.segments)
        return lastSegPool

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

        # 如上一步中添加了baseline选项
        # 那么这里的排序需要先排除baseline选项
        # gspp = GCStripePoolPlot(tempSP)
        # gspp.plot()

        tempSP.stripes.sort(key = lambda item: int(item.tag))

        for idx, sp in enumerate(tempSP.stripes):
            tempSP.stripes[idx].id = idx

        return tempSP

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
            robjects.r("GRL[[{0}]]=tempGR".format(str(idx+1)))
            GRL = robjects.r["GRL"]

        # 此处由于list中保存的是指向目标Seg的指针，所以更新nonBLSegs即可
        nonBlSegs = list(set(self._segPoolL[-1].segments) - set(blSegsL[-1]))
        chromNames = StrVector([seg.chromName for seg in nonBlSegs])
        starts = IntVector([seg.start for seg in nonBlSegs])
        ends = IntVector([seg.end for seg in nonBlSegs])
        nonBlGR = GR.GRanges(seqnames = chromNames, ranges=IR.IRanges(starts, ends))

        fo = IR.findOverlaps(nonBlGR, GRL)
        globalenv["fo"] = fo
        robjects.reval("fom <- as.matrix(fo)")
        overlapIdx = np.array(list(robjects.r.fom)).reshape(tuple(reversed(robjects.r.fom.dim))) - 1
        # [[2, 2, 3, 3],
        # [1, 2, 1, 2]]
        #
        print overlapIdx

        for index in set(overlapIdx[0,]):
            yIdxes = np.where(overlapIdx[0,]==index)[0]
            ts = np.max(overlapIdx[1,yIdxes]+1)
            nonBlSegs[index].tag = str(ts)

    def _load_segs(self, readFromBed=True):
        """
        load segments for each tumor sample
        """
        assert len(self._tBamNameL) == len(self._bedNameL)
        assert len(self._tBamNameL) == len(self.__subcloneNumberL)

        for tBamName, bedName, coverage, subcloneNumber, i in zip(self._tBamNameL,
            self._bedNameL, self.__coverageL, self.__subcloneNumberL, range(len(self._tBamNameL))):
            print >> sys.stdout, 'Loading segments from bam file:\n{0}\n'.format(tBamName)
            print >> sys.stdout, 'and bed file with gc:\n{0}\n'.format(bedName)
            tempSP = SegmentPool(i, self.__maxCopyNumber, coverage)
            if not readFromBed:
                nBam = pysam.Samfile(self._nBamName, 'rb')
                tBam = pysam.Samfile(tBamName, 'rb')
                tempSP.load_seg_bam(nBam, tBam, bedName)
                nBam.close()
                tBam.close()
            else:
                tempSP.load_seg_bed (bedName)
            self._segPoolL.append(tempSP)

    def _correct_bias(self, method="auto"):
        """
        correct bias of each tumor sample
        """
        assert len(self._segPoolL) == len(self.__subcloneNumberL)
        for segPool, subcloneNumber in zip(self._segPoolL,
                                               self.__subcloneNumberL):
            if "auto" == method:
                self._MCMC_GC_C(segPool, subcloneNumber)
            elif "visual" == method:
                self._V_GC_C(segPool, len(segPool.segments))

    def _get_baseline(self, mergeSeg=False):
        """
        get the baseline segments
        calculate baseline of each SegmentPool
        return: the baseline segments list of each SegmentPool
        """

        blSegsL = []

        for segPool, idx in zip(self._segPoolL, range(len(self._segPoolL))):
            tempBL = segPool.get_baseline(self.__maxCopyNumber ,
                                          self.__subcloneNumberL[idx],
                                          self.__baselineThredLOH,
                                          self.__baselineThredAPM,
                                          mergeSeg,
                                          isPreprocess=True,
                                          index=idx)
            print len(tempBL)
            blSegsL.append(tempBL)
            # self.visualize(segPool)

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
        y = np.array(map(lambda seg: np.log(seg.tReadNum + 1) -
                         np.log(seg.nReadNum + 1), segPool.segments))

        yCorrected = self._correct(x, y, m, c)

        for i in range(len(yCorrected)):
            segPool.segments[i].tReadNum = np.exp(
                yCorrected[i] +
                np.log(segPool.segments[i].nReadNum + 1)
            ) - 1
            segPool.segments[i].log_ratio = np.log(
                (yCorrected[i] + 1.0) /
                (segPool.segments[i].nReadNum + 1.0)
            )

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(m, c)

    def _correct(self, x, y, slope, intercept):
        K = np.percentile(y, 50)
        A = slope * x + intercept
        return y - A + K

    def visualize(self, segPool):
        gsp = GCStripePlot(segPool.segments, len(segPool.segments))
        print "total number: {}".format(len(segPool.segments))
        gsp.plot()
        _, _, m, c = gsp.output()
        print "m, c"
        print m, c

    def _V_GC_C(self, segPool, sampleNumber=10000):
        gsp = GCStripePlot(segPool.segments, 10000)
        print >> sys.stdout, "total number: {}".format(len(segPool.segments))
        gsp.plot()
        # print >> sys.stdout, "x, y, m, c"
        # print >> sys.stdout, gsp.output()

        _, _, m, c = gsp.output()
        print >> sys.stdout, m, c

        x = np.array(map(lambda seg: seg.gc, segPool.segments))
        y = np.array(map(lambda seg: np.log(seg.tReadNum + 1) -
                         np.log(seg.nReadNum + 1), segPool.segments))
        yCorrected = self._correct(x, y, m, c)

        for i in range(len(yCorrected)):
            segPool.segments[i].tReadNum = np.exp( yCorrected[i] +
                np.log(segPool.segments[i].nReadNum + 1)) - 1
            segPool.segments[i].log_ratio = np.log(
                (yCorrected[i] + 1.0) /
                (segPool.segments[i].nReadNum + 1.0)
            )

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(m, c)

    def _get_counts(self, tBamName, segPool):
        """
        get allele counts of target bam file
        save the counts into segPool
        """

        segNum = len(segPool.segments)
        processNum = self.__processNum
        print "processNum = {}".format(processNum)

        if processNum > segNum:
            processNum = segNum

        pool = Pool(processes=processNum)

        argsL = []

        for j in range(0, segNum):
            segName = segPool.segments[j].name
            chromName = segPool.segments[j].chromName
            chromIdx = segPool.segments[j].chromIdx
            start = segPool.segments[j].start
            end = segPool.segments[j].end

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
            pairedCountsJ, BAFCountsJ = countsTL[j]
            segPool.segments[j].pairedCounts = pairedCountsJ
            segPool.segments[j].BAFCounts = BAFCountsJ

# ===============================================================================
#  Function
# ===============================================================================


def process_by_segment(argsT):
    segName, chromName, chromIdx, start, end, nBamName,\
        tBamName, refFaName, minDepth, minBqual,\
        minMqual = argsT

    print 'Preprocessing segment {0}...'.format(segName)
    sys.stdout.flush()

    nBam = pysam.Samfile(nBamName, 'rb')
    tBam = pysam.Samfile(tBamName, 'rb')
    refFasta = pysam.Fastafile(refFaName)

    normalPileupIter = nBam.pileup(chromName, start, end)
    tumorPileupIter = tBam.pileup(chromName, start, end)

    pairedPileupIter = PairedPileupIterator(
        normalPileupIter, tumorPileupIter, start, end)
    pairedCountsIter = PairedCountsIterator(
        pairedPileupIter,
        refFasta,
        chromName,
        chromIdx,
        minDepth,
        minBqual,
        minMqual)

    pairedCountsJ, BAFCountsJ = iterator_to_counts(pairedCountsIter)
    countsTuple_j = (pairedCountsJ, BAFCountsJ)

    nBam.close()
    tBam.close()
    refFasta.close()

    return countsTuple_j


def iterator_to_counts(pairedCountsIter):
    buff = 100000

    pairedCountsJ = np.array([[], [], [], [], [], []], dtype=int).transpose()
    BAFCountsJ = np.zeros((100, 100))
    buffCounts = []
    i = 0

    for counts in pairedCountsIter:
        buffCounts.append(counts)
        i = i + 1

        if i < buff:
            continue

        buffCounts = np.array(buffCounts)

        if buffCounts.shape[0] != 0:
            BAFCountsBuff = get_BAF_counts(buffCounts)
            BAFCountsJ += BAFCountsBuff

        buffCountsFiltered = normal_heterozygous_filter(buffCounts)

        if buffCountsFiltered.shape[0] != 0:
            pairedCountsJ = np.vstack((pairedCountsJ, buffCountsFiltered))

        buffCounts = []
        i = 0

    buffCounts = np.array(buffCounts)

    if buffCounts.shape[0] != 0:
        BAFCountsBuff = get_BAF_counts(buffCounts)
        BAFCountsJ += BAFCountsBuff

    buffCountsFiltered = normal_heterozygous_filter(buffCounts)

    if buffCountsFiltered.shape[0] != 0:
        pairedCountsJ = np.vstack((pairedCountsJ, buffCountsFiltered))

    return (pairedCountsJ, BAFCountsJ)

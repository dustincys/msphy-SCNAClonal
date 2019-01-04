#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: data.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-11-24 13:53:26
#       History:
# =============================================================================
'''
import sys
from collections import Counter

import numpy as np
import scipy.cluster.hierarchy as hcluster

import phySCNAClonal.constants as constants
from phySCNAClonal.preprocess.data.elements.segment import Segment
from phySCNAClonal.preprocess.data.elements.stripe import Stripe
from phySCNAClonal.preprocess.data.pools.stripePool import StripePool
from phySCNAClonal.preprocess.utils import (BEDnParser, BEDParser,
                                            chrom_idx_to_name,
                                            chrom_name_to_idx,
                                            get_APM_frac_MAXMIN,
                                            get_APM_status, get_chrom_format,
                                            get_chrom_lens_idxs,
                                            get_cn_allele_config, get_LOH_frac,
                                            get_LOH_status, get_segment_name)


class SegmentPool:

    def __init__(self, idx = 0, maxCopyNumber = 6, coverage = 30):
        self.idx = idx

        self.maxCopyNumber = maxCopyNumber
        self.coverage = coverage

        self._alleleConfig = get_cn_allele_config(maxCopyNumber)

        self.baseline = -1  #此处baseline是对数值
        self.segments = []

    def load_seg_bed(self, bedName, containsReadNum=True):
        """
        load seg from bed, here, bed could be contains the gc corrected read
        number
        """

        if containsReadNum:
            bedChromL, bedStartL, bedEndL, tReadNumL, nReadNumL, gcL =\
                BEDnParser(bedName)
        else:
            bedChromL, bedStartL, bedEndL, gcL = BEDnParser(bedName)

        get_chrom_format(bedChromL)
        bedNum = len(bedChromL)

        for i in range(0, bedNum):
            chromIdx = chrom_name_to_idx(bedChromL[i])
            segName = get_segment_name(
                bedChromL[i], bedStartL[i], bedEndL[i])

            tempSeg = Segment()
            tempSeg.name = segName
            tempSeg.chromIdx = chromIdx
            tempSeg.chromName = bedChromL[i]
            tempSeg.start = bedStartL[i]
            tempSeg.end = bedEndL[i]

            if containsReadNum:
                tempSeg.nReadNum = nReadNumL[i]
                tempSeg.tReadNum = tReadNumL[i]

            tempSeg.gc = gcL[i]

            self.segments.append(tempSeg)

    def load_seg_bam(self, nBam, tBam, bedName):
        chromIdxL = constants.CHROM_IDX_LIST
        chromStart = constants.CHROM_START

        samSQ = nBam.header['SQ']
        samChromFormat = get_chrom_format(map(lambda x: x['SN'], samSQ))
        chromLenL, chromIdxL = get_chrom_lens_idxs(chromIdxL, samSQ)

        bedChromL, bedStartL, bedEndL, gcL = BEDParser(bedName)
        get_chrom_format(bedChromL)
        bedNum = len(bedChromL)

        for i in range(0, bedNum):
            chromIdx = chrom_name_to_idx(bedChromL[i])
            chromName = chrom_idx_to_name(chromIdx, samChromFormat)
            segName = get_segment_name(chromName, bedStartL[i], bedEndL[i])

            if chromIdx not in chromIdxL:
                print "Chromsome {0} not found, segment {1} excluded...".format(bedChromL[i], segName)
                sys.stdout.flush()
                continue

            chromListIdx = chromIdxL.index(chromIdx)

            if bedStartL[i] < chromStart or bedEndL[
                    i] > chromLenL[chromListIdx]:
                print 'Out of range chromsome {0}, segment {1} excluded...'.\
                    format(bedChromL[i], segName)
                sys.stdout.flush()
                continue

            nReadNum = nBam.count(
                chromName, bedStartL[i], bedEndL[i])
            tReadNum = tBam.count(
                chromName, bedStartL[i], bedEndL[i])

            tempSeg = Segment()
            tempSeg.name = segName
            tempSeg.chromIdx = chromIdx
            tempSeg.chromName = chromName
            tempSeg.start = bedStartL[i]
            tempSeg.end = bedEndL[i]
            tempSeg.nReadNum = nReadNum
            tempSeg.tReadNum = tReadNum
            tempSeg.gc = gcL[i]

            self.segments.append(tempSeg)

    def get_baseline(self, maxCopyNumber, subcloneNum, baselineThredLOH,
                     baselineThredAPM, mergeSeg=False, isPreprocess=False, index=-1):
        """
        return the baseline segment list
        """
        self._get_LOH_frac()
        self._get_LOH_status(baselineThredLOH, isPreprocess)
        self._get_APM_frac()
        self._get_APM_status(baselineThredAPM)

        self.baseline = self._calc_baseline_clusterd(maxCopyNumber,
                                                     subcloneNum,
                                                     isPreprocess)

        #debug
        debug = True
        if debug:
            if index == 0:
                self.baseline = -0.32
            else:
                self.baseline = -0.23

            deltaUp = 0.08
            deltaDown = 0.1
            blSegL = []

            for seg in self.segments:
                loga = np.log(seg.tReadNum + 1) - np.log(seg.nReadNum + 1)
                if loga > self.baseline - deltaDown and\
                        loga < self.baseline + deltaUp:
                    seg.tag = "BASELINE"
                    blSegL.append(seg)
                else:
                    seg.tag = "0"

            return blSegL

        else:
            if mergeSeg:
                return self._get_baseline_segs_from_stripe()
            else:
                return self._get_baseline_segs_by_label()


    def get_seg_by_tag(self, tag="BASELINE"):
        """
        return seg sub set by tag
        """
        return filter(lambda item: item.tag==tag, self.segments)

    def _get_baseline_segs_by_label(self):
        """
        return baseline segments by baseline label
        """
        return filter(lambda item: item.baselineLabel == 'TRUE', self.segments)

    def _get_baseline_segs_from_stripe(self):
        # 此处应该获取基线的条带，使用StripePool中的功能获取条带
        # def __init__(self, segPool, baseline=0.0, yDown, yUp, stripeNum, noiseStripeNum=2):
        # 获取yDown 和yUp 等相应的参数值，此处应该使用手动输入
        # print self.idx

        yDown = constants.YDOWNL[self.idx]
        yUp = constants.YUPL[self.idx]
        # 此处应该近似为最大拷贝数与亚克隆数的乘积，作为手工输入也可以
        stripeNum = constants.STRIPENUML[self.idx]
        noiseStripeNum = constants.NOISESTRIPENUML[self.idx]

        tempSP = StripePool(self, self.baseline, yDown, yUp, stripeNum,
                            noiseStripeNum)
        return tempSP.get_baseline_segs()

    def _calc_baseline(self, maxCopyNumber, subcloneNum, isPreprocess=False):
        """
        compute the Lambda S, no clustering
        """
        if not isPreprocess:
            print >> sys.stdout, "compute_Lambda_S function called from model"
            return

        rdRatioLog = []
        for j in range(0, len(self.segments)):
            if self.segments[j].APMStatus == 'TRUE' and\
                    self.segments[j].LOHStatus == 'FALSE':
                ratio = self.segments[j].tReadNum*1.0/\
                    self.segments[j].nReadNum
                rdRatioLog.append(np.log(ratio))

        rdRatioLog = np.array(rdRatioLog)
        if rdRatioLog.shape[0] == 0:
            print >> sys.stderr, 'Error: no APM-LOH position found, existing...'
            print >> sys.stderr, 'Either the baselineThredAPM is too large, or\
            the constants APM_N_MIN is too large; Or, the baseline_thred_LOH is\
            too small'
            sys.exit(1)

        rdrMinLog = np.min(rdRatioLog)

        for j in range(0, len(self.segments)):
            if self.segments[j].APMStatus == 'TRUE' and\
                    self.segments[j].LOHStatus == 'FALSE':
                ratio = self.segments[j].tReadNum*1.0/\
                    self.segments[j].nReadNum
                if rdrMinLog == np.log(ratio):
                    self.segments[j].baselineLabel = 'TRUE'
                else:
                    self.segments[j].baselineLabel = 'FALSE'
            else:
                self.segments[j].baselineLabel = 'FALSE'

        return rdrMinLog


    def _calc_baseline_clusterd(self, maxCopyNumber, subcloneNum, isPreprocess=False):
        """
        compute the Lambda S, through hierarchy clustering
        """
        if not isPreprocess:
            print >> sys.stdout, "compute_Lambda_S function called from model"
            return

        thresh = constants.HC_THRESH

        rdRatioLog = []
        for j in range(0, len(self.segments)):
            if self.segments[j].APMStatus == 'TRUE' and\
                    self.segments[j].LOHStatus != 'TRUE':
                ratio = self.segments[j].tReadNum*1.0/\
                    self.segments[j].nReadNum
                rdRatioLog.append(np.log(ratio))

        rdRatioLog = np.array(rdRatioLog)
        if rdRatioLog.shape[0] == 0:
            print >> sys.stderr, 'Error: no APM-LOH position found, existing...'
            print >> sys.stderr, 'Either the baselineThredAPM is too large, or\
            the constants APM_N_MIN is too large; Or, the baseline_thred_LOH is\
            too small'
            sys.exit(1)
        if rdRatioLog.shape[0] < 5:
            for j in range(0, len(self.segments)):
                if self.segments[j].APMStatus == 'TRUE' and\
                        self.segments[j].LOHStatus != 'TRUE':
                    self.segments[j].baselineLabel = 'TRUE'
                else:
                    self.segments[j].baselineLabel = 'FALSE'

            return np.mean(rdRatioLog)

        rdRatioLog = rdRatioLog.reshape(rdRatioLog.shape[0], 1)
        y = np.ones(rdRatioLog.shape)
        rdRatioLog = np.hstack((rdRatioLog, y))
        clusters = hcluster.fclusterdata(rdRatioLog,
                                         thresh,
                                         criterion="distance")
        mccs = Counter(clusters).most_common(maxCopyNumber * subcloneNum)

        rdrMinLog = float('Inf')
        clusterMin = -1
        for i in range(0, len(mccs)):
            clusterTemp = mccs[i][0]
            print >> sys.stdout, "cluster temp : {}".format(clusterTemp)
            tempRdrLog = np.mean(rdRatioLog[clusters == clusterTemp], axis=0)[0]
            print >> sys.stdout, "tempRdrLog"
            print >> sys.stdout, "log: {}".format(tempRdrLog)
            if rdrMinLog > tempRdrLog:
                rdrMinLog = tempRdrLog
                clusterMin = clusterTemp

        print >> sys.stdout, mccs
        print >> sys.stdout, "log baseline: {}".format(rdrMinLog)
        sys.stdout.flush()

        clusterFlag = (clusters == clusterMin)
        baselineNum = 0
        rdrIdx = 0
        for j in range(0, len(self.segments)):
            if self.segments[j].APMStatus == 'TRUE' and\
                    self.segments[j].LOHStatus != 'TRUE':
                if clusterFlag[rdrIdx]:
                    self.segments[j].baselineLabel = 'TRUE'
                    baselineNum = baselineNum + 1
                else:
                    self.segments[j].baselineLabel = 'FALSE'
                rdrIdx = rdrIdx + 1
            else:
                self.segments[j].baselineLabel = 'FALSE'

        print >> sys.stdout, "baselineNum: {}".format(baselineNum)

        if baselineNum == 0:
            print >> sys.stderr, 'Error: No diploid segments found, existing...'
            sys.exit(1)

        return rdrMinLog

    def _get_LOH_frac(self):
        for j in range(0, len(self.segments)):
            self.segments[j].LOHFrac = get_LOH_frac(
                self.segments[j].pairedCounts)

    def _get_APM_frac(self):
        cpbaf = constants.COVERAGE_PROPERATION_BAF
        for j in range(0, len(self.segments)):
            self.segments[j].APMFrac = get_APM_frac_MAXMIN(
                self.segments[j].pairedCounts, self.coverage,  cpbaf)

    def _get_LOH_status(self, baseThred, isPreprocess=False):
        if isPreprocess:
            LOHNum = 0
            FLOHNum = 0
            for j in range(0, len(self.segments)):
                self.segments[j].LOHStatus = get_LOH_status(
                    self.segments[j].LOHFrac, baseThred)
                if self.segments[j].LOHStatus == "TRUE":
                    LOHNum = LOHNum + 1
                elif self.segments[j].LOHStatus != "TRUE":
                    FLOHNum = FLOHNum + 1

            print >> sys.stdout, "LOHNum/segNum = {0}/{1}".format(LOHNum, len(self.segments))
            print >> sys.stdout, "FLOHNum/segNum = {0}/{1}".format(FLOHNum, len(self.segments))
        else:
            print >> sys.stdout, "get_LOH_status function called from model."

    def _get_APM_status(self, baselineThredAPM):
        APMNum = 0
        for j in range(0, len(self.segments)):
            self.segments[j].APMStatus = get_APM_status(
                self.segments[j].APMFrac, baselineThredAPM)
            if self.segments[j].APMStatus == "TRUE":
                APMNum = APMNum + 1

        print "APMNum/segNum = {0}/{1}".format(APMNum, len(self.segments))


    def output_txt(self, outFileName):
        """
        output txt file as input file of mh model
        """
        with open(outFileName, 'w') as outFile:
            outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                "name", "sid", "segsIdxL", "pairedCounts", "tReadNum",
                "nReadNum", "tag", "fixedC"))

            for s in self.segments:
                aT = s.pairedCounts[:,2]
                bT = s.pairedCounts[:,3]
                aTstrl = np.array_str(aT).strip("[]").split()
                bTstrl = np.array_str(bT).strip("[]").split()

                print s.name
                outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
                    "{0}_{1}_{2}".format(s.chromName, str(s.start), str(s.end)),
                    "{0}_{1}_{2}".format(s.chromName, str(s.start), str(s.end)),
                    ",",
                    "{0}|{1}".format(",".join(aTstrl), ",".join(bTstrl)),
                    s.tReadNum,
                    s.nReadNum,
                    s.tag,
                    s.fixedC))


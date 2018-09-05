#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: stripePool.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-03-02 18:20:21
#       History:
# =============================================================================
'''
from collections import Counter

import numpy as np
import scipy.cluster.hierarchy as hierarchy
from sklearn.cluster import MeanShift, estimate_bandwidth

from phySCNAClonal.preprocess.data.elements.stripe import Stripe
from phySCNAClonal.preprocess.plotGC import (BAFPlot, GCStripePlot,
                                             GCStripePoolPlot, plotMeanShift)
from phySCNAClonal.preprocess.utils import getBAFofSeg


import phySCNAClonal.constants as constants


class StripePool(object):
    """The stripe objects, including load, property operations"""

    def __init__(self, segPool, baseline, yDown, yUp, stripeNum, noiseStripeNum=2):
        """import segPool object

        :segPool: TODO

        """
        self.segPool = segPool
        self.baseline = baseline
        self._yDown = yDown
        self._yUp = yUp
        self.stripeNum = stripeNum
        self.noiseStripeNum = noiseStripeNum

        self.stripes = []  # stripes

    def get(self, byTag=False, plot=False):
        """
        if byTag, the output of Stripe should contains tag too.
        """
        self._aggregate(self._yDown, self._yUp, self.stripeNum, self.noiseStripeNum, byTag)
        if plot:
            gspp = GCStripePoolPlot(self)
            gspp.plot()

    def output_txt(self, outFileName):
        with open(outFileName, 'w') as outFile:
            outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(
                "name", "sid", "segsIdxL", "pairedCounts", "tReadNum",
                "nReadNum", "tag"))

            for s in self.stripes:
                aT = s.pairedCounts[:,2]
                bT = s.pairedCounts[:,3]
                aTstrl = np.array_str(aT).strip("[]").split()
                bTstrl = np.array_str(bT).strip("[]").split()

                print s.segsIdxL
                outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(
                    s.name,
                    s.sid,
                    ",".join(map(str, s.segsIdxL)),
                    "{0}|{1}".format(",".join(aTstrl), ",".join(bTstrl)),
                    s.tReadNum,
                    s.nReadNum,
                    s.tag))

    def _aggregate(self, yDown, yUp, stripeNum, noiseStripeNum=2, byTag=False,
                   plot=False):
        """The aggregation operations for segments in data
        :returns: stripes data structure
        """
        # function for out data and visualization
        # def writeToFile(self,cluster):
            # segments = self.segPool.segments
            # x = np.array(map(lambda seg:seg.gc,segments))
            # y = np.array(map(lambda seg: np.log(seg.tReadNum + 1) -
                                          # np.log(seg.nReadNum + 1), segments))
            # if not os.path.exists("plot_data/aggregate_data"):
                # os.makedirs("plot_data/aggregate_data")
            # outFile = open("plot_data/aggregate_data/aggregate_data.txt", "wr")
            # for i in range(0,len(cluster)):
                # outFile.write("{0}\t{1}\t{2}\n".format(x[i],y[i],cluster[i]))
            # outFile.close()

        assert stripeNum > 0

        rdRaioLog = []

        # here should keep idx
        ycV = np.array([
            np.log(seg.tReadNum + 1) - np.log(seg.nReadNum + 1)
            for seg in self.segPool.segments
        ])

        # 记录是否是outlier
        statusYcV = np.logical_and(ycV > yDown, ycV < yUp)

        yFcd = ycV.reshape(ycV.shape[0], 1)
        clusters = hierarchy.fclusterdata(
            yFcd, stripeNum + noiseStripeNum, criterion="maxclust",
            method="complete")

        mostClosedCluster, _ = self.__get_baseline_from_stripe(clusters, ycV)

        if plot:
            writeToFile(self, clusters)

        # 此处应该只获取最大和最小值之间的条带，且要保留原始位置，以方便索引
        # 此处获取最小和最大值之间的条带的方法是：直接去除这些位置不列入计算范围
        # 此处应该是去除了outlier之后的Counter

        mccs = Counter(
            clusters[statusYcV]).most_common(stripeNum + noiseStripeNum)

        for cId, _ in mccs:
            # 对每一个条带进行裂解操作，生成子条带, return
            if cId == mostClosedCluster:
                # 此处去除baseline
                continue
            self._decompose(cId, clusters, statusYcV, byTag, plot)

    def get_baseline_segs(self):
        """
        Get baseline segments, given the baseline value
        """

        yDown, yUp, stripeNum, noiseStripeNum =self._yDown, self._yUp,\
            self.stripeNum, self.noiseStripeNum

        assert stripeNum > 0

        rdRaioLog = []

        # here should keep idx
        ycV = np.array([
            np.log(seg.tReadNum + 1) - np.log(seg.nReadNum + 1)
            for seg in self.segPool.segments
        ])

        # 记录是否是outlier
        statusYcV = np.logical_and(ycV > yDown, ycV < yUp)

        yFcd = ycV.reshape(ycV.shape[0], 1)
        clusters = hierarchy.fclusterdata(
            yFcd, stripeNum + noiseStripeNum, criterion="maxclust",
            method="complete")

        _, blSegL = self.__get_baseline_from_stripe(clusters, ycV)
        # writeToFile(self,clusters)

        return blSegL

    def __get_baseline_from_stripe(self, clusters, ycV):
        """
        Get baseline segments in the aggregation step
        for the sake of
        the more the stripe approaches to baseline stripe, the smaller the gap
        it is
        """
        def get_cluster_centers(clusters, ycv):
            clustersUnique = np.unique(clusters)
            clusterCenters = {}
            for i in clustersUnique:
                index = np.where(clusters == i)
                clusterCenters[i] = np.mean(ycv[index[0]])
            return clusterCenters

        clusterCenters = get_cluster_centers(clusters, ycV)

        #寻找和baseline最接近的那一个center，标记其中的segment为baseline
        mostClosedCluster = 0
        minDis = float("Inf")

        for key in clusterCenters.keys():
            tempDis = np.abs(self.baseline - clusterCenters[key])
            if minDis > tempDis:
                minDis = tempDis
                mostClosedCluster = key

        index = np.where(clusters == mostClosedCluster)
        for i in index[0]:
            self.segPool.segments[i].tag = "BASELINE"

        blSegL = filter(lambda item:item.tag == "BASELINE",
                        self.segPool.segments)

        return mostClosedCluster, blSegL


    def _decompose(self, cId, clusters, statusYcV, byTag=False, plot=False):
        """The decomposition operations for segments in data

        :parameters: TODO
        :returns: TODO

        """
        # 获得该类别的所有结点idx：
        # 即，clusters 中与cId相等且，在statusYcV中的位置
        ca = np.argwhere(clusters == cId).flatten()
        sa = np.argwhere(statusYcV).flatten()
        mSIdx = np.intersect1d(ca, sa)

        # 这里的基于BAF的归类处理分为3个步骤

        # 首先进行所有seg的BAF的密度估计，然后获得峰值    类别定位
        # 然后对每一个seg进行归类，按照内部投票的方式     Seg归类
        # 然后返回

        # 这里需要有一个记录原始向量中位置的向量
        # 此处需要对空seg 去除空pairedCounts

        segsIDL = [(self.segPool.segments[idx], idx) for idx in mSIdx if
                   len(self.segPool.segments[idx].pairedCounts) > 1]

        segsLnoPCIDL = [(self.segPool.segments[idx], idx) for idx in mSIdx if
                        len(self.segPool.segments[idx].pairedCounts) == 0]
        # 此处需要进行每一个seg中的BAF统一投票之后在进行聚类 #
        # 此处决定是否是用最大最小限制
        # status_p_T_v = np.logical_and(pT > p_T_min, pT < p_T_max).flatten()
        if len(segsLnoPCIDL) > 0:
            if not byTag:
                tempStripe = Stripe()
                tempStripe.sid = "{0}".format(str(cId))
                subSegL, subSegIdxL = map(list, zip(*segsLnoPCIDL))
                tempStripe.init_segs(subSegL, subSegIdxL)
                self.stripes.append(tempStripe)
                for segIdx in subSegIdxL:
                    self.segPool.segments[segIdx].stripeID = tempStripe.sid
            else:
                tempTags = set([seg.tag for seg, idx in segsLnoPCIDL])
                for tempTag in tempTags:
                    tagIdx = 0
                    if "BASELINE" == tempTag:
                        continue

                    tempL = [(seg, idx) for seg, idx in segsLnoPCIDL if seg.tag == tempTag]
                    subSegL, subSegIdxL = map(list, zip(*tempL))
                    tempStripe = Stripe()
                    tempStripe.name = "{0}_{1}_{2}".format(str(cId),
                                                            tempTag,
                                                            str(tagIdx))
                    tempStripe.sid = "{0}_{1}_{2}".format(str(cId),
                                                            tempTag,
                                                            str(tagIdx))

                    tempStripe.init_segs(subSegL, subSegIdxL)
                    tempStripe.tag = tempTag
                    self.stripes.append(tempStripe)
                    for segIdx in subSegIdxL:
                        self.segPool.segments[segIdx].stripeID = tempStripe.sid
                    tagIdx = tagIdx + 1

        if len(segsIDL) > 0:

            #############################################
            #  We only decompose large number segments  #
            #############################################
            decomposeNumThresh = constants.DECOMPOSE_NUMBER_THRESHOLD

            segL, segIdxL = map(list, zip(*segsIDL))
            segsBAFL = map(getBAFofSeg, segL)

            if len(segL) < decomposeNumThresh:
                if not byTag:
                    tempStripe = Stripe()
                    tempStripe.sid = "{0}".format(str(cId))
                    tempStripe.init_segs(segL, segIdxL)
                    self.stripes.append(tempStripe)
                    for segIdx in segIdxL:
                        self.segPool.segments[segIdx].stripeID = tempStripe.sid
                else:
                    tempTags = set([seg.tag for seg in segL])
                    for tempTag in tempTags:
                        tagIdx = 0
                        if "BASELINE" == tempTag:
                            continue

                        tempL = [(seg, idx) for seg, idx in zip(segL, segIdxL)
                                 if seg.tag == tempTag]
                        subSegL, subSegIdxL = map(list, zip(*tempL))
                        tempStripe = Stripe()
                        tempStripe.name = "{0}_{1}_{2}".format(str(cId),
                                                                tempTag,
                                                                str(tagIdx))
                        tempStripe.sid = "{0}_{1}_{2}".format(str(cId),
                                                                tempTag,
                                                                str(tagIdx))

                        tempStripe.init_segs(subSegL, subSegIdxL)
                        tempStripe.tag = tempTag
                        self.stripes.append(tempStripe)
                        for segIdx in subSegIdxL:
                            self.segPool.segments[segIdx].stripeID = tempStripe.sid
                        tagIdx = tagIdx + 1
            else:
                pT = np.array(segsBAFL)
                pT = pT.reshape(pT.shape[0], 1)
                y = np.ones(pT.shape)
                pTy = np.hstack((pT, y))
                pTy = pTy[~np.isnan(pTy).any(axis=1)]

                bandwidth = estimate_bandwidth(pTy, quantile=0.2, n_samples=500)
                if bandwidth == 0:
                    bandwidth = 1
                ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
                ms.fit(pTy)
                labels = ms.labels_
                clusterCenters = ms.cluster_centers_

                if plot:
                    plotMeanShift(cId, pTy, labels, clusterCenters)

                segLabelL = [
                    self._getSegLabl(seg, clusterCenters) for seg in segL
                ]

                tempL = zip(segL, segIdxL, segLabelL)

                # 注意此处要对目标stripe的seg进行tag分类,分类之后在生成条带
                for label in set(segLabelL):
                    if label == -1:
                        continue

                    subTempL = [(seg, idx) for seg, idx, l in tempL if l == label]

                    subSegL, subSegIdxL = map(list, zip(*subTempL))

                    if not byTag:
                        tempStripe = Stripe()
                        tempStripe.sid = "{0}_{1}".format(str(cId), str(label))
                        tempStripe.init_segs(subSegL, subSegIdxL)
                        self.stripes.append(tempStripe)
                        for segIdx in subSegIdxL:
                            self.segPool.segments[segIdx].stripeID = tempStripe.sid
                    else:
                        tempTags = set([seg.tag for seg in subSegL])
                        for tempTag in tempTags:
                            tagIdx = 0
                            if "BASELINE" == tempTag:
                                continue
                            subSubTempL = [(seg, idx) for seg, idx in subTempL if
                                        seg.tag == tempTag ]
                            subSubSegL, subSubSegIdxL = map(list, zip(*subSubTempL))
                            tempStripe = Stripe()
                            tempStripe.name = "{0}_{1}_{2}_{3}".format(str(cId),
                                                                    label,
                                                                    tempTag,
                                                                    str(tagIdx))
                            tempStripe.sid = "{0}_{1}_{2}_{3}".format(str(cId),
                                                                    label,
                                                                    tempTag,
                                                                    str(tagIdx))

                            tempStripe.init_segs(subSubSegL, subSubSegIdxL)
                            tempStripe.tag = tempTag
                            self.stripes.append(tempStripe)
                            for segIdx in subSubSegIdxL:
                                self.segPool.segments[segIdx].stripeID = tempStripe.sid
                            tagIdx = tagIdx + 1

        # merge baseline, or not baseline in the stripe? toggle
        #  TODO: check out
        if False and byTag:
            blSegL = [seg for seg in self.segPool.segments if "BASELINE" ==
                    seg.tag]
            blSegIdxL = [idx for idx, seg in enumerate(self.segPool.segments) if "BASELINE" ==
                    seg.tag]
            tempStripe = Stripe()
            tempStripe.sid = "{0}_{1}_{2}".format(str(-1), str(idx), "BASELINE")
            tempStripe.init_segs(blSegL, blSegIdxL)
            tempStripe.tag = "-999999"
            self.stripes.append(tempStripe)

    def _getSegLabl(self, seg, clusterCenters):
        if seg.pairedCounts is None:
            return -1

        aTseg = seg.pairedCounts[:, 2]
        bTseg = seg.pairedCounts[:, 3]
        dTseg = aTseg + bTseg
        lTseg = np.min(seg.pairedCounts[:, 2:4], axis=1)
        pTseg = lTseg * 1.0 / dTseg

        disSeg = np.abs(pTseg[:, None] - clusterCenters[:, 0])
        labelsSeg = np.argmin(disSeg, axis=1)

        print Counter(labelsSeg)
        print Counter(labelsSeg).most_common(1)
        if len(Counter(labelsSeg).most_common(1)) == 0:
            print "counter equals 0!"

        return Counter(labelsSeg).most_common(1)[0][0]

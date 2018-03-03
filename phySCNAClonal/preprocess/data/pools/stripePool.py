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

from phySCNAClonal.preprocess.elements.stripe import Stripe


class StripePool(object):
    """The stripe objects, including load, property operations"""

    def __init__(self, segPool, baseline, yDown, yUp, stripeNum, noiseStripeNum=2):
        """import segPool object

        :segPool: TODO

        """
        self._segPool = segPool
        self._baseline = baseline
        self._yDown = yDown
        self._yUp = yUp
        self._stripeNum = stripeNum
        self._noiseStripeNum = noiseStripeNum

        self.stripes = []  # stripes

        self.baseline = baseline

    def get(self, byTag=False):
        """
        if byTag, the output of Stripe should contains tag too.
        """
        self._aggregate(self._yDown, self._yUp, self._stripeNum, self._noiseStripeNum, byTag)

    def output_txt(self, outFileName):
        with open(outFileName, 'w') as outFile:
            outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                "id", "segsIdxL", "pairedCounts", "tReadNum",
                "nReadName", "tag"))

            for s in self.stripes:
                aT = s.pairedCounts[:,2]
                bT = s.pairedCounts[:,3]
                aTstrl = np.array_str(aT).strip("[]").split()
                bTstrl = np.array_str(bT).strip("[]").split()

                outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    s.id,
                    ",".join(s.segsIdxL),
                    "{0}|{1}".format(",".join(aTstrl), ",".join(bTstrl)),
                    s.tReadNum,
                    s.nReadName,
                    s.tag)

    def _aggregate(self, yDown, yUp, stripeNum, noiseStripeNum=2, byTag=False):
        """The aggregation operations for segments in data
        :returns: stripes data structure
        """
        assert stripeNum > 0

        rdRaioLog = []

        # here should keep idx
        ycV = np.array([
            np.log(seg.tReadNum + 1) - np.log(seg.nReadName + 1)
            for seg in self._segPool.segments
        ])

        # 记录是否是outlier
        statusYcV = np.logical_and(ycV > yDown, ycV < yUp)

        yFcd = ycV.reshape(ycV.shape[0], 1)
        clusters = hierarchy.fclusterdata(
            yFcd, stripeNum + noiseStripeNum, criterion="distance")

        # 此处应该只获取最大和最小值之间的条带，且要保留原始位置，以方便索引
        # 此处获取最小和最大值之间的条带的方法是：直接去除这些位置不列入计算范围
        # 此处应该是去除了outlier之后的Counter

        mccs = Counter(
            clusters[statusYcV]).most_common(stripeNum + noiseStripeNum)

        for cId, _ in mccs:
            # 对每一个条带进行裂解操作，生成子条带, return
            self._decompose(cId, clusters, statusYcV, byTag)

    def _decompose(self, cId, clusters, statusYcV, byTag=False):
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
        segsL = [self._segPool.segments[idx] for idx in mSIdx]

        pairedCountsAll = np.array(
            [[], [], [], [], [], []], dtype=int).transpose()
        for seg in segsL:
            pairedCountsAll = np.vstack((pairedCountsAll,
                                           seg.pairedCounts))

        aT = pairedCountsAll[:, 2]
        bT = pairedCountsAll[:, 3]
        dT = aT + bT
        lT = np.min(pairedCountsAll[:, 2:4], axis=1)
        pT = lT * 1.0 / dT

        # 此处决定是否是用最大最小限制
        # status_p_T_v = np.logical_and(pT > p_T_min, pT < p_T_max).flatten()

        y = np.ones(pT.shape)
        pTy = np.hstack((pT, y))
        bandwidth = estimate_bandwidth(pTy, quantile=0.2, n_samples=500)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(pTy)
        labels = ms.labels_
        clusterCenters = ms.clusterCenters
        # labelsUnique = np.unique(labels)
        # nClusters = len(labelsUnique)

        segLabelL = [
            self._getSegLabl(seg, clusterCenters) for seg in segsL
        ]

        # 注意此处要对目标stripe的seg进行tag分类,分类之后在生成条带
        for label in set(segLabelL):
            if label == -1:
                continue
            subSegL = [
                seg for idx, seg in enumerate(segsL)
                if segLabelL[idx] == label
            ]
            subSegIdxL = [
                mSIdx[idx] for idx, seg in enumerate(segsL)
                if segLabelL[idx] == label
            ]
            if not byTag:
                tempStripe = Stripe()
                tempStripe.id = "{0}".format(str(cId))
                tempStripe.init_segs(subSegL, subSegIdxL)
                self.stripes.append(tempStripe)
            else:
                tempTags = set([seg.tag for seg in subSegL])
                for tempTag in tempTags:
                    if "BASELINE" == tempTag:
                        continue

                    subSubSegL = [seg for seg in subSegL if seg.tag ==
                                     tempTag]
                    subSubSegIdxL = [ mSIdx[idx] for idx, seg in enumerate(segsL)
                                    if segLabelL[idx] == label and
                                    seg.tag == tempTag ]
                    tempStripe = Stripe()
                    tempStripe.id = "{0}_{1}".format(str(cId), tempTag)
                    tempStripe.init_segs(subSubSegL, subSubSegIdxL)
                    # if byTag, stripe contains tag too
                    tempStripe.tag = tempTag
                    self.stripes.append(tempStripe)

        # merge baseline
        if byTag:
            blSegL = [seg for seg in self._segPool.segments if "BASELINE" ==
                    seg.tag]
            blSegIdxL = [idx for idx, seg in enumerate(self._segPool.segments) if "BASELINE" ==
                    seg.tag]
            tempStripe = Stripe()
            tempStripe.id = "{0}_{1}".format(str(-1), str(idx), "BASELINE")
            tempStripe.init_segs(blSegL, blSegIdxL)
            tempStripe.tag = "BASELINE"
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

        return Counter(labelsSeg).most_common(1)[0][0]

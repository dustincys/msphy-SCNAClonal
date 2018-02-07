'''
# =============================================================================
#      FileName: stripe.py
#          Desc: Stripe aggregation and decomposition
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-11-24 13:53:39
#       History:
# =============================================================================
'''

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import heapq
import sys
from collections import Counter
from random import randint

import numpy as np
from scipy.cluster import hierarchy
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
from scipy.stats.mstats import gmean
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.datasets.samples_generator import make_blobs

import constants
from pydp.densities import Density, log_poisson_pdf
from utils import (get_cn_allele_config, get_loga, get_mu_E_joint,
                   log_binomial_likelihood, mad_based_outlier)

# import json

class Stripe:
    def __init__(self):
        # 记录对应原始seg的索引
        self.name = ""
        self.id = ""
        self.segIdL = None

        self.pairedCounts = None

        self.nReadName = -1.0
        self.tReadNum = -1.0

        # 似乎用不到，实际上就是:
        # self.tReadNum /self.nReadName
        # self.rdr = -1.0

        self.tag = "1"

        # 这两个应该放在这里
        self.copyNumber = -1
        self.genotype = ""

        # phi应该放在node结点中
        # self.phi = 0.0

        self.tssb = None
        self.node = None  # this is the node where the datum resides

    def init_segs(self, segList, segIdL):
        self.segIdL = segIdL

        self._init_RD(segList)
        self._init_BAF(segList)

    def _init_RD(self, segList):
        # 获取几何平均值
        tReadNum = [seg.tReadNum for seg in segList]
        nReadName = [seg.nReadName for seg in segList]

        self.tReadNum = gmean(tReadNum)
        self.nReadName = gmean(nReadName)

    def _init_BAF(self, segList):
        self.pairedCounts = np.array(
            [[], [], [], [], [], []], dtype=int).transpose()

        for seg in segList:
            self.pairedCounts = np.vstack((self.pairedCounts,
                                            seg.pairedCounts))

    def _log_likelihood(self, phi, update_tree=True):
        if update_tree:
            ##################################################
            # some useful info about the tree,
            # used by CNV related computations,
            u.set_node_height(self.tssb)
            u.set_path_from_root_to_node(self.tssb)
            u.map_datum_to_node(self.tssb)
            ##################################################

        # 注意：此处应该不受CN\genotype的限制，但不记录目标的CN和genotype
        # 因为此处的parameter不是准确的phi,所以无法达到最优，但为了能够抽样得到
        # 最佳结构。此处要使CN和genotype自由发挥

        # 在时间上的先后顺序能够明确影响测序数据的位置，只有重叠位置的时候才会
        # 发生

        return self.__log_likelihood_RD_BAF(phi)


    def __log_likelihood_RD_BAF(self, phi, baseline):
        # 此处是否添加记录
        copyNumbers = None
        if seg.tag == "True":
            copyNumbers = [2]
        elif get_loga(seg) > baseline:
            copyNumbers = range(2, self._maxCopyNumber + 1)
        else:
            copyNumbers = range(0, 2 + 1)

        llPiS = [self._getLLStripe(copyNumber, phi) for copyNumber in
                   copyNumbers]
        ll, pi = max(llPiS, key=lambda x: x[0])
        cn = llPiS.index((ll, pi))

        self.copyNumber = cn
        self.genotype = pi

        return ll


    def _getLLStripe(self, copyNumber, phi):
        rdWeight = constants.RD_WEIGHT_TSSB

        llStripe = 0
        llRd = self._getRD(copyNumber, phi)
        alleleTypes = self._allele_config[copyNumber]

        # remove the weak baf point
        self._augBAF(copyNumber)

        if 0 == self.pairedCounts.shape[0]:
            llBAF = 0
            pi = "*"
        else:
            llBAF, pi = self._getBAF(self, copyNumber, alleleTypes, phi)

        llStripe = llRd * rdWeight + llBAF * (1 - rdWeight)

        return llStripe, pi


    def _augBAF(self, copyNumber):
        # todo: remove the baf point is not a good idea
        threshold = constants.BAF_THRESHOLD * self._coverage

        if copyNumber > 2:
            dTj = np.sum(self.pairedCounts[:, 2:4], axis=1)
            idxRm = tuple(np.where(dTj < threshold)[0])
            self.pairedCounts = np.delete(self.pairedCounts, idxRm, axis=0)
        else:
            pass

    def _getRD(self, copyNumber, phi):
        cN = constants.COPY_NUMBER_NORMAL

        barC = phi * copyNumber + (1.0 - phi) * cN

        lambdaPossion = (
            barC / cN) * self._baseline * (seg.nReadName + 1)
        if lambdaPossion < 0:
            lambdaPossion = 0

        llRD = log_poisson_pdf(seg.tReadNum, lambdaPossion)
        return llRD

    def _getBAF(self, copyNumber, alleleTypes, phi):
        cN = constants.COPY_NUMBER_NORMAL
        muN = constants.MU_N

        muG = np.array(alleleTypes.values())
        muE = get_mu_E_joint(muN, muG, cN, copyNumber, phi)

        if seg.pairedCounts.shape[0] > 1:
            bTj = np.min(seg.pairedCounts[:, 2:4], axis=1)
            dTj = np.sum(seg.pairedCounts[:, 2:4], axis=1)
            baf = bTj * 1.0 / dTj
            outlier = mad_based_outlier(baf)
            BAF = np.delete(seg.pairedCounts, list(outlier.astype(int)), axis=0)
            bTj = np.min(BAF[:, 2:4], axis=1)
            dTj = np.sum(BAF[:, 2:4], axis=1)

        else:
            bTj = np.min(seg.pairedCounts[:, 2:4], axis=1)
            dTj = np.sum(seg.pairedCounts[:, 2:4], axis=1)
            pass

        ll = log_binomial_likelihood(bTj, dTj, muE)
        llBAFs = ll.sum(axis=0)
        idxMax = llBAFs.argmax(axis=0)
        llBAF = llBAFs[idxMax]
        pi = alleleTypes[alleleTypes.keys()[idxMax]]
        return llBAF, pi


class StripePool(object):
    """The stripe objects, including load, property operations"""

    def __init__(self, segmentPool):
        """import segmentPool object

        :segmentPool: TODO

        """
        self._segmentPool = segmentPool
        self.stripes = []  # stripes

        self.baseline = -1

    def get(self, yDown, yUp, stripeNum, noiseStripeNum=2):
        """TODO: Docstring for get.
        :returns: TODO

        """
        self._aggregate(yDown, yUp, stripeNum, noiseStripeNum=2)

    def output_txt(self, outFileName):
        with open(outFileName, 'w') as outFile:
            outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                "id", "segIdL", "pairedCounts", "tReadNum",
                "nReadName", "tag"))

            for s in self.stripes:
                aT = s.pairedCounts[:,2]
                bT = s.pairedCounts[:,3]
                aTstrl = np.array_str(aT).strip("[]").split()
                bTstrl = np.array_str(bT).strip("[]").split()

                outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    s.id,
                    ",".join(s.segIdL),
                    "{0}|{1}".format(",".join(aTstrl), ",".join(bTstrl)),
                    s.tReadNum,
                    s.nReadName,
                    s.tag)
            pass

    def _aggregate(self, yDown, yUp, stripeNum, noiseStripeNum=2):
        """The aggregation operations for segments in data

        :returns: stripes data structure

        """
        assert stripeNum > 0

        rdRaioLog = []

        # here should keep idx
        ycV = np.array([
            np.log(seg.tReadNum + 1) - np.log(seg.nReadName + 1)
            for seg in self._segmentPool.segments
        ])

        # 记录是否是outlier
        statusYcV = np.logical_and(ycV > y_min, ycV < y_max)

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
            self._decompose(cId, clusters, statusYcV)

    def _decompose(self, cId, clusters, statusYcV):
        """The decomposition operations for segments in data

        :parameters: TODO
        :returns: TODO

        """
        # 获得该类别的所有结点idx：
        # 即，clusters 中与cId相等且，在statusYcV中的位置
        ca = np.argwhere(clusters == cId).flatten()
        sa = np.argwhere(statusYcV).flatten()
        mSIdx = np.intersectid(ca, sa)

        # 这里的基于BAF的归类处理分为3个步骤

        # 首先进行所有seg的BAF的密度估计，然后获得峰值    类别定位
        # 然后对每一个seg进行归类，按照内部投票的方式     Seg归类
        # 然后返回

        # 这里需要有一个记录原始向量中位置的向量
        segList = [self._segmentPool.segments[idx] for idx in mSIdx]

        pairedCountsAll = np.array(
            [[], [], [], [], [], []], dtype=int).transpose()
        for seg in segList:
            pairedCountsAll = np.vstack((pairedCountsAll,
                                           seg.pairedCounts))

        aT = pairedCountsAll[:, 2]
        bT = pairedCountsAll[:, 3]
        dT = aT + bT
        lT = np.min(pairedCountsAll[:, 2:4], axis=1)
        pT = lT * 1.0 / dT

        # status_p_T_v = np.logical_and(pT > p_T_min, pT < p_T_max).flatten()

        y = np.ones(pT.shape)
        pTy = np.hstack((pT, y))
        bandwidth = estimate_bandwidth(pTy, quantile=0.2, n_samples=500)
        ms.fit(X)
        labels = ms.labels_
        clusterCenters = ms.clusterCenters
        labelsUnique = np.unique(labels)
        nClusters = len(labelsUnique)

        segLabel = [
            self._getSegLabl(seg, clusterCenters) for seg in segList
        ]

        for label in set(segLabel):
            if label == -1:
                continue
            subSegList = [
                seg for seg, idx in enumerate(segList)
                if segLabel[idx] == label
            ]
            subSegIdx = [
                mSIdx[idx] for seg, idx in enumerate(segList)
                if segLabel[idx] == label
            ]
            tempStripe = Stripe()
            tempStripe.id = "{0}_{1}".format(str(cId), str(idx))
            tempStripe.init_segs(subSegList, subSegIdx)
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

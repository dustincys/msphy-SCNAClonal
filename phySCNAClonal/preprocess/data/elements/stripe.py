#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

import numpy as np
from scipy.stats.mstats import gmean

import phySCNAClonal.constants as constants
from phySCNAClonal.model.util2 import (map_datum_to_node, set_node_height,
                                       set_path_from_root_to_node)
from phySCNAClonal.preprocess.utils import (get_loga, get_mu_E_joint,
                                            log_binomial_likelihood,
                                            mad_based_outlier)
from pydp.densities import log_poisson_pdf


# import json

class Stripe:
    def __init__(self):
        # 记录对应原始seg的索引
        self.id = -1
        self.name = ""
        self.sid = ""
        self.segsIdxL = None

        self.pairedCounts = None

        self.nReadNum = -1.0
        self.tReadNum = -1.0

        # 似乎用不到，实际上就是:
        # self.tReadNum /self.nReadNum
        # self.rdr = -1.0

        self.tag = "0"

        # 这两个应该放在这里
        self.copyNumber = -1
        self.genotype = ""

        # phi应该放在node结点中
        self.phi = 0.0
        self.fixedC = -1

        self.tssb = None
        self.node = None  # this is the node where the datum resides

    def init_segs(self, segsL, segsIdxL):
        self.segsIdxL = segsIdxL

        self._init_RD(segsL)
        self._init_BAF(segsL)

    def _init_RD(self, segsL):
        # 获取几何平均值
        tReadNum = [seg.tReadNum for seg in segsL]
        nReadNum = [seg.nReadNum for seg in segsL]

        self.tReadNum = gmean(tReadNum)
        self.nReadNum = gmean(nReadNum)

    def _init_BAF(self, segsL):
        self.pairedCounts = np.array(
            [[], [], [], [], [], []], dtype=int).transpose()

        for seg in segsL:
            self.pairedCounts = np.vstack((self.pairedCounts,
                                            seg.pairedCounts))

    def _log_likelihood(self,
                        phi,
                        alleleConfig,
                        baseline,
                        maxCopyNumber,
                        update_tree=True):
        if update_tree:
            ##################################################
            # some useful info about the tree,
            # used by CNV related computations,
            set_node_height(self.tssb)
            set_path_from_root_to_node(self.tssb)
            map_datum_to_node(self.tssb)
            ##################################################

        # 注意：此处应该不受CN\genotype的限制，但不记录目标的CN和genotype
        # 因为此处的parameter不是准确的phi,所以无法达到最优，但为了能够抽样得到
        # 最佳结构。此处要使CN和genotype自由发挥

        # 在时间上的先后顺序能够明确影响测序数据的位置，只有重叠位置的时候才会
        # 发生
        # 注意此处可以设置默认参数

        return self.__log_likelihood_RD_BAF(phi,
                                            alleleConfig,
                                            baseline,
                                            maxCopyNumber,fixedC)


    def __log_likelihood_RD_BAF(self, phi, alleleConfig, baseline,
                                maxCopyNumber, fixedC=None):
        # 此处是否添加记录
        copyNumbers = None
        # 此处需要确认是否是使用默认的baseline为tag
        if fixedC is not  None:
            copyNumbers = [fixedC]
        else:
            if self.tag == "BASELINE":
                copyNumbers = [2]
            elif get_loga(self) > baseline:
                copyNumbers = range(3, maxCopyNumber + 1)
            else:
                copyNumbers = range(0, 2)

        llPiS = [self._getLLStripe(copyNumber, phi, baseline, alleleConfig) for copyNumber in
                   copyNumbers]
        ll, pi = max(llPiS, key=lambda x: x[0])
        cn = llPiS.index((ll, pi))

        self.copyNumber = cn
        self.genotype = pi

        return ll

    def _getLLStripe(self, copyNumber, phi, baseline, alleleConfig):
        rdWeight = constants.RD_WEIGHT_TSSB

        llStripe = 0
        llRd = self._getRD(copyNumber, phi, baseline)
        alleleTypes = alleleConfig[copyNumber]

        # remove the weak baf point
        self._augBAF(copyNumber)

        if 0 == self.pairedCounts.shape[0]:
            llBAF = 0
            pi = "*"
        else:
            llBAF, pi = self._getBAF(copyNumber, alleleTypes, phi)

        llStripe = llRd * rdWeight + llBAF * (1 - rdWeight)

        return llStripe, pi

    def _augBAF(self, copyNumber):
        # todo: remove the baf point is not a good idea
        threshold = constants.BAF_THRESHOLD * constants.COVERAGE

        if copyNumber > 2:
            dTj = np.sum(self.pairedCounts[:, 2:4], axis=1)
            idxRm = tuple(np.where(dTj < threshold)[0])
            self.pairedCounts = np.delete(self.pairedCounts, idxRm, axis=0)
        else:
            pass

    def _getRD(self, copyNumber, phi, baseline):
        cN = constants.COPY_NUMBER_NORMAL
        cMIN = constants.MINIMUM_POSITIVE

        barC = phi * copyNumber + (1.0 - phi) * cN

        lambdaPossion = (barC / cN) * np.exp(baseline) * (self.nReadNum + 1)
        if lambdaPossion <= 0:
            lambdaPossion = cMIN

        # print (self.tReadNum, lambdaPossion)
        llRD = log_poisson_pdf(self.tReadNum, lambdaPossion)
        return llRD

    def _getBAF(self, copyNumber, alleleTypes, phi):
        cN = constants.COPY_NUMBER_NORMAL
        muN = constants.MU_N

        muG = np.array(alleleTypes.values())
        muE = get_mu_E_joint(muN, muG, cN, copyNumber, phi)

        if self.pairedCounts.shape[0] > 1:
            bTj = np.min(self.pairedCounts[:, 2:4], axis=1)
            dTj = np.sum(self.pairedCounts[:, 2:4], axis=1)
            baf = bTj * 1.0 / dTj
            outlier = mad_based_outlier(baf)
            BAF = np.delete(self.pairedCounts, np.where(outlier)[0], axis=0)
            bTj = np.min(BAF[:, 2:4], axis=1)
            dTj = np.sum(BAF[:, 2:4], axis=1)

        else:
            bTj = np.min(self.pairedCounts[:, 2:4], axis=1)
            dTj = np.sum(self.pairedCounts[:, 2:4], axis=1)
            pass

        ll = log_binomial_likelihood(bTj, dTj, muE)
        llBAFs = ll.sum(axis=0)
        idxMax = llBAFs.argmax(axis=0)
        llBAF = llBAFs[idxMax]
        pi = alleleTypes.keys()[idxMax]
        return llBAF, pi

#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: stripenode.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-03-08 15:21:59
#       History:
# =============================================================================
'''
import scipy.stats as stat
from numpy.random import rand
from scipy.misc import comb
from scipy.stats import beta, binom

import phySCNAClonal.constants as constants
from phySCNAClonal.model.node import Node
# 在节点中保存varphiR\piR\remainR
from phySCNAClonal.model.usegsampler.segsupportive import MultiRangeSampler


class DataNode(Node):

    initMean = 0.5
    minConc = 0.01
    maxConc = 0.1

    def __init__(self, parent=None, tssb=None, conc=0.1):
        super(DataNode, self).__init__(parent=parent, tssb=tssb)

        # pi is a first-class citizen
        self.pi = 0.0
        self.param = 0.0
        self.param1 = 0.0
        self.pi1 = 0.0  # used in MH	to store old state

        self.path = None  # set of nodes from root to this node
        self.ht = 0.0

        if parent is None:
            self._conc = conc
            self.pi = 1.0
            self.param = 1.0

        else:
            self.pi = rand(1)*parent.pi
            parent.pi = parent.pi - self.pi
            self.param = self.pi

        self.varphiR = MultiRangeSampler(0,1)
        self.piR = MultiRangeSampler(0,1)
        self.epsilon = ""

    def conc(self):
        if self.parent() is None:
            return self._conc
        else:
            return self.parent().conc()

    def kill(self):
        if self._parent is not None:
            self._parent._children.remove(self)
        self._parent.pi = self._parent.pi + self.pi
        self._parent = None
        self._children = None

    def logprob(self, x, alleleConfig, baseline, maxCopyNumber):
        # 此处需要进一步处理
        return x[0]._log_likelihood(self.param,
                                    alleleConfig,
                                    baseline,
                                    maxCopyNumber)

    def logprob_restricted(self, x, alleleConfig, baseline, maxCopyNumber):
        # 此处添加time stamp 的限制
        # 在投掷过程中，如果出现一下情况，需要重新投掷
        #
        # 如果与该结点内的tag不一致为log(0)
        #
        # 如果tag小于其所有父结点任意一个或大于子结点任意一个为log(0)
        #
        # 如果出现结点内部的stripe的gap过于靠近
        #
        # The data x is added into the node, by default
        #
        # 注意此处，不必须限制同一个节点之内的拷贝数必须不相同
        # 可能有相同拷贝数但不同的基因型的状况
        # if self.__is_good_gaps(x):
        #
        # cmpTag = self.__cmp_tags(x)
        # if cmpTag == -1:
            # return -float('Inf')
        # elif cmpTag == 1:
            # return -float('Inf')
        # else:

        cmpPedigree = self.__cmp_pedigree(x)
        if cmpPedigree == -1:
            return -float('Inf')
        elif cmpPedigree == 1:
            return -float('Inf')
        else:
            return self.logprob(x, alleleConfig, baseline, maxCopyNumber)

        #######################################################
        #  here -float("Inf") means the situation restricted  #
        #######################################################

    def __cmp_pedigree(self, x):
        """is the time tag descending along the pedigree?

        @return:  Flag
        @rtype :  int
        """
        #  TODO: check out the len(n.data) #
        #  here, should not return itself

        tagCurX = int(x[0].tag)

        children = self.children()
        if 0 < len(children):
            offsprings = []
            for child in children:
                offsprings.extend(child.get_offsprings())
            timeTagOfs = [int(item.tag) for n in offsprings for item in
                       n.get_data()]
            if 0 < len(timeTagOfs):
                if tagCurX > min(timeTagOfs):
                    return 1

        return 0

    def __cmp_tags(self, x):
        # 此处为什么不考虑当前数据的tag
        # or what about there is no datum in this node?, the answer is no.
        # current data is added into this node already
        # testFile = open('./testfile.dot', 'w')
        # self.tssb.print_graph(testFile)
        # testFile.close()

        datums = self.get_data()
        timeTag = [int(datum.tag) for datum in datums]
        if int(x[0].tag) < max(timeTag):
            return -1
        elif int(x[0].tag) > min(timeTag):
            return 1
        else:
            return 0


    def __is_good_gaps(self, x):
        lowerNode, upperNode = self.__find_neighbor_datum_n(x)

        lFlag = True
        uFlag = True

        if lowerNode is not None:
            lFlag = self.__is_good_gap(lowerNode, x, "lower")
        else:
            lFlag = True
        if upperNode is not None:
            uFlag = self.__is_good_gap(lowerNode, x, "upper")
        else:
            uFlag = True

        return lFlag and uFlag

    def data_log_likelihood(self, alleleConfig, baseline, maxCopyNumber):
        return self.complete_logprob(alleleConfig, baseline, maxCopyNumber)

    def complete_logprob(self, alleleConfig, baseline, maxCopyNumber):
        return sum([self.logprob([data], alleleConfig, baseline, maxCopyNumber)
                    for data in self.get_data()])

    def __find_neighbor_datum_n(self, x):
        # 对当前node中的stripe进行排序，计算gap
        datums = self.get_data()
        # it seems useless here
        if x[0] not in datums:
            datums.append(x[0])

        # here, it is possible that this node contains this data only
        # so return None

        if len(datums) < 1:
            raise Exception('Datum number error!')
        elif len(datums) == 1:
            return (None, None)
        else:
            datumsSortedL = sorted(datums,
                key=lambda item: 1.0*item.tReadNum/item.nReadNum)
            idx = datumsSortedL.index(x)
            if 0 == idx:
                return (None, datumsSortedL[1])
            elif len(datumsSortedL) - 1 == idx:
                return (datumsSortedL[idx-1], None)
            else:
                return (datumsSortedL[idx-1], datumsSortedL[idx+1])

    def __is_good_gap(self, lowerNode, upperNode, position):
        varpi = constants.VARPI

        rdrLower = 1.0*lowerNode.tReadNum/lowerNode.nReadNum
        rdrUpper = 1.0*upperNode.tReadNum/upperNode.nReadNum
        L = rdrUpper / rdrLower

        if "lower" == position:
            cn = lowerNode.copyNumber
        elif "upper" == position:
            cn = upperNode.copyNumber - 1

        if cn < 0:
            return False
        else:
            return L >= varpi * (1.0 + (self.param /
                    (cn * self.param + 2 * (1 - self.param))))

#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: tssb.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-12-06 15:39:03
#       History: phylowgs
# =============================================================================
'''
import sys
from copy import deepcopy
from os.path import commonprefix
from time import *

from numpy import *
from numpy.random import *

import scipy.stats
from gwpy.segments import Segment, SegmentList
from phySCNAClonal.model.usegsampler.segsupportive import MultiRangeSampler
from phySCNAClonal.preprocess.utils import get_cn_allele_config
from util import betapdfln, boundbeta, logsumexp, sticks_to_edges


# from phySCNAClonal.model.printo import show_tree_structure3



class TSSB(object):
    minDpAlpha = 1.0
    maxDpAlpha = 50.0
    minDpGamma = 1.0
    maxDpGamma = 10.0
    minAlphaDecay = 0.05
    maxAlphaDecay = 0.80

    def __init__(self,
                 dpAlpha=1.0,
                 dpGamma=1.0,
                 rootNode=None,
                 data=None,
                 minDepth=0,
                 maxDepth=50,
                 alphaDecay=1.0,
                 baseline=0,
                 maxCopyNumber=6):

        self.maxCopyNumber = maxCopyNumber
        # 此处的baseline是指什么
        self.baseline = baseline
        self.alleleConfig = get_cn_allele_config(maxCopyNumber)

        if rootNode is None:
            raise Exception("Root node must be specified.")

        self.minDepth = minDepth
        self.maxDepth = maxDepth
        self.dpAlpha = dpAlpha
        self.dpGamma = dpGamma
        self.alphaDecay = alphaDecay
        self.data = data
        self.dataNum = 0 if data is None else len(
            data)  # data.shape[0] #shankar
        self.root = {
            'node': rootNode,
            'main': boundbeta(1.0, dpAlpha) if self.minDepth == 0 else 0.0,
            'sticks': empty((0, 1)),
            'children': [],
            'tag': False
        }
        rootNode.tssb = self

        if False:
            dataU = rand(self.dataNum)
            self.assignments = []
            for n in range(self.dataNum):
                (c, path) = self.find_node(dataU[n])
                c.add_datum(n)
                self.assignments.append(c)
        else:
            self.assignments = []
            for n in range(self.dataNum):
                self.root['node'].add_datum(n)
                self.assignments.append(self.root['node'])

    def add_data(self, data):
        (weights, nodes) = self.get_mixture()
        newDataNum = len(data)  # data.shape[0] #shankar
        for n in range(newDataNum):
            logprobs = []
            for k, node in enumerate(nodes):
                logprobs.append(log(weights[k]) + node.logprob(data[n],
                                                               self.alleleConfig,
                                                               self.baseline,
                                                               self.maxCopyNumber))
            logprobs = array(logprobs)
            probs = exp(logprobs - logsumexp(logprobs))
            best_k = sum(rand() > cumsum(probs))
            nodes[best_k].add_datum(n + self.dataNum)
            self.assignments.append(nodes[best_k])
        self.data = vstack([self.data, data])
        self.dataNum += newDataNum

        # shankar

    #    def clear_data(self):
    #        dims = self.data.shape[1]
    #        for n in range(self.dataNum):
    #            self.assignments[n].remove_datum(n)
    #        self.assignments = []
    #        self.data        = empty((0,dims))
    #        self.dataNum    = 0

    def resample_node_params(self, iters=1):
        for iter in range(iters):

            def descend(root):
                for index, child in enumerate(root['children']):
                    descend(child)
                root['node'].resample_params()

            descend(self.root)

    def resample_assignments(self, tag):
        def path_lt(path1, path2):
            if len(path1) == 0 and len(path2) == 0:
                return 0
            elif len(path1) == 0:
                return 1
            elif len(path2) == 0:
                return -1
            s1 = "".join(map(lambda i: "%03d" % (i), path1))
            s2 = "".join(map(lambda i: "%03d" % (i), path2))

            return cmp(s2, s1)

        epsilon = finfo(float64).eps
        # this is not useful
        lengths = []

        # change data range
        # 需要保存索引号码

        tagL = array([int(item.tag) for item in self.data])
        timeTagSeq = sort(unique(tagL))
        currentTimeTagIdx = where(timeTagSeq == tag)[0][0]
        uSegL = MultiRangeSampler(0,1)
        if currentTimeTagIdx == 0:
            self.reset_time_tag()
        else:
            self.mark_negative_space(timeTagSeq[currentTimeTagIdx - 1])
            uNegtive = self.get_u_segL()
            uSegL.remove(uNegtive)


        for n in where(tagL == tag)[0]:
            # change to segment operation
            tempUSegL = deepcopy(uSegL)
            minU = tempUSegL.lowerBoundary
            maxU = tempUSegL.upperBoundary

            llhMapD = {}
            # Get an initial uniform variate.
            ancestors = self.assignments[n].get_ancestors()
            current = self.root
            indices = []
            for anc in ancestors[1:]:
                index = map(lambda c: c['node'],
                            current['children']).index(anc)
                current = current['children'][index]
                indices.append(index)

            # Here indices could be used as path to locate v stick

            currentVStick = self._locate_v_stick(indices)
            isInNegativeSpace = currentVStick['tag']

            llhS = -float('Inf')
            if isInNegativeSpace:
                oldLlh = -float('Inf')
                llhS = -float('Inf')
            else:
                oldLlh = self.assignments[n].logprob(self.data[n:n + 1],
                                                     self.alleleConfig,
                                                     self.baseline,
                                                     self.maxCopyNumber)
                llhMapD[self.assignments[n]] = oldLlh
                llhS = log(rand()) + oldLlh

            while True:
                newU = tempUSegL.sample()
                (newNode, newPath) = self.find_node(newU)
                if newNode.parent() is None:
                    # shankar: to make root node empty
                    newNode = newNode.children()[0]
                    newPath = [0]
                oldNode = self.assignments[n]
                oldNode.remove_datum(n)
                newNode.add_datum(n)
                self.assignments[n] = newNode
                if newNode in llhMapD:
                    newLlh = llhMapD[newNode]
                else:
                    ####################################
                    #  Record most likely copy number  #
                    ####################################
                    newLlh = newNode.logprob(self.data[n:n + 1],
                                             self.alleleConfig, self.baseline,
                                             self.maxCopyNumber)
                    if -float('Inf') == newLlh:
                        print >> sys.stderr, "Slice sampler weirdness"
                    llhMapD[newNode] = newLlh
                if newLlh > llhS:
                    break
                elif abs(maxU - minU) < epsilon:
                    if -float('Inf') == newLlh:
                        raise Exception("Slice sampler weirdness.")

                    newNode.remove_datum(n)
                    oldNode.add_datum(n)
                    self.assignments[n] = oldNode
                    print >> sys.stderr, "Slice sampler shrank down.  Keep current state."
                    break
                else:
                    newNode.remove_datum(n)
                    oldNode.add_datum(n)
                    self.assignments[n] = oldNode
                    pathComp = path_lt(indices, newPath)
                    if pathComp < 0:
                        tempUSegL.removeLeft(newU)
                    elif pathComp >= 0:  # temporary fix only!!!!!!
                        tempUSegL.removeRight(newU)
                    else:
                        raise Exception("Slice sampler weirdness.")
                    minU = tempUSegL.lowerBoundary
                    maxU = tempUSegL.upperBoundary
            lengths.append(len(newPath))
        lengths = array(lengths)

    def resample_assignments_crossing(self, timeOrderL, negativeSD, phiDL, srtree):
        def path_lt(path1, path2):
            if len(path1) == 0 and len(path2) == 0:
                return 0
            elif len(path1) == 0:
                return 1
            elif len(path2) == 0:
                return -1
            s1 = "".join(map(lambda i: "%03d" % (i), path1))
            s2 = "".join(map(lambda i: "%03d" % (i), path2))

            return cmp(s2, s1)

        epsilon = finfo(float64).eps
        # this is not useful
        lengths = []

        # change data range

        # 经过crossing 变换的时序数据的抽样方法如下：
        # 根据时序标记进行排序，先对未标时序的数据进行抽样。
        # 然后根据每一个规则内的时序进行分批次抽样。
        # 不同批次之间没有负空间依赖。
        # 在同一个批次内，当前时序的负空间包括所有已经抽样过得数据中小于当前时序
        # 的数据。
        # 每次进行节点搜索时，需要返回目标节点和所有祖先节点对应的R空间，以及当前节点所在的\varphi
        # 说对应的R。
        # 如果当前节点满足抽样条件，则要更新历史树对应的phi，如果有节点依赖当前节点的负空间，则要
        # 生成目标负空间。
        # 如果出现两个相同的phi值，需要进行当前节点内判断，如果已经有节点抽样过且其在所有的样本中
        # 的phi值都相同，那么不进行更新负空间。
        # 对上次抽样获得的节点需要对其是否在负空间中进行判断，如果在返回负无穷
        # 如果不在，但不满足历史限制条件，同样返回负无穷，但对其varphi所在R进行限制。


        #################################
        #  Sort data list in time order  #
        #################################

        nonOrderedDataIdLL = self._get_non_ordered_data(phiDL)
        totalTagedDataIdLL = timeOrderL + nonOrderedDataIdLL

        for timeOrderIndex, dataL in enumerate(totalTagedDataIdLL):
            if len(dataL) == 0:
                continue

            isNonOrderedData = False
            if dataL[0][1] == -1:
                isNonOrderedData = True

            uSegL = MultiRangeSampler(0,1)
            lastUL = MultiRangeSampler(0,0)
            lastt = -1

            for n, t in dataL:
                # change to segment operation

                tempUSegL = deepcopy(uSegL)
                uNegtive = None
                if not isNonOrderedData:
                    self.mark_specific_time_tag(negativeSD[n])
                    uNegtive = self.get_u_segL()
                    tempUSegL.remove(uNegtive)
                else:
                    self.mark_specific_time_tag(set([]))

                minU = tempUSegL.lowerBoundary
                maxU = tempUSegL.upperBoundary

                llhMapD = {}
                # Get an initial uniform variate.
                ancestors = self.assignments[n].get_ancestors()
                current = self.root
                indices = []
                for anc in ancestors[1:]:
                    index = map(lambda c: c['node'],
                                current['children']).index(anc)
                    current = current['children'][index]
                    indices.append(index)

                # Here indices could be used as path to locate v stick
                currentVStick = self._locate_v_stick(indices)
                isInNegativeSpace = currentVStick['tag']

                isSummingOK = True
                if isInNegativeSpace:
                    isSummingOK = False
                elif not isNonOrderedData:
                    currentPath = "".join(map(lambda i: "%03d" % (i), indices))
                    isSummingOK = srtree.is_good_path(timeOrderIndex, n, currentPath)

                llhS = -float('Inf')
                if isInNegativeSpace or not isSummingOK:
                    oldLlh = -float('Inf')
                    llhS = -float('Inf')
                else:
                    oldLlh = self.assignments[n].logprob(self.data[n:n + 1],
                                                        self.alleleConfig,
                                                        self.baseline,
                                                        self.maxCopyNumber)
                    llhMapD[self.assignments[n]] = oldLlh
                    llhS = log(rand()) + oldLlh

                timesRMVarphiR = 0
                lastVarphiRIndex = 0
                lastVarphiRL = srtree.get_lastVarphiRL(timeOrderIndex, n)
                while True:
                    newU = tempUSegL.sample()

                    (newNode, newPath, varphiR)  = self.find_node_varphi_range(newU)
                    newPathS = "".join(map(lambda i: "%03d" % (i), newPath))

                    if not isNonOrderedData and not srtree.is_good_path(
                        timeOrderIndex, n, newPathS):
                        if tempUSegL.remove(varphiR) is False or timesRMVarphiR > 500:
                            if lastVarphiRIndex >= len(lastVarphiRL):
                                raise Exception("Varphi range error!")
                            tempUSegL.assign_supportive(lastVarphiRL[lastVarphiRIndex])
                            tempUSegL.remove(uNegtive)
                            lastVarphiRIndex = lastVarphiRIndex + 1
                        timesRMVarphiR = timesRMVarphiR + 1
                        continue

                    if newNode.parent() is None:
                        # shankar: to make root node empty
                        newNode = newNode.children()[0]
                        newPath = [0]

                    oldNode = self.assignments[n]
                    oldNode.remove_datum(n)
                    newNode.add_datum(n)
                    self.assignments[n] = newNode
                    if newNode in llhMapD:
                        newLlh = llhMapD[newNode]
                    else:
                        ####################################
                        #  Record most likely copy number  #
                        ####################################
                        newLlh = newNode.logprob(self.data[n:n + 1],
                                                self.alleleConfig, self.baseline,
                                                self.maxCopyNumber)
                        if -float('Inf') == newLlh:
                            print >> sys.stderr, "Slice sampler weirdness"
                        llhMapD[newNode] = newLlh
                    if newLlh > llhS:
                        break
                    elif abs(maxU - minU) < epsilon:
                        if -float('Inf') == newLlh:
                            raise Exception("Slice sampler weirdness.")

                        newNode.remove_datum(n)
                        oldNode.add_datum(n)
                        self.assignments[n] = oldNode
                        print >> sys.stderr, "Slice sampler shrank down.  Keep current state."
                        break
                    else:
                        newNode.remove_datum(n)
                        oldNode.add_datum(n)
                        self.assignments[n] = oldNode
                        pathComp = path_lt(indices, newPath)
                        if pathComp < 0:
                            tempUSegL.removeLeft(newU)
                        elif pathComp >= 0:  # temporary fix only!!!!!!
                            tempUSegL.removeRight(newU)
                        else:
                            raise Exception("Slice sampler weirdness.")
                        minU = tempUSegL.lowerBoundary
                        maxU = tempUSegL.upperBoundary

                if not isNonOrderedData:
                    srtree.update_path_varphiR(timeOrderIndex, n, newPathS, varphiR)
                lengths.append(len(newPath))
        lengths = array(lengths)

    def resample_assignments_singlecell(self, orderMatrix, srtree):
        def path_lt(path1, path2):
            if len(path1) == 0 and len(path2) == 0:
                return 0
            elif len(path1) == 0:
                return 1
            elif len(path2) == 0:
                return -1
            s1 = "".join(map(lambda i: "%03d" % (i), path1))
            s2 = "".join(map(lambda i: "%03d" % (i), path2))

            return cmp(s2, s1)

        epsilon = finfo(float64).eps
        # this is not useful
        lengths = []

        # filter out data not in orderMatrix
        scData = set(orderMatrix[0,])
        nonOrderedData = set(range(len(self.data))) - orderedData

        # 对单细胞测序数据中每一个样本，即矩阵中每一行路径搜索

        # 每一个单细胞测序样本路径搜索完毕后需要记录在已经搜索的内容中
        # 用字典类型表示
        # 此处每一个搜索到的节点可以放到节点类中保存，用来方便的更新其祖先节点的
        # {dataIdx: Node}

        scDataFoundD = {}

        lastDataIdx = -1
        for orderVector in orderMatrix[1:,]:

            # 根据当前单细胞测序样本中的变异Stage进行其中包含-1
            currentSampleStageS = sorted(set(orderVector))

            # 此处应该设置为树根节点的varphiR - piR
            # 该变量的初始化应该放在路径搜索的起点
            # 每一个路径中有多个Stage
            lastStageRemainR = MultiRangeSampler(0,1)

            # 记录当前路径信息
            currentStageStatus = {}
            currentStageStatus["lowest_epsilon"] = ""
            currentStageStatus["lowest_remain_r"] = MultiRangeSampler(0,1)
            currentStageStatus["lowest_idx"] = -1

            for stage in currentSampleStageS:
                if stage == -1
                    # 当前单细胞测序中不含有该变异
                    continue
                # 获取当前stage的数据id
                idxL = orderMatrix[0, where(orderVector==stage)[0]]

                lastStageLowestEpsilon = ""

                for n in idxL:
                    if n in scDataFoundD.keys():
                        # 此处需要搜索最下方的节点
                        # 一个stage
                        if len(lastStageLowestEpsilon) < len(scDataFoundD[n].epsilon):
                            lastStageLowestEpsilon = scDataFoundD[n].epsilon
                            lastStageRemainR = deepcopy(scDataFoundD[n].remainR)
                            currentStageStatus["lowest_remain_r"] = deepcopy(scDataFoundD[n].remainR)
                            # 跟据最下方节点的索引，搜索树，获得目标空间
                            lastDataIdx = n
                        continue
                    else:
                        # 需要搜索当前stage 的搜索空间
                        #
                        # 第一步，确定当前的搜索空间，LastStage中　所有 piR的交集
                        # 即，位于最下方的节点的后代
                        # 第二步，获取当前stage中所有的祖先节点的R
                        # 第三步，获取所有的piR的交集与所有祖先节点的交集
                        # 获取当前状态位于最下方节点的piR，与上面的交集取并集
                        tempUSegL = deepcopy(lastStageRemainR)

                        if currentStageStatus["lowest_idx"] != -1:
                            tempcss = deepcopy(lastStageRemainR)
                            self.mark_specific_time_tag([currentStageStatus["lowest_idx"]])
                            ancestorsR = self.get_u_segL()
                            tempcss.remove(ancestorsR)
                            tempcss.remove(currentStageStatus["lowest_remain_r"])
                            tempUSegL.remove(tempcss)

                        minU = tempUSegL.lowerBoundary
                        maxU = tempUSegL.upperBoundary

                        llhMapD = {}
                        # Get an initial uniform variate.
                        ancestors = self.assignments[n].get_ancestors()
                        current = self.root
                        indices = []
                        for anc in ancestors[1:]:
                            index = map(lambda c: c['node'],
                                        current['children']).index(anc)
                            current = current['children'][index]
                            indices.append(index)

                        targetEpsilon = "".join(map(lambda i: "%03d" % (i), indices))
                        # Here indices could be used as path to locate v stick
                        # 此处需要判断上次抽样的节点是否落入当前path中
                        # 如果没有落入当前path中，赋予最小概率
                        # 判断方法应该使用epsilon 判断
                        isInPath = self.__is_in_path(lastStageLowestEpsilon,
                                                     currentStageStatus["lowest_epsilon"],
                                                     targetEpsilon)
                        llhS = -float('Inf')
                        if not isInPath:
                            oldLlh = -float('Inf')
                        else:
                            oldLlh = self.assignments[n].logprob(self.data[n:n + 1],
                                                                self.alleleConfig,
                                                                self.baseline,
                                                                self.maxCopyNumber)
                            llhMapD[self.assignments[n]] = oldLlh
                            llhS = log(rand()) + oldLlh

                        while True:
                            newU = tempUSegL.sample()
                            (tnode, newNode, newPath, varphiR, piR)  = self.find_node_varphi_pi_range(newU)
                            newPathEpsilon = "".join(map(lambda i: "%03d" % (i), newPath))

                            if newNode.parent() is None:
                                # shankar: to make root node empty
                                newNode = newNode.children()[0]
                                newPath = [0]

                            oldNode = self.assignments[n]
                            oldNode.remove_datum(n)
                            newNode.add_datum(n)
                            self.assignments[n] = newNode
                            if newNode in llhMapD:
                                newLlh = llhMapD[newNode]
                            else:
                                ####################################
                                #  Record most likely copy number  #
                                ####################################
                                newLlh = newNode.logprob(self.data[n:n + 1],
                                                        self.alleleConfig, self.baseline,
                                                        self.maxCopyNumber)
                                if -float('Inf') == newLlh:
                                    print >> sys.stderr, "Slice sampler weirdness"
                                llhMapD[newNode] = newLlh
                            if newLlh > llhS:
                                break
                            elif abs(maxU - minU) < epsilon:
                                if -float('Inf') == newLlh:
                                    raise Exception("Slice sampler weirdness.")

                                newNode.remove_datum(n)
                                oldNode.add_datum(n)
                                self.assignments[n] = oldNode
                                print >> sys.stderr, "Slice sampler shrank down.  Keep current state."
                                break
                            else:
                                newNode.remove_datum(n)
                                oldNode.add_datum(n)
                                self.assignments[n] = oldNode
                                pathComp = path_lt(indices, newPath)
                                if pathComp < 0:
                                    tempUSegL.removeLeft(newU)
                                elif pathComp >= 0:  # temporary fix only!!!!!!
                                    tempUSegL.removeRight(newU)
                                else:
                                    raise Exception("Slice sampler weirdness.")
                                minU = tempUSegL.lowerBoundary
                                maxU = tempUSegL.upperBoundary

                        # 此处需要更新当前状态信息，和所有已经抽样的数据的状态信
                        # 息
                        if len(currentStageStatus["lowest_epsilon"]) < len(newPathEpsilon):
                            currentStageStatus["lowest_idx"] = n
                            currentStageStatus["lowest_epsilon"] = newPathEpsilon
                            currentStageStatus["lowest_remain_r"] = varphiR.remove(piR)

                        # 此处需要更新当前节点的remain r
                        # 由于有空节点的存在，需要对树结构进行迭代
                        # 此处对全树搜索开销太大
                        # 只对当前路径最新抽样且最为低节点进行更新。
                        # 更新方法是通过对当前路径进行逆向迭代搜索更新remain r
                        # 输入数据：当前路径已经搜索到的所有节点

                        self.assignments[n].varphiR = varphiR
                        self.assignments[n].piR = piR
                        # 此处需要判断是否是最低节点，
                        # 如果不是最低节点，需要更新该节点的
                        lastDataIdx = n

                        lengths.append(len(newPath))
                lengths = array(lengths)

    def _find_path_init_R(self, targetNode, varphiR, piR, negativeDataS):
        # 此处判断当前节点中是否含有目标数据的id，如果含有该数据
        # 则对当前数据进行抽样

        def descend(rootNode, negativeDataS):
            if len(set.intersection(rootNode.data(), negativeDataS)) > 0:
                return True
            else:
                reFlag = False
                for childNode in rootNode.children():
                    reFlag = reFlag or descend(childNode, negativeDataS)
                    return reFlag
                return False

        reRange = deepcopy(varphiR)
        reRange.remove(piR)
        for childNode in targetNode['node'].children():
            if descend(childNode):
                reRange.remove(childNode.varphiR)

        return reRange


    def __is_in_path(self, lastStageLowestEpsilon, currentStageLoestEpsilon,
                     targetEpsilon):
        ltcp = commonprefix([lastStageLowestEpsilon, targetEpsilon])
        ctcp = commonprefix([currentStageLoestEpsilon, targetEpsilon])
        if ltcp == lastStageLowestEpsilon:
            if ctcp == currentStageLoestEpsilon:
                return True
            elif ctcp == targetEpsilon:
                return True
            else:
                return False
        else:
            retrun False


    def _get_non_ordered_data(self, phiDL):
        orderedData = set([j for phiD in phiDL for j in phiD.keys()])
        nonOrderedData = set(range(len(self.data))) - orderedData
        nonOrderedDataIdL = [(item, -1) for item in nonOrderedData]
        return [nonOrderedDataIdL]


    def _locate_v_stick(self, indices):
        """locate and return v stick by indices

        :indices: The indices to current data node
        :returns: v stick

        """
        target = self.root
        for index in indices:
            target = target['children'][index]

        return target

    def cull_tree(self):
        def descend(root):
            counts = array(map(lambda child: descend(child), root['children']))
            keep = len(trim_zeros(counts, 'b'))

            for child in root['children'][keep:]:
                child['node'].kill()
                del child['node']

            root['sticks'] = root['sticks'][:keep]
            root['children'] = root['children'][:keep]

            return sum(counts) + root['node'].num_local_data()

        descend(self.root)

    def resample_sticks(self):
        def descend(root, depth=0):

            dataDown = 0
            indices = range(len(root['children']))
            indices.reverse()
            for i in indices:
                child = root['children'][i]
                childData = descend(child, depth + 1)
                postAlpha = 1.0 + childData
                postBeta = self.dpGamma + dataDown
                root['sticks'][i] = boundbeta(
                    postAlpha,
                    postBeta) if depth != 0 else .999999  # shankar
                dataDown += childData

            # Resample the main break.
            dataHere = root['node'].num_local_data()
            postAlpha = 1.0 + dataHere
            postBeta = (self.alphaDecay**depth) * self.dpAlpha + dataDown
            root['main'] = boundbeta(
                postAlpha, postBeta) if self.minDepth <= depth else 0.0
            if depth == 0:
                root['main'] = 1e-30  # to make root node empty (shankar)

            return dataHere + dataDown

        descend(self.root)

    def resample_stick_orders(self):
        def descend(root, depth=0):
            if not root['children']:
                return

            newOrder = []
            represented = set(
                filter(lambda i: root['children'][i]['node'].has_data(),
                       range(len(root['children']))))
            # 每一phi stick 的实际长度
            allWeights = diff(hstack([0.0, sticks_to_edges(root['sticks'])]))
            while True:
                if not represented:
                    break

                u = rand()
                while True:
                    subIndices = filter(lambda i: i not in newOrder,
                                         range(root['sticks'].shape[0]))
                    # 此处添加了剩余空间的长度
                    subWeights = hstack(
                        [allWeights[subIndices], 1.0 - sum(allWeights)])
                    # 每一个空间所占用的比重
                    subWeights = subWeights / sum(subWeights)
                    # 随机获得一个位置，此位置之前完整空间的个数
                    index = sum(u > cumsum(subWeights))

                    if index == len(subIndices):
                        root['sticks'] = vstack(
                            [root['sticks'],
                             boundbeta(1, self.dpGamma)])
                        root['children'].append({
                            'node': root['node'].spawn(),
                            'main': boundbeta(1.0, (self.alphaDecay ** (depth + 1)) * self.dpAlpha) if self.minDepth <= (depth + 1) else 0.0, # 此处minDepth 应该是手动控制的
                            'sticks': empty((0, 1)),
                            'children': [],
                            'tag': False
                        })
                        allWeights = diff(
                            hstack([0.0, sticks_to_edges(root['sticks'])]))
                    else:
                        index = subIndices[index]
                        break
                newOrder.append(index)
                represented.discard(index)

            newChildren = []
            for k in newOrder:
                child = root['children'][k]
                newChildren.append(child)
                descend(child, depth + 1)

            for k in filter(lambda k: k not in newOrder,
                            range(root['sticks'].shape[0])):
                root['children'][k]['node'].kill()
                del root['children'][k]['node']

            root['children'] = newChildren
            root['sticks'] = zeros((len(root['children']), 1))

        descend(self.root)

        # Immediately resample sticks.
        self.resample_sticks()

    def resample_hypers(self, dpAlpha=True, alphaDecay=True, dpGamma=True):
        def dp_alpha_llh(dpAlpha, alphaDecay):
            def descend(dpAlpha, root, depth=0):
                llh = betapdfln(root['main'], 1.0, (alphaDecay**depth) *
                                dpAlpha) if self.minDepth <= depth else 0.0
                for child in root['children']:
                    llh += descend(dpAlpha, child, depth + 1)
                return llh

            return descend(dpAlpha, self.root)

        if dpAlpha:
            upper = self.maxDpAlpha
            lower = self.minDpAlpha
            llhS = log(rand()) + dp_alpha_llh(self.dpAlpha, self.alphaDecay)
            while True:
                newDpAlpha = (upper - lower) * rand() + lower
                newLlh = dp_alpha_llh(newDpAlpha, self.alphaDecay)
                if newLlh > llhS:
                    break
                elif newDpAlpha < self.dpAlpha:
                    lower = newDpAlpha
                elif newDpAlpha > self.dpAlpha:
                    upper = newDpAlpha
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.dpAlpha = newDpAlpha

        if alphaDecay:
            upper = self.maxAlphaDecay
            lower = self.minAlphaDecay
            llhS = log(rand()) + dp_alpha_llh(self.dpAlpha, self.alphaDecay)
            while True:
                newAlphaDecay = (upper - lower) * rand() + lower
                newLlh = dp_alpha_llh(self.dpAlpha, newAlphaDecay)
                if newLlh > llhS:
                    break
                elif newAlphaDecay < self.alphaDecay:
                    lower = newAlphaDecay
                elif newAlphaDecay > self.alphaDecay:
                    upper = newAlphaDecay
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.alphaDecay = newAlphaDecay

        def dp_gamma_llh(dpGamma):
            def descend(dpGamma, root):
                llh = 0
                for i, child in enumerate(root['children']):
                    llh += betapdfln(root['sticks'][i], 1.0, dpGamma)
                    llh += descend(dpGamma, child)
                return llh

            return descend(dpGamma, self.root)

        if dpGamma:
            upper = self.maxDpGamma
            lower = self.minDpGamma
            llhS = log(rand()) + dp_gamma_llh(self.dpGamma)
            while True:
                newDpGamma = (upper - lower) * rand() + lower
                newLlh = dp_gamma_llh(newDpGamma)
                if newLlh > llhS:
                    break
                elif newDpGamma < self.dpGamma:
                    lower = newDpGamma
                elif newDpGamma > self.dpGamma:
                    upper = newDpGamma
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.dpGamma = newDpGamma

    def draw_data(self, dataNum=1, **args):
        self.data = []
        self.assignments = []
        for n in range(dataNum):
            u = rand()
            (node, path) = self.find_node(u)
            self.data.append(node.sample(args))
            self.assignments.append(node)
            node.add_datum(n)
            self.dataNum += 1
        self.data = concatenate(self.data)
        return self.data

    def resample_data(self, **args):
        for n in range(self.dataNum):
            u = rand()
            (node, path) = self.find_node(u)
            self.assignments[n].remove_datum(n)
            node.add_datum(n)
            self.assignments[n] = node
            self.data[n] = node.sample(args)[0]

    def mark_negative_space(self, tag):
        return self.mark_time_tag2(tag)

    def mark_specific_time_tag(self, Q):
        """mark negative node

        :Q: TODO
        :returns: TODO

        """
        def descend(root):
            if 0 == len(root['children']):
                for item in root['node'].data:
                    if item in Q:
                        root['tag'] = True
                        return True
                root['tag'] = False
                return False
            else:
                if 0 < sum([descend(child) for child in root['children']]):
                    root['tag'] = True
                    return True
                else:
                    for item in root['node'].data:
                        if item in Q:
                            root['tag'] = True
                            return True
                    root['tag'] = False
                    return False

        descend(self.root)

    def mark_time_tag2(self, tag):
        """
        mark each node's time tag status.
        All the nodes in the tssb tree, that contains offspings data with tag <=
        given tag.
        """
        def descend(root, tag):
            if 0 == len(root['children']):
                timeTags = [int(item.tag) for item in root['node'].get_data() if
                            int(item.tag) <= tag]
                if 0 < len(timeTags):
                    root['tag'] = True
                    return True
                else:
                    root['tag'] = False
                    return False
            else:
                if 0 < sum([descend(child, tag) for child in root['children']]):
                    root['tag'] = True
                    return True
                else:
                    timeTags = [int(item.tag) for item in root['node'].get_data() if
                                int(item.tag) <= tag]
                    if 0 < len(timeTags):
                        root['tag'] = True
                        return True
                    else:
                        root['tag'] = False
                        return False

        descend(self.root, tag)

    def mark_time_tag(self, tag):
        """
        mark each node's time tag status.
        All the nodes in the tssb tree, that contains offspings data with tag <=
        given tag.
        """
        def descend(root, tag):
            if 0 == len(root['children']):
                root['tag'] = False

                timeTags = [int(item.tag) for item in root['node'].get_data() if
                            int(item.tag) <= tag]
                if 0 < len(timeTags):
                    return True
                else:
                    return False
            else:
                if 0 < sum([descend(child, tag) for child in root['children']]):
                    root['tag'] = True
                    return True
                else:
                    root['tag'] = False
                    timeTags = [int(item.tag) for item in root['node'].get_data() if
                                int(item.tag) <= tag]
                    if 0 < len(timeTags):
                        return True
                    else:
                        return False

        descend(self.root, tag)

    def reset_time_tag(self):
        def descend(root):
            root['tag'] = False
            for child in root['children']:
                descend(child)

        descend(self.root)

    def find_node(self, u):
        def descend(root, u, depth=0):
            if depth >= self.maxDepth:
                # print >>sys.stderr, "WARNING: Reached maximum depth."
                return (root['node'], [])
            elif u < root['main']:
                if root['tag']:
                     print >>sys.stderr, "Negative space located!!."

                return (root['node'], [])
            else:
                # Rescale the uniform variate to the remaining interval.
                u = (u - root['main']) / (1.0 - root['main'])

                # Perhaps break sticks out appropriately.
                if depth > 0:
                    while not root['children'] or (
                            1.0 - prod(1.0 - root['sticks'])) < u:
                        root['sticks'] = vstack([
                            root['sticks'],
                            # 注意此处为右边界
                            boundbeta(1, self.dpGamma) if depth != 0 else .999
                        ])  # shankar
                        root['children'].append({
                            'node':
                            root['node'].spawn(),
                            'main':
                            boundbeta(1.0, (self.alphaDecay**
                                            (depth + 1)) * self.dpAlpha)
                            if self.minDepth <= (depth + 1) else 0.0,
                            'sticks':
                            empty((0, 1)),
                            'children': [],
                            'tag': False
                        })

                    edges = 1.0 - cumprod(1.0 - root['sticks'])
                    index = sum(u > edges)
                    edges = hstack([0.0, edges])
                    u = (u - edges[index]) / (edges[index + 1] - edges[index])

                    (node, path) = descend(root['children'][index], u,
                                           depth + 1)
                else:
                    index = 0
                    (node, path) = descend(root['children'][index], u,
                                           depth + 1)

                path.insert(0, index)

                return (node, path)

        return descend(self.root, u)

    def find_node_varphi_pi_range(self, u):
        # 此处需要返回当前抽样节点的祖先节点的对应的
        def descend(root, u, varphiR, depth=0):
            if depth >= self.maxDepth:
                # print >>sys.stderr, "WARNING: Reached maximum depth."
                return (root, [], varphiR, [varphiR[0], piEnd])
            elif u < root['main']:
                if root['tag']:
                     print >>sys.stderr, "Negative space located!!."
                return (root, [], varphiR, [varphiR[0], piEnd])
            else:
                # Rescale the uniform variate to the remaining interval.
                u = (u - root['main']) / (1.0 - root['main'])
                varphiR[0] = (varphiR[1] - varphiR[0]) * root['main'] + varphiR[0]

                # Perhaps break sticks out appropriately.
                if depth > 0:
                    while not root['children'] or (
                            1.0 - prod(1.0 - root['sticks'])) < u:
                        root['sticks'] = vstack([
                            root['sticks'],
                            # 注意此处为右边界
                            boundbeta(1, self.dpGamma) if depth != 0 else .999
                        ])  # shankar
                        root['children'].append({
                            'node':
                            root['node'].spawn(),
                            'main':
                            boundbeta(1.0, (self.alphaDecay**
                                            (depth + 1)) * self.dpAlpha)
                            if self.minDepth <= (depth + 1) else 0.0,
                            'sticks':
                            empty((0, 1)),
                            'children': [],
                            'tag': False
                        })

                    edges = 1.0 - cumprod(1.0 - root['sticks'])
                    index = sum(u > edges)
                    edges = hstack([0.0, edges])
                    u = (u - edges[index]) / (edges[index + 1] - edges[index])

                    varphiR = [
                        (varphiR[1]-varphiR[0])*edges[index]+varphiR[0],
                        (varphiR[1]-varphiR[0])*edges[index+1]+varphiR[0]]

                    (tnode, path, varphiR, piR) = descend(root['children'][index], u,
                                                    varphiR, depth + 1)
                else:
                    index = 0
                    (tnode, path, varphiR, piR) = descend(root['children'][index], u,
                                                    varphiR, depth + 1)

                path.insert(0, index)

                return (tnode, path, varphiR, piR)

        tn, p, vR, piR = descend(self.root, u, [0, 1])
        return tn, tn['node'], p, SegmentList([Segment(vR[0], vR[1])]), SegmentList([Segment(piR[0], piR[1])])


    def find_node_varphi_range(self, u):
        def descend(root, u, varphiR, depth=0):
            if depth >= self.maxDepth:
                # print >>sys.stderr, "WARNING: Reached maximum depth."
                return (root['node'], [], varphiR)
            elif u < root['main']:
                if root['tag']:
                     print >>sys.stderr, "Negative space located!!."
                return (root['node'], [], varphiR)
            else:
                # Rescale the uniform variate to the remaining interval.
                u = (u - root['main']) / (1.0 - root['main'])
                varphiR[0] = (varphiR[1] - varphiR[0]) * root['main'] + varphiR[0]

                # Perhaps break sticks out appropriately.
                if depth > 0:
                    while not root['children'] or (
                            1.0 - prod(1.0 - root['sticks'])) < u:
                        root['sticks'] = vstack([
                            root['sticks'],
                            # 注意此处为右边界
                            boundbeta(1, self.dpGamma) if depth != 0 else .999
                        ])  # shankar
                        root['children'].append({
                            'node':
                            root['node'].spawn(),
                            'main':
                            boundbeta(1.0, (self.alphaDecay**
                                            (depth + 1)) * self.dpAlpha)
                            if self.minDepth <= (depth + 1) else 0.0,
                            'sticks':
                            empty((0, 1)),
                            'children': [],
                            'tag': False
                        })

                    edges = 1.0 - cumprod(1.0 - root['sticks'])
                    index = sum(u > edges)
                    edges = hstack([0.0, edges])
                    u = (u - edges[index]) / (edges[index + 1] - edges[index])

                    varphiR = [
                        (varphiR[1]-varphiR[0])*edges[index]+varphiR[0],
                        (varphiR[1]-varphiR[0])*edges[index+1]+varphiR[0]]

                    (node, path, varphiR) = descend(root['children'][index], u,
                                                    varphiR, depth + 1)
                else:
                    index = 0
                    (node, path, varphiR) = descend(root['children'][index], u,
                                                    varphiR, depth + 1)

                path.insert(0, index)

                return (node, path, varphiR)

        n, p, vR = descend(self.root, u, [0, 1])
        return n, p, SegmentList([Segment(vR[0], vR[1])])

    def get_u_segL(self):
        """
        return the non-supportive range of u in the find_node function.
        here, all the target node is pre-marked in the mark_time_tag function.
        """
        def descend(root):
            if not root['tag']:
                return None, None
            elif 0 == len(root['children']):
                return array([0]), array(root['main'])
            else:
                # get edges
                starts = array([])
                ends = array([])

                edges = 1.0 - cumprod(1.0 - root['sticks'])
                edges = hstack([0.0, edges])
                for index, child in zip(range(len(root['children'])),
                                        root['children']):
                    startsChild, endsChild = descend(child)
                    if startsChild is None or endsChild is None:
                        continue
                    else:
                        startsChild = startsChild *\
                            (edges[index + 1] - edges[index]) + edges[index]
                        endsChild = endsChild *\
                            (edges[index + 1] - edges[index]) + edges[index]
                        starts = r_[starts, startsChild]
                        ends = r_[ends, endsChild]

                starts = starts * (1 - root['main']) + root['main']
                ends = ends * (1 - root['main']) + root['main']

                starts = r_[starts, 0]
                ends = r_[ends, root['main']]
                return starts, ends

        starts, ends = descend(self.root)
        sl = SegmentList([])

        if starts is not None:
            assert(len(starts) == len(ends))
            for i in range(len(starts)):
                sl.append(Segment(starts[i], ends[i]))
            sl.coalesce()
        else:
            assert(ends == None)

        return sl


    def get_nodes(self):
        def descend(root):
            node = [root['node']]
            for child in root['children']:
                child_nodes = descend(child)
                node.extend(child_nodes)
            return node

        return descend(self.root)

    def get_mixture(self):
        def descend(root, mass):
            weight = [mass * root['main']]
            node = [root['node']]
            edges = sticks_to_edges(root['sticks'])
            weights = diff(hstack([0.0, edges]))

            for i, child in enumerate(root['children']):
                (child_weights, child_nodes) = descend(
                    child, mass * (1.0 - root['main']) * weights[i])
                weight.extend(child_weights)
                node.extend(child_nodes)
            return (weight, node)

        # 返回两个向量
        return descend(self.root, 1.0)

    def complete_data_log_likelihood(self):
        weights, nodes = self.get_mixture()
        llhs = []
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.num_local_data() * log(weights[i]) +
                            node.data_log_likelihood(self.alleleConfig,
                                                     self.baseline,
                                                     self.maxCopyNumber))
        return sum(array(llhs))

    def complete_log_likelihood(self):
        weights, nodes = self.get_mixture()
        llhs = [
            self.dp_alpha_llh(self.dpAlpha, self.alphaDecay),
            self.dp_gamma_llh(self.dpGamma)
        ]
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.data_log_likelihood(self.alleleConfig,
                                                     self.baseline,
                                                     self.maxCopyNumber))
        return sum(array(llhs))

    def dp_alpha_llh(self, dpAlpha, alphaDecay):
        def descend(dpAlpha, root, depth=0):
            llh = betapdfln(root['main'], 1.0, (alphaDecay**depth) *
                            dpAlpha) if self.minDepth <= depth else 0.0
            for child in root['children']:
                llh += descend(dpAlpha, child, depth + 1)
            return llh

        return descend(dpAlpha, self.root)

    def dp_gamma_llh(self, dpGamma):
        def descend(dpGamma, root):
            llh = 0
            for i, child in enumerate(root['children']):
                llh += betapdfln(root['sticks'][i], 1.0, dpGamma)
                llh += descend(dpGamma, child)
            return llh

        return descend(dpGamma, self.root)

    def print_graph(self, fh, base_width=5000, min_width=5):
        print >> fh, """graph: { title:            "TSSB Graph"  \
                                portsharing:      no            \
                                smanhattanedges:  yes           \
                                equalydist:       yes           \
                                layout_algorithm: tree          \
                                node.fontname:    "helvR8"      \
                                node.height:      25 """


        print >> fh, """node: { label:"%0.5f" title:"%s" width:%d}""" % (
            self.root['main'], "X",
            max(int(self.root['main'] * base_width), min_width))

        def descend(root, name, mass):
            total = 0.0
            edges = sticks_to_edges(root['sticks'])
            weights = diff(hstack([0.0, edges]))
            for i, child in enumerate(root['children']):
                childName = "%s-%d" % (name, i)
                childMass = mass * weights[i] * child['main']
                print >> fh, """node: {  label:"%0.5f" title:"%s" width:%d}""" % (
                    childMass, childName,
                    max(int(childMass * base_width), min_width))
                print >> fh, """edge: { source:"%s" target:"%s" anchor:1}""" % (
                    name, childName)
                total += childMass + descend(child, childName,
                                              mass * weights[i] *
                                              (1.0 - child['main']))
            return total

        print >> fh, """}"""

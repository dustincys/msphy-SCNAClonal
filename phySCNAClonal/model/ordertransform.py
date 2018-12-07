#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: ordertransform.py
#          Desc: transform crossing rule and summing rule into time ordered rule
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-11-26 20:54:40
#       History:
# =============================================================================
'''
from collections import deque
from copy import deepcopy
from os.path import commonprefix

import numpy as np


class CS2T:

    """
    transform crossing rule to time ordered rule and keep the summing rule
    The configuration file in following format:

    # Current complete rule
    3|0.7   8|0.6   7|0.5   4|0.3   6|0.3   5|0.2   0|0.1   2|0.1   1|0.05
    # rule 1
    1|0.8   2|0.7   5|0.4   6|0.3   8|0.2
    5   8   6   1   2
    # rule 2
    3   4   7   0

    Each digit is the index j of variant.
    Each line mean the phi_j reverse order of each sample $\theVector{\varkappa}
    _{k,l}$.
    """

    def __init__(self, configFilePath):
        self._configFilePath = configFilePath

    def transform(self):
        varkappaMatrixArray, phiDL = self.__read_config_file()
        timeOrderL, negativeSD = self.__toTimeOrder(varkappaMatrixArray, phiDL)
        return timeOrderL, negativeSD, phiDL

    def __toTimeOrder(self, varkappaMatrixArray, phiDL):
        """ tranform to time order
        :returns: time order

        """
        negativeSD = {}
        totalTimeOrderL = []

        for ruleIdx, varkappaMatrix in enumerate(varkappaMatrixArray):
            currentVarkappa = np.array(varkappaMatrix[-1])
            currentData = [(idx, phiDL[ruleIdx][idx][-1]) for idx in currentVarkappa]

            for si in range(len(varkappaMatrix)-1):
                targetVarkappa = np.array(varkappaMatrix[-1])
                targetData = [(idx, phiDL[ruleIdx][idx][si]) for idx in targetVarkappa]

                for kj in range(len(currentVarkappa[:])):
                    currentSmallS = set([i for i, n in filter(lambda item: item[1] <= phiDL[ruleIdx][currentVarkappa[kj]][-1], currentData)])
                    targetSmallS = set([i for i, n in filter(lambda item: item[1] < phiDL[ruleIdx][currentVarkappa[kj]][si], targetData)])
                    negativeS = currentSmallS - targetSmallS

                    if currentVarkappa[kj] not in negativeSD.keys():
                        negativeSD[currentVarkappa[kj]] = negativeS
                    else:
                        negativeSD[currentVarkappa[kj]] = set.union(negativeS, negativeSD[currentVarkappa[kj]])

            tempTimeOrder = self.__idxOrder_to_timeOrder(negativeSD, currentData)
            for idx, time in tempTimeOrder:
                preTimeS = set([ i for i, t in filter(lambda item: item[1] < time, tempTimeOrder)])
                if idx not in negativeSD.keys():
                    negativeSD[idx] = preTimeS
                else:
                    negativeSD[idx] = set.union(negativeSD[idx], preTimeS)

            totalTimeOrderL.append(tempTimeOrder)

        return totalTimeOrderL, negativeSD

    def __idxOrder_to_timeOrder(self, negativeSD, currentData):
        tQueue = deque()
        tStack = deque()

        currentVarkappa = sorted(currentData, key=lambda item: item[1], reverse=True)

        nsL = [(idx, deepcopy(negativeSD[idx]), phi) if idx in
                negativeSD.keys() else (idx, set(), phi) for idx, phi in
                currentVarkappa]

        stackNegativeSet = set()

        for idx, ns, phi in nsL:
            if len(stackNegativeSet) == 0:
                if len(tStack) > 0:
                    tQueue.append(tStack)
                tStack = deque()
            tStack.append([idx, ns, phi])
            stackNegativeSet = set.union(ns, stackNegativeSet)
            if idx in stackNegativeSet:
                stackNegativeSet.remove(idx)
        if len(tStack) > 0:
            tQueue.append(tStack)

        timeOrder = []
        currentTime = -1

        while len(tQueue) != 0:
            tStack = tQueue.popleft()
            tStack = deque(sorted(tStack, key=lambda item: (len(item[1]), item[2]),
                                  reverse=True))
            lastTimeIdxS = set()
            lastPhi = -1
            while len(tStack) != 0:
                if len(tStack[-1][1]) != 0:
                    for item in tStack:
                        item[1] = item[1] - lastTimeIdxS
                    tStack = deque(sorted(tStack, key=lambda item: (len(item[1]), item[2]),
                                          reverse=True))
                    lastTimeIdxS.clear()
                    currentTime = currentTime + 1
                elif lastPhi != tStack[-1][2]:
                    currentTime = currentTime + 1
                else:
                    pass

                idx, nsl, phi = tStack.pop()
                timeOrder.append((idx, currentTime))
                lastPhi = phi
                lastTimeIdxS.add(idx)

        return timeOrder

    def __read_config_file(self):
        """Read config file: index|phi, phi descreasing order
        # rule 1
        1|0.8	2|0.7	5|0.4	6|0.3	8|0.2
        5|0.7	8|0.5	6|0.1	1|0.05	2|0.02
        8|0.6	6|0.3	5|0.2	2|0.1	1|0.05
        # rule 2
        3|0.8	4|0.5	7|0.3	0|0.2
        3|0.7	7|0.5	4|0.3	0|0.1
        """
        with open(self._configFilePath) as configFile:
            nextNewVarkappa = False
            nextCurrent = False

            varkappaMatrix = []
            varkappaMatrixArray = []
            phiDL = []
            phiD = {}
            currentVarkappa = []
            currentPhiD = {}

            for line in configFile:
                line = line.strip().lower()
                if line == "":
                    continue
                if line.startswith("#"):
                    nextNewVarkappa = True
                    if varkappaMatrix != []:
                        varkappaMatrixArray.append(varkappaMatrix)
                        varkappaMatrix = []
                    if not phiD == {}:
                        phiDL.append(phiD)
                        phiD = {}
                else:
                    if not nextNewVarkappa:
                        nextNewVarkappa = False
                    arrayLine = [int(item.split("|")[0]) for item in
                                    line.split("\t")]
                    varkappaMatrix.append(arrayLine)
                    for item in line.split("\t"):
                        spItem = item.split("|")
                        if len(spItem) > 1:
                            if int(spItem[0]) not in phiD.keys():
                                phiD[int(spItem[0])] = [float(spItem[1])]
                            else:
                                phiD[int(spItem[0])].append(
                                    float(spItem[1]))

            if varkappaMatrix != []:
                varkappaMatrixArray.append(varkappaMatrix)
            if not phiD == {}:
                phiDL.append(phiD)

        print "____>>> read_config: varkappaMatrixArray____"
        print varkappaMatrixArray
        print "_________end read_config:varkappaMatrixArray______________"
        print "____>>> read_config: phiDL____"
        print phiDL
        print "_________end read_config:phiDL______________"

        return varkappaMatrixArray, phiDL


class SRTree(object):

    """Summing rule while sampling tree structure"""

    def __init__(self, phiDL, timeOrderL):
        """initialize the self._srtree data structure"""
        self._srtree = []
        self._srtree_root = []
        for phiD, timeOrder in zip(phiDL, timeOrderL):
            tempTree = []
            for idx, t in timeOrder:
                # data index \t phi list \t time order \t path
                tempTree.append([idx, np.array(phiD[idx]), t, ""])
            self._srtree.append(sorted(tempTree, key=lambda item: item[2]))
            self._srtree_root.append(np.array([1] * len(phiD[idx])))


    def is_good_path(self, timeOrderIdx, dataIndex, newPath):
        """

        :timeOrderIdx: TODO
        :dataIndex: TODO
        :path: TODO
        :returns: TODO

        """

        index = -1
        idxL = [idx for idx, (j,phiL,t,path) in
                            enumerate(self._srtree[timeOrderIdx]) if j ==
                            dataIndex]
        if len(idxL) == 0:
            # 目标数据不在当前时序中
            return True
        else:
            # 目标数据在当前时序中，索引为index
            index = idxL[0]

        if index == 0:
            # 此处更新首节点的path
            # 此处添加一个根节点
            if np.prod(self._srtree_root[timeOrderIdx] -
                       self._srtree[timeOrderIdx][index][1]) > 0:
                return True
            else:
                return False
        else:
            cpl, pL, ancestorIndex = self.__load_most_recent_ancestor(timeOrderIdx, index, newPath)
            if cpl > 0:
                if np.prod(self._srtree[timeOrderIdx][ancestorIndex][1] -
                           self._srtree[timeOrderIdx][index][1]) > 0:
                    return True
                else:
                    return False
            else:
                if np.prod(self._srtree_root[timeOrderIdx] -
                        self._srtree[timeOrderIdx][index][1]) > 0:
                    return True
                else:
                    return False



    def update_path(self, timeOrderIdx, dataIndex, newPath):
        """Update path dictionary
        """
        # 获取所有比当前时序小的数据集合，在其中寻找最近祖先
        # 此处需要判断当前节点内是否有已经进行搜索过的节点

        index = -1
        idxL = [idx for idx, (j,phiL,t,path) in
                            enumerate(self._srtree[timeOrderIdx]) if j ==
                            dataIndex]
        if len(idxL) == 0:
            # 目标数据不在当前时序中
            return True
        else:
            # 目标数据在当前时序中，索引为index
            index = idxL[0]

        if index == 0:
            # 此处更新首节点的path
            # 此处添加一个根节点
            if np.prod(self._srtree_root[timeOrderIdx] -
                       self._srtree[timeOrderIdx][index][1]) > 0:
                self._srtree_root[timeOrderIdx] =\
                    self._srtree_root[timeOrderIdx] -\
                    self._srtree[timeOrderIdx][index][1]
                self._srtree[timeOrderIdx][index][3] = newPath
                return True
            else:
                return False
        else:
            cpl, pL, ancestorIndex = self.__load_most_recent_ancestor(timeOrderIdx, index, newPath)
            if cpl == len(newPath) and pL == self._srtree[timeOrderIdx][index][1]:
                # 且有相同的phi，实际上这种情况不可能出现，在无限位点假设条件下
                self._srtree[timeOrderIdx][index][3] = newPath
                return True
            elif cpl == len(newPath):
                # 且有相同的phi，实际上这种情况不可能出现，在无限位点假设条件下
                return False
            elif cpl > 0:
                if np.prod(self._srtree[timeOrderIdx][ancestorIndex][1] -
                           self._srtree[timeOrderIdx][index][1]) > 0:
                    self._srtree[timeOrderIdx][ancestorIndex][1] =\
                        self._srtree[timeOrderIdx][ancestorIndex][1] -\
                        self._srtree[timeOrderIdx][index][1]
                    self._srtree[timeOrderIdx][index][3] = newPath
                    return True
                else:
                    return False
            else:
                if np.prod(self._srtree_root[timeOrderIdx] -
                        self._srtree[timeOrderIdx][index][1]) > 0:
                    self._srtree_root[timeOrderIdx] =\
                        self._srtree_root[timeOrderIdx] -\
                        self._srtree[timeOrderIdx][index][1]
                    self._srtree[timeOrderIdx][index][3] = newPath
                    return True
                else:
                    return False

    def __load_most_recent_ancestor(self, timeOrderIdx, index, path):
        '''
        load most recent ancestor
        '''
        # Obtain maximum prefix
        #
        # EpsilonD should in the time order primarily and phi reverse order
        # secondly

        cplL = [(len(commonprefix([e, path])), pL, idx) if
                len(commonprefix([e, path])) == len(e) else (0, pL, idx)
                for idx, (j, pL, t, e) in
                enumerate(self._srtree[timeOrderIdx][0:index])]

        cpl, pL, idx = max(cplL, key=lambda item: item[0])

        if cpl > 0:
            return cpl, pL, idx
        else:
            return cpl, pL, -1

def main():
    # c2t = C2T("./crossing.text.example")
    # print c2t.toTimeOrder()
    cs2t = CS2T("./summingandcrossing.example.txt")
    timeOrder, negativeSD, phiDL = cs2t.transform()

    print timeOrder, negativeSD, phiDL

if __name__ == "__main__":
    main()

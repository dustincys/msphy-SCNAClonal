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


class C2T(object):

    """transform crossing rule to time ordered rule
    The configuration file in following format:

    # rule 1
    1   2   5   6   8
    5   8   6   1   2
    8   6   5   2   1
    # rule 2
    3   4   7   0
    3   7   4   0
    3   7   4   0

    Each digit is the index j of variant.
    Each line mean the phi_j reverse order of each sample $\theVector{\varkappa}
    _{k,l}$.
    """

    def __init__(self, configFilePath):
        """initialize parameter

        :configFilePath: The phi order configuration

        """
        self._configFilePath = configFilePath
        # self.varkappaMatrixArray = {}
        # self.negativeSD = {}

    def toTimeOrder(self, varkappaMatrixArray):
        """
        :returns: TODO

        """
        for varkappaMatrix in varkappaMatrixArray:
            for ki in range(len(varkappaMatrix)-1):
                tArray = [len(
                    set(varkappaMatrix[-1][kj:]) - set(
                        varkappaMatrix[ki][np.where(
                            np.array(varkappaMatrix[ki]) ==\
                            varkappaMatrix[-1][kj])[0][0]:]))
                    for kj in range(len(varkappaMatrix[-1][:]))]

                for ti in range(len(varkappaMatrix[-1])):
                    if varkappaMatrix[-1][ti] not in\
                            negativeSD.keys():
                        negativeSD[varkappaMatrix[-1][ti]] = set()
                    for tj in range(ti+1, len(varkappaMatrix[-1])):
                        if tArray[tj] < tArray[ti]:
                            if varkappaMatrix[-1][ti] not in\
                                    negativeSD.keys():
                                negativeSD[varkappaMatrix[-1][ti]] =\
                                    set([varkappaMatrix[-1][tj]])
                            else:
                                negativeSD[varkappaMatrix[-1][ti]].add(
                                    varkappaMatrix[-1][tj])
        idxOrder = sorted(negativeSD.keys())
        timeOrder = [(idx, len(negativeSD[idx])) for idx in idxOrder]

        return timeOrder, negativeSD

    def read_config_file(self):
        """TODO: Docstring for __read_config_file.

        :uu: TODO
        :returns: TODO

        """
        with open(self._configFilePath) as configFile:
            nextNewVarkappa = False

            varkappaMatrix = []
            varkappaMatrixArray = []

            for line in configFile:
                line = line.strip()
                if line == "":
                    continue

                if line.startswith("#"):
                    nextNewVarkappa = True

                    if varkappaMatrix != []:
                        varkappaMatrixArray.append(varkappaMatrix)
                        varkappaMatrix = []

                else:
                    if not nextNewVarkappa:
                        nextNewVarkappa = False

                    arrayLine = [int(item) for item in line.split("\t")]
                    varkappaMatrix.append(arrayLine)
                pass
            pass

            if varkappaMatrix != []:
                varkappaMatrixArray.append(varkappaMatrix)

        return varkappaMatrixArray

class CS2T(C2T):

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
        C2T.__init__(self, configFilePath)

        self.varkappaMatrixArray, self.currentVarkappa, self.phiDL,\
            self.currentPhiD = self.read_config_file()

        # print self.varkappaMatrixArray
        # print self.currentVarkappa
        # print self.phiDL
        # print self.currentPhiD

        self.timeOrder, self.negativeSD = self.toTimeOrder(
            self.varkappaMatrixArray, self.currentVarkappa)

        self.epsilonDL = []
        # print self.phiDL
        for phiDIndex in range(len(self.phiDL)):
            epsilonD = {}
            for key in self.phiDL[phiDIndex].keys():
                epsilonD[key] = ""
            self.epsilonDL.append(epsilonD)

    def toTimeOrder(self, varkappaMatrixArray, currentVarkappa):
        """ tranform to time order
        :returns: time order

        """
        negativeSD = {}
        totalTimeOrderD = {}

        for varkappaMatrix in varkappaMatrixArray:
            for ki in range(len(varkappaMatrix)):
                currentTempVarkappa = np.array([currentVarkappa[item] for item
                                                in range(len(currentVarkappa))
                                                if currentVarkappa[item] in
                                                varkappaMatrix[ki]])
                negativeSetList = [set(currentTempVarkappa[kj:]) -
                                   set(varkappaMatrix[ki][ np.where(
                                       np.array(varkappaMatrix[ki]) ==
                                       currentTempVarkappa[kj])[0][0]:]) for kj
                                   in range(len(currentTempVarkappa[:]))]

                for kj in range(len(currentTempVarkappa[:])):
                    if currentTempVarkappa[kj] not in negativeSD.keys():
                        negativeSD[currentTempVarkappa[kj]] =\
                            negativeSetList[kj]
                    else:
                        negativeSD[currentTempVarkappa[kj]] = set.union(
                            negativeSetList[kj],
                            negativeSD[currentTempVarkappa[kj]])

            tempTimeOrder = self.__idxOrder_to_timeOrder(negativeSD,
                                                          varkappaMatrix[0],
                                                          currentVarkappa)
            for tto in tempTimeOrder:
                totalTimeOrderD[tto[0]] = tto[1]

        timeOrder = [(idx, totalTimeOrderD[idx]) if idx in
                     totalTimeOrderD.keys() else (idx, 0) for idx in
                     currentVarkappa]

        return timeOrder, negativeSD


    def __idxOrder_to_timeOrder(self, negativeSD, targetData,
                                currentVarkappa):
        tQueue = deque()
        tStack = deque()
        currentTempVarkappa = np.array([currentVarkappa[item] for item in
                                        range(len(currentVarkappa)) if
                                        currentVarkappa[item] in targetData])
        nsL = [(idx, deepcopy(negativeSD[idx])) if idx in
                negativeSD.keys() else (idx, set()) for idx in
                currentTempVarkappa]
        stackNegativeSet = set()

        for idx, ns in nsL:
            if len(stackNegativeSet) == 0:
                if len(tStack) > 0:
                    tQueue.append(tStack)
                tStack = deque()
            tStack.append([idx, ns])
            stackNegativeSet = set.union(ns, stackNegativeSet)
            if idx in stackNegativeSet:
                stackNegativeSet.remove(idx)
        if len(tStack) > 0:
            tQueue.append(tStack)

        timeOrder = []
        currentTime = 0

        while len(tQueue) != 0:
            tStack = tQueue.popleft()
            tStack = deque(sorted(tStack, key=lambda item: len(item[1]),
                                  reverse=True))
            lastTimeIdxS = set()
            while len(tStack) != 0:
                if len(tStack[-1][1]) != 0:
                    for item in tStack:
                        item[1] = item[1] - lastTimeIdxS
                    tStack = deque(sorted(tStack, key=lambda item: len(item[1]),
                                          reverse=True))
                    lastTimeIdxS.clear()
                    currentTime = currentTime + 1

                idx, nsl = tStack.pop()
                timeOrder.append((idx, currentTime))
                lastTimeIdxS.add(idx)

        return timeOrder

    def update_epsilon(self, dataIndex, epsilon):
        """Update epsilon dictionary
        """
        indexes = [i for i, epsilonD in enumerate(self.epsilonDL) if dataIndex
                   in epsilonD.keys()]
        if len(indexes) == 0:
            return True
        else:
            index = indexes[0]
            ancestorIndex = self.__load_most_recent_ancestor(
                self.epsilonDL[index], dataIndex, epsilon)
            if ancestorIndex != -1:
                if np.prod(self.phiDL[index][ancestorIndex] -
                           self.phiDL[index][dataIndex]) > 0:
                    self.phiDL[index][ancestorIndex] =\
                        self.phiDL[index][ancestorIndex] -\
                        self.phiDL[index][dataIndex]
                    return True
                else:
                    return False
            else:
                return True

    def __load_most_recent_ancestor(self, epsilonD, dataIndex, epsilon):
        ''' load most recent ancestor
        '''
        cplL = [(key, len(commonprefix([epsilonD[key], epsilon]))) if key !=
                dataIndex else (dataIndex, 0) for key in epsilonD ]
        _, (ancestorIndex, cpl) = max(enumerate(cplL), key=lambda v: v[1][1])
        if cpl > 0:
            return ancestorIndex
        else:
            return -1

    def read_config_file(self):
        """Read config file: index|phi, phi descreasing order
        # Current
        3|0.7   8|0.6   7|0.5   4|0.3   6|0.3   5|0.2   0|0.1   2|0.1   1|0.05
        # rule 1
        1|0.8   2|0.7   5|0.4   6|0.3   8|0.2
        5   8   6   1   2
        # rule 2
        3   4   7   0
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
                    if line.find("current") >= 0:
                        nextCurrent = True
                    else:
                        nextNewVarkappa = True
                        if varkappaMatrix != []:
                            varkappaMatrixArray.append(varkappaMatrix)
                            varkappaMatrix = []
                        if not phiD == {}:
                            phiDL.append(phiD)
                            phiD = {}
                else:
                    if nextCurrent:
                        currentVarkappa = [int(item.split("|")[0]) for item in
                                           line.split("\t")]
                        # print "currentVarkappa"
                        # print currentVarkappa
                        for item in line.split("\t"):
                            spItem = item.split("|")
                            if len(spItem) > 1:
                                if int(spItem[0]) not in phiD.keys():
                                    currentPhiD[int(spItem[0])] =\
                                        [float(spItem[1])]
                                else:
                                    currentPhiD[int(spItem[0])].append(
                                        float(spItem[1]))
                        nextCurrent = False
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

            if nextCurrent:
                currentVarkappa = [int(item.split("|")[0]) for item in
                                    line.split("\t")]
                for item in line.split("\t"):
                    spItem = item.split("|")
                    if len(spItem) > 1:
                        if int(spItem[0]) not in phiD.keys():
                            currentPhiD[int(spItem[0])] =\
                                [float(spItem[1])]
                        else:
                            currentPhiD[int(spItem[0])].append(
                                float(spItem[1]))

        # print "____>>> read_config: varkappaMatrixArray____"
        # print varkappaMatrixArray
        # print "_________end read_config:varkappaMatrixArray______________"
        # print "____>>> read_config: currentVarkappa____"
        # print currentVarkappa
        # print "_________end read_config:currentVarkappa______________"
        # print "____>>> read_config: phiDL____"
        # print phiDL
        # print "_________end read_config:phiDL______________"

        return varkappaMatrixArray, currentVarkappa, phiDL, currentPhiD

def main():
    # c2t = C2T("./crossing.text.example")
    # print c2t.toTimeOrder()
    cs2t = CS2T("./summingandcrossing.example.txt")
    print cs2t.timeOrder

if __name__ == "__main__":
    main()

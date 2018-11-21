#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: crossing.py
#          Desc: use crossing rule to find topological restriction
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-11-15 16:37:09
#       History:
# =============================================================================
'''
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
        self._varkappaMatrixArray = self.read_config_file(self._configFilePath)
        self._negativeSetDict = {}

    def toTimeOrder(self):
        """
        :returns: TODO

        """
        for varkappaMatrix in self._varkappaMatrixArray:
            for ki in range(len(varkappaMatrix)-1):
                tArray = [len(
                    set(varkappaMatrix[-1][kj:]) - set(
                        varkappaMatrix[ki][np.where(
                            np.array(varkappaMatrix[ki]) == varkappaMatrix[-1][kj])[0][0]:]))
                    for kj in range(len(varkappaMatrix[-1][:]))]

                for ti in range(len(varkappaMatrix[-1])):
                    if varkappaMatrix[-1][ti] not in\
                            self._negativeSetDict.keys():
                        self._negativeSetDict[varkappaMatrix[-1][ti]] = set()
                    for tj in range(ti+1, len(varkappaMatrix[-1])):
                        if tArray[tj] < tArray[ti]:
                            if varkappaMatrix[-1][ti] not in\
                                    self._negativeSetDict.keys():
                                self._negativeSetDict[varkappaMatrix[-1][ti]] =\
                                    set([varkappaMatrix[-1][tj]])
                            else:
                                self._negativeSetDict[varkappaMatrix[-1][ti]].add(
                                    varkappaMatrix[-1][tj])
        idxOrder = sorted(self._negativeSetDict.keys())
        timeOrder = [(idx, len(self._negativeSetDict[idx])) for idx in idxOrder]

        return timeOrder, self._negativeSetDict


    def read_config_file(self, configFilePath):
        """TODO: Docstring for __read_config_file.

        :uu: TODO
        :returns: TODO

        """
        with open(configFilePath) as configFile:
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
        pass


def main():
    c2t = C2T("./crossing.text")
    print c2t.toTimeOrder()

if __name__ == "__main__":
    main()

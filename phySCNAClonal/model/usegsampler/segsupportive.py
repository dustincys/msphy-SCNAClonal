#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: segsupportive.py
#          Desc: uniform sampler for multiple supportive ranges
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-06-08 09:48:22
#       History:
# =============================================================================
'''

from gwpy.segments import Segment, SegmentList
import numpy as np
import sys


class MultiRangeSampler(object):

    """Sample from multiple supportive ranges"""

    def __init__(self, minU=0, maxU=1):
        """initialized by default"""
        self._minU = minU
        self._maxU = maxU
        self._supportiveRanges = self._get_supportive_ranges()
        self._cumLens = self._get_cumulative_lens()

    def _get_cumulative_lens(self):
        lens = np.array([seg[1] - seg[0] for seg in self._supportiveRanges])
        # here need to transform to float
        lens = lens * 1.0 / sum(lens)
        return np.cumsum(lens)

    def _get_supportive_ranges(self):
        return SegmentList([Segment(self._minU, self._maxU)])

    @property
    def lowerBoundary(self):
        return self._minU

    @property
    def upperBoundary(self):
        return self._maxU


    def assign_supportive(self, ranges):
        """

        :ranges: TODO
        :returns: TODO

        """
        self._supportiveRanges = ranges
        self._supportiveRanges.coalesce()

        if 0 < len(self._supportiveRanges):
            self._minU = self._supportiveRanges[0][0]
            self._maxU = self._supportiveRanges[-1][1]
            self._cumLens = self._get_cumulative_lens()
        else:
            # this is not possible
            self._minU = 0
            self._maxU = 0
            self._cumLens = np.array([])
            raise Exception("remove range error!")


    def remove(self, ranges):
        self._supportiveRanges = self._supportiveRanges - ranges
        self._supportiveRanges.coalesce()

        if 0 < len(self._supportiveRanges):
            self._minU = self._supportiveRanges[0][0]
            self._maxU = self._supportiveRanges[-1][1]
            self._cumLens = self._get_cumulative_lens()
            return True
        else:
            # this is not possible
            self._minU = 0
            self._maxU = 0
            self._cumLens = np.array([])
            return False


    def removeLeft(self, boundary):
        slRm = SegmentList([Segment(self._minU, boundary)])
        self._supportiveRanges = self._supportiveRanges - slRm
        self._supportiveRanges.coalesce()

        if 0 < len(self._supportiveRanges):
            self._minU = self._supportiveRanges[0][0]
            self._cumLens = self._get_cumulative_lens()
        else:
            self._minU = self._maxU
            self._cumLens = np.array([])

    def removeRight(self, boundary):
        slRm = SegmentList([Segment(boundary, self._maxU)])
        self._supportiveRanges = self._supportiveRanges - slRm
        self._supportiveRanges.coalesce()

        if 0 < len(self._supportiveRanges):
            self._maxU = self._supportiveRanges[-1][1]
            self._cumLens = self._get_cumulative_lens()
        else:
            self._maxU = self._minU
            self._cumLens = np.array([])

    def sample(self):
        if 0 == len(self._supportiveRanges):
            return self._maxU
        else:
            index = self._getIndex()
            # Samples are uniformly distributed over the half-open interval
            # [low, high) (includes low, but excludes high). I

            leftBorder = self._supportiveRanges[index][0] +\
                sys.float_info.min * (self._supportiveRanges[index][1] -
                                      self._supportiveRanges[index][0])

            return np.random.uniform(leftBorder,
                                     self._supportiveRanges[index][1])

    def _getIndex(self):
        prn = np.random.uniform(0, 1)
        return np.min(np.where(self._cumLens > prn)[0])


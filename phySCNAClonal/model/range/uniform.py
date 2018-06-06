#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: uniform.py
#          Desc: uniform sampler for multiple supportive ranges
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-06-06 12:39:43
#       History:
# =============================================================================
'''

from gwpy.segments import Segment, SegmentList
import numpy as np


class MultiRangeSampler(object):

    """Sample from multiple supportive ranges"""

    def __init__(self, minU=0, maxU=1):
        """initialized by default"""
        self._minU = minU
        self._maxU = maxU
        self._supportiveRanges = self._get_supportive_ranges()
        self._cumLens = self._get_cumulative_lens()

    def _get_cumulative_lens(self):
        lens = np.array([seg[1] - seg[0] for seg in selg._supportiveRanges])
        lens = lens / sum(lens)
        return np.cumsum(lens)

    def _get_supportive_ranges(self):
        return SegmentList([Segment(self._minU, self._maxU)])

    def remove(self, ranges):
        self._supportiveRanges = self._supportiveRanges - ranges
        self._supportiveRanges.coalesce()
        self._cumLens = self._get_cumulative_lens()

    def sample(self, number):
        index = self._getIndex()
        return np.random.uniform(self._supportiveRanges[index][0],
                                 self._supportiveRanges[index][1])

    def _getIndex(self):
        prn = np.random.uniform(0, 1)
        return np.min(np.where(self._cumLens > prn)[0])


#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: mcmc.py
#          Desc: The PYMC model for gc correction
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-11-17 12:00:33
#       History:
# =============================================================================
'''

import pymc
import numpy as np
from scipy.stats import gaussian_kde

import constants
from scipy.signal import argrelextrema

import heapq


class MCMCLM(object):

    """The MCMC model for linear regression, return the slope and inlier"""

    def __init__(self, data, n, tau, max_copynumber):
        """Initialize the MCMCLM model

        :data: the segment data object
        :n: the sampling number
        :tau: the subclone number
        :max_copynumber: maximum copy number
        """
        self._data = data
        self._n = n
        self._tau = tau
        self._max_copynumber = max_copynumber

        # parameters
        self._downGCBoundaryPercentile = constants.DOWN_GC_BOUNDARY_PERCENTILE
        self._upGCBoundaryPercentile = constants.UP_GC_BOUNDARY_PERCENTILE

        self._downLOGABoundaryPercentile = \
            constants.DOWN_LOGA_BOUNDARY_PERCENTILE
        self._upLOGABoundaryPercentile = \
            constants.UP_LOGA_BOUNDARY_PERCENTILE

        self._slope_range = constants.SLOPE_RANGE

        self._zoom_p = constants.ZOOM_P
        self._x_zoom_in_factor = constants.X_ZOOM_IN_FACTOR

        self._y, self._x = self._getSampledData()
        self._x_zoommed = self._zoomx()
        if self._x_zoommed:
            self._x = self._x * self._x_zoom_in_factor

    def _zoomx(self):
        spanp = (max(self._y) - min(self._y)) / (max(self._x) - min(self._x))
        if spanp > self._zoom_p:
            return True
        else:
            return False

    def run(self):
        """Correct Y
        return: the corrected Y
        """
        slope_best, intercept_best = self._getMCPosterior(self._y, self._x)
        if self._x_zoommed:
            slope_best = slope_best * self._x_zoom_in_factor

        return slope_best, intercept_best

    def _getMCPrior(self, y, x):
        x_down_ceil = np.percentile(x, self._downGCBoundaryPercentile)
        x_median = np.percentile(x, 50)
        x_up_floor = np.percentile(x, self._upGCBoundaryPercentile)

        y_up = y[np.logical_and(x < x_up_floor, x > x_median)]
        y_down = y[np.logical_and(x > x_down_ceil, x < x_median)]
        x_up = x[np.logical_and(x < x_up_floor, x > x_median)]
        x_down = x[np.logical_and(x > x_down_ceil, x < x_median)]

        y_up_y = np.percentile(y_up, 50)
        x_up_x = np.percentile(x_up, 50)

        y_down_y = np.percentile(y_down, 50)
        x_down_x = np.percentile(x_down, 50)

        m = (y_up_y - y_down_y) * 1.0 / (x_up_x - x_down_x)
        c = y_down_y - m * x_down_x

        return m, c

    def _correctY(self, y, x, slope, intercept):
        K = np.percentile(y, 50)
        A = slope * x + intercept
        return y - A + K

    def _getMCPosterior(self, y_with_outlier, x):
        m, c = self._getMCPrior(y_with_outlier, x)
        slope = pymc.Uniform('slope', m-self._slope_range, m+self._slope_range)

        def log_posterior_likelihood_of_slope(
                y_with_outlier,
                slope
                ):
            y_corrected = self._correctY(y_with_outlier, x, slope, 0)
            y_density = gaussian_kde(y_corrected)

            y_down = min(y_corrected)
            y_up = max(y_corrected)
            y_xs = np.linspace(y_down, y_up,
                               1000*self._tau*self._max_copynumber)
            y_ys = y_density(y_xs)
            peaks = argrelextrema(y_ys, np.greater)

            prob = sum(heapq.nlargest(self._max_copynumber,  y_ys[peaks[0]]))

            return prob

        slope_distribution = pymc.stochastic_from_dist(
            'slope_distribution',
            logp=log_posterior_likelihood_of_slope,
            dtype=np.float,
            mv=True)

        slope_dist = slope_distribution('slope_dist',
                                        slope=slope,
                                        observed=True,
                                        value=y_with_outlier)

        model = dict(slope_dist=slope_dist,
                     slope=slope
                     )

        M = pymc.MAP(model)
        M.fit()

        slope_best = M.slope.value
        y_median = np.percentile(y_with_outlier, 50)
        x_median = x[sum(y_with_outlier < y_median)]
        intercept_best = y_median - slope_best * x_median

        return slope_best, intercept_best

    def _getSampledData(self):
        if 0 == self._n:
            sampledSegs = self._data.segments
        else:
            sampledSegs = np.random.choice(self._data.segments, self._n)

        print "all sample: {}".format(len(sampledSegs))

        x0 = np.array(map(lambda seg: seg.gc, sampledSegs))
        y0 = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                          np.log(seg.normal_reads_num + 1), sampledSegs))
        l = sorted(zip(y0, x0), reverse=True)
        y0, x0 = [list(t) for t in zip(*l)]

        return np.array(y0), np.array(x0)

    def getPeakRange(self, slope_best):
        """TODO: Docstring for getPeakRange.

        :y: TODO
        :x: TODO
        :slope: TODO
        :returns: TODO

        """
        if self._x_zoommed:
            slope = slope_best / self._x_zoom_in_factor
        y_corrected = self._correctY(self._y, self._x, slope, 0)
        y_density = gaussian_kde(y_corrected)

        y_down = min(y_corrected)
        y_up = max(y_corrected)
        y_xs = np.linspace(y_down, y_up, 1000*self._tau*self._max_copynumber)
        y_ys = y_density(y_xs)
        peaks = argrelextrema(y_ys, np.greater)
        y_ys_nl = np.array(heapq.nlargest(self._max_copynumber,
                                          y_ys[peaks[0]]))
        idx = np.nonzero(y_ys_nl[:,None] == y_ys)[1]
        y_xs_nl = np.sort(y_xs[idx])
        pr = 0
        for i in range(len(y_xs_nl) - 1):
            pr = pr + y_xs_nl[i+1] - y_xs_nl[i]
        pr = pr * 1.0 / (len(y_xs_nl) - 1)

        print "idx = "
        print idx
        print "y_xs_nl = "
        print y_xs_nl

        return pr

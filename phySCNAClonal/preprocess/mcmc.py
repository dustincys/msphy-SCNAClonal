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
#    LastChange: 2018-03-03 08:55:12
#       History:
# =============================================================================
'''

import heapq

import numpy as np
import pymc
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde

import phySCNAClonal.constants as constants


class MCMCLM(object):

    """The MCMC model for linear regression, return the slope and inlier"""

    def __init__(self, data, n, tau, maxCopyNumber):
        """Initialize the MCMCLM model

        :data: the segment data object
        :n: the sampling number
        :tau: the subclone number
        :maxCopyNumber: maximum copy number
        """
        self._data = data
        self._n = n
        self._tau = tau
        self._maxCopyNumber = maxCopyNumber

        # parameters
        self._downGCBoundaryPercentile = constants.DOWN_GC_BOUNDARY_PERCENTILE
        self._upGCBoundaryPercentile = constants.UP_GC_BOUNDARY_PERCENTILE

        self._downLOGABoundaryPercentile = \
            constants.DOWN_LOGA_BOUNDARY_PERCENTILE
        self._upLOGABoundaryPercentile = \
            constants.UP_LOGA_BOUNDARY_PERCENTILE

        self._slopeRange = constants.SLOPE_RANGE

        self._zoomP = constants.ZOOM_P
        self._xZoomInFactor = constants.X_ZOOM_IN_FACTOR

        self._y, self._x = self._getSampledData()
        self._xZoommed = self._zoomx()
        if self._xZoommed:
            self._x = self._x * self._xZoomInFactor

    def _zoomx(self):
        spanp = (max(self._y) - min(self._y)) / (max(self._x) - min(self._x))
        if spanp > self._zoomP:
            return True
        else:
            return False

    def run(self):
        """Correct Y
        return: the corrected Y
        """
        slopeBest, interceptBest = self._getMCPosterior(self._y, self._x)
        if self._xZoommed:
            slopeBest = slopeBest * self._xZoomInFactor

        return slopeBest, interceptBest

    def _getMCPrior(self, y, x):
        xDownCeil = np.percentile(x, self._downGCBoundaryPercentile)
        xMedian = np.percentile(x, 50)
        xUpFloor = np.percentile(x, self._upGCBoundaryPercentile)

        yUp = y[np.logical_and(x < xUpFloor, x > xMedian)]
        yDown = y[np.logical_and(x > xDownCeil, x < xMedian)]
        xUp = x[np.logical_and(x < xUpFloor, x > xMedian)]
        xDown = x[np.logical_and(x > xDownCeil, x < xMedian)]

        yUpY = np.percentile(yUp, 50)
        xUpX = np.percentile(xUp, 50)

        yDownY = np.percentile(yDown, 50)
        xDownX = np.percentile(xDown, 50)

        m = (yUpY - yDownY) * 1.0 / (xUpX - xDownX)
        c = yDownY - m * xDownX

        return m, c

    def _correctY(self, y, x, slope, intercept):
        K = np.percentile(y, 50)
        A = slope * x + intercept
        return y - A + K

    def _getMCPosterior(self, yWithOutlier, x):
        m, c = self._getMCPrior(yWithOutlier, x)
        slope = pymc.Uniform('slope', m-self._slopeRange, m+self._slopeRange)

        def log_posterior_likelihood_of_slope(
                yWithOutlier,
                slope
                ):
            yCorrected = self._correctY(yWithOutlier, x, slope, 0)
            yDensity = gaussian_kde(yCorrected)

            yDown = min(yCorrected)
            yUp = max(yCorrected)
            yXs = np.linspace(yDown, yUp,
                               1000*self._tau*self._maxCopyNumber)
            yYs = yDensity(yXs)
            peaks = argrelextrema(yYs, np.greater)

            prob = sum(heapq.nlargest(self._maxCopyNumber,  yYs[peaks[0]]))

            return prob

        slope_distribution = pymc.stochastic_from_dist(
            'slope_distribution',
            logp=log_posterior_likelihood_of_slope,
            dtype=np.float,
            mv=True)

        slope_dist = slope_distribution('slope_dist',
                                        slope=slope,
                                        observed=True,
                                        value=yWithOutlier)

        model = dict(slope_dist=slope_dist, slope=slope)

        M = pymc.MAP(model)
        M.fit()

        slopeBest = M.slope.value
        y_median = np.percentile(yWithOutlier, 50)
        xMedian = x[sum(yWithOutlier < y_median)]
        interceptBest = y_median - slopeBest * xMedian

        return slopeBest, interceptBest

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

    def getPeakRange(self, slopeBest):
        """TODO: Docstring for getPeakRange.

        :y: TODO
        :x: TODO
        :slope: TODO
        :returns: TODO

        """
        if self._xZoommed:
            slope = slopeBest / self._xZoomInFactor
        yCorrected = self._correctY(self._y, self._x, slope, 0)
        yDensity = gaussian_kde(yCorrected)

        yDown = min(yCorrected)
        yUp = max(yCorrected)
        yXs = np.linspace(yDown, yUp, 1000*self._tau*self._maxCopyNumber)
        yYs = yDensity(yXs)
        peaks = argrelextrema(yYs, np.greater)
        yYsNl = np.array(heapq.nlargest(self._maxCopyNumber,
                                          yYs[peaks[0]]))
        idx = np.nonzero(yYsNl[:,None] == yYs)[1]
        yXsNl = np.sort(yXs[idx])
        pr = 0
        for i in range(len(yXsNl) - 1):
            pr = pr + yXsNl[i+1] - yXsNl[i]
        pr = pr * 1.0 / (len(yXsNl) - 1)

        print "idx = "
        print idx
        print "yXsNl = "
        print yXsNl

        return pr

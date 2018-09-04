#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: plotGC.py
#          Desc: plot gc stripes and interactively adjust the lm
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-11 15:55:22
#       History:
# =============================================================================
'''

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import colorConverter
from matplotlib.path import Path
from matplotlib.widgets import Button, LassoSelector, Slider

import seaborn as sns
from itertools import cycle


class GCStripePlot():

    """plot the gc stripe and interactively adjust the lm line"""

    def __init__(self, segments, n):
        """init
        segments:  all the segments of data, need to be sampled
        n:  the sample number

        return none
        """

        self.segments = segments
        self.n = n

#       parameters for plot
        self.m0 = 0
        self.c0 = 0
        self.alpha0 = 0.02
        self.area0 = 10

        self.m = 0
        self.c = 0

        self.m1 = 0
        self.c1 = 0

        self.alpha = 1
        self.area = 10

        self.x = None
        self.y = None

        self.colorin = colorConverter.to_rgba('red', self.alpha0)
        self.colorout = colorConverter.to_rgba('blue', self.alpha0)

    def output(self):
        """
        :returns: TODO

        """
        return self.x, self.y, self.m, self.c

    def _updateM(self, mt):
        """

        :m: TODO
        :returns: TODO

        """

        self.m = np.tan(np.arctan(self.m1) + mt)

    def _updateC(self, ct):
        """

        :ct: TODO
        :returns: TODO

        """
        self.c = self.c1 + ct

    def sampleln(self, n_list, m):
        """
        :n: sampling number
        :m: iteration number
        :returns: TODO

        """
        with open("/home/dustin/temp/sampling_result.txt", "w") as outFile:
            outFile.write("n\tm\ta\tb\n")
            for n in n_list:
                for i in range(m):
                    sampledSegs = np.random.choice(self.segments, n)
                    x0 = np.array(map(lambda seg: seg.gc, sampledSegs))
                    y0 = np.array(
                        map(lambda seg: np.log(seg.tReadNum + 1) -
                            np.log(seg.nReadNum + 1), sampledSegs))
                    A = np.vstack([x0, np.ones(len(x0))]).T
                    a, b = np.linalg.lstsq(A, y0)[0]
                    outFile.write("{0}\t{1}\t{2}\t{3}\n".format(n, m, a, b))

        print "finished sampling"

    def plot(self):

        sampledSegs = np.random.choice(self.segments, self.n)
        x0 = np.array(map(lambda seg: seg.gc, sampledSegs))
        y0 = np.array(map(lambda seg: np.log(seg.tReadNum + 1) -
                          np.log(seg.nReadNum + 1), sampledSegs))

        self.x = x0
        self.y = y0

        fig, ax = plt.subplots()

        pts = ax.scatter(x0, y0, s=self.area0, alpha=self.alpha0, color="b")
        ax.set_xlabel("GC content")
        ax.set_ylabel(r'$\log(D^T/D^N)$')
        plt.subplots_adjust(bottom=0.45)

        A = np.vstack([x0, np.ones(len(x0))]).T
        self.m0, self.c0 = np.linalg.lstsq(A, y0)[0]
        self.m1, self.c1 = self.m0, self.c0
        self.m, self.c = self.m0, self.c0

        # xseq = np.arange(min(x0), max(x0), (max(x0) - min(x0))/100.0)
        # fl, = ax.plot(xseq, self.m0*xseq + self.c0, 'r', label='Fitted line')
        # hl = ax.axhline(y=np.median(y0), linewidth=1, color='black')
        # vl = ax.axvline(x=np.median(x0), linewidth=1, color='black')

        axcolor = 'lightgoldenrodyellow'
        axm = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
        axc = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
        axalpha = plt.axes([0.25, 0.2, 0.65, 0.03], axisbg=axcolor)
        axarea = plt.axes([0.25, 0.25, 0.65, 0.03], axisbg=axcolor)

        sm = Slider(axm, 'slope_delta', -np.pi/4, np.pi/4, valinit=0)
        sc = Slider(axc, 'interception_delta', -(max(self.y) - min(self.y)) / 4,
                    (max(self.y) - min(self.y)) / 4, valinit=0)
        salpha = Slider(axalpha, 'alpha', 0, 0.8, valinit=self.alpha0)
        sarea = Slider(axarea, 'area', 1, 50, valinit=self.area0)

        def update_m(val):
            """
            :returns: TODO

            """
            x_axis = np.median(self.x)
            y_axis = x_axis * self.m + self.c
            self._updateM(sm.val)
            self.c = y_axis - self.m * x_axis
            # fl.set_ydata(self.m*xseq + self.c)
            fig.canvas.draw_idle()

        def update_c(val):
            """
            :returns: TODO

            """
            self._updateC(sc.val)
            # fl.set_ydata(self.m*xseq + self.c)
            fig.canvas.draw_idle()

        def update_alpha(val):
            """
            :returns: TODO

            """
            self.alpha = salpha.val
            pts.set_alpha(self.alpha)
            # fig.canvas.draw_idle()

        def update_area(val):
            """
            :returns: TODO

            """
            self.area = sarea.val
            pts.set_sizes(np.ones(len(self.x)) * self.area)
            fig.canvas.draw_idle()

        sm.on_changed(update_m)
        sc.on_changed(update_c)
        salpha.on_changed(update_alpha)
        sarea.on_changed(update_area)

        select_col = lambda array, colNum: map(lambda x: x[colNum], array)

        def onselect(verts):
            path = Path(verts)
            xys = pts.get_offsets()
            ind = np.nonzero([path.contains_point(xy) for xy in xys])[0]
            self.x = np.array(select_col(xys[ind], 0))
            self.y = np.array(select_col(xys[ind], 1))

            self.colorin = colorConverter.to_rgba('red', self.alpha0)
            self.colorout = colorConverter.to_rgba('blue', self.alpha0)
            facecolors = np.array([self.colorout for xy in xys])
            facecolors[ind] = self.colorin
            pts.set_facecolors(facecolors)

            A = np.vstack([self.x, np.ones(len(self.x))]).T
            self.m, self.c = np.linalg.lstsq(A, self.y)[0]
            self.m1, self.c1 = self.m, self.c
            # fl.set_ydata(self.m*xseq + self.c)
            # hl.set_ydata(np.median(self.y))

            sm.reset()
            sc.reset()

            fig.canvas.draw_idle()

        lasso = LassoSelector(ax, onselect)

        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

        def reset(event):
            sm.reset()
            sc.reset()
            self.colorout = colorConverter.to_rgba('blue', self.alpha0)
            facecolors = np.array([self.colorout for y in y0])
            pts.set_facecolors(facecolors)
            self.y = y0
            self.x = x0
            A = np.vstack([self.x, np.ones(len(self.x))]).T
            self.m, self.c = np.linalg.lstsq(A, self.y)[0]
            self.m1, self.c1 = self.m, self.c
            # fl.set_ydata(self.m*xseq + self.c)
            # hl.set_ydata(np.median(self.y))

        button.on_clicked(reset)

        resetax = plt.axes([0.5, 0.025, 0.1, 0.04])
        button_exit = Button(resetax, 'Ok', color=axcolor, hovercolor='0.975')

        def ok_exit(event):
            plt.close()

        button_exit.on_clicked(ok_exit)

        plt.show()


class GCStripePoolPlot(GCStripePlot):

    """plot the stripe pool, add different color to each stripe"""

    def __init__(self, stripePool):
        """ Initialize stripePool plot """
        GCStripePlot.__init__(self, stripePool.segPool.segments,
                              len(stripePool.segPool.segments))

        self._stripePool = stripePool
        self.alpha0 = 1

    def getColorVec(self):

        colorVec = np.ones((1, len(self._stripePool.segPool.segments)))[0] * (-1)
        print "____>>> getColorVec: colorVec____"
        print colorVec
        print "_________end getColorVec:colorVec______________"

        for idx, sp in enumerate(self._stripePool.stripes):
            if sp.tag == "-999999":
                continue
            tempIdxl = np.array(sp.segsIdxL)
            colorVec[tempIdxl] = idx

        return colorVec

    def plot(self):
        # sampledSegs = np.random.choice(self.segments, self.n)
        sampledSegs = self.segments
        x0 = np.array(map(lambda seg: seg.gc, sampledSegs))
        y0 = np.array(map(lambda seg: np.log(seg.tReadNum + 1) -
                          np.log(seg.nReadNum + 1), sampledSegs))

        self.x = x0
        self.y = y0

        fig, ax = plt.subplots()

        # here it needs to be filtered out segments not in the stripes
        colorVec = self.getColorVec()

        inIndices = np.where(colorVec != -1)[0]
        pts = ax.scatter(x0[inIndices], y0[inIndices], s=self.area0,
                         alpha=self.alpha0, c=colorVec[inIndices])

        spBaselines = filter(lambda sp:sp.tag == "-999999",
                            self._stripePool.stripes)
        if spBaselines != []:
            baselineIndices = np.array(spBaselines[0].segsIdxL)
            pts2 = ax.scatter(x0[baselineIndices], y0[baselineIndices], s=100,
                            alpha=1, marker='s', c="0000")

        ax.set_xlabel("GC content")
        ax.set_ylabel(r'$\log(D^T/D^N)$')
        plt.subplots_adjust(bottom=0.45)

        A = np.vstack([x0, np.ones(len(x0))]).T
        self.m0, self.c0 = np.linalg.lstsq(A, y0)[0]
        self.m1, self.c1 = self.m0, self.c0
        self.m, self.c = self.m0, self.c0

        axcolor = 'lightgoldenrodyellow'
        axm = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
        axc = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
        axalpha = plt.axes([0.25, 0.2, 0.65, 0.03], axisbg=axcolor)
        axarea = plt.axes([0.25, 0.25, 0.65, 0.03], axisbg=axcolor)

        sm = Slider(axm, 'slope_delta', -np.pi/4, np.pi/4, valinit=0)
        sc = Slider(axc, 'interception_delta', -(max(self.y) - min(self.y)) / 4,
                    (max(self.y) - min(self.y)) / 4, valinit=0)
        salpha = Slider(axalpha, 'alpha', 0, 0.8, valinit=self.alpha0)
        sarea = Slider(axarea, 'area', 1, 50, valinit=self.area0)

        def update_m(val):
            """
            :returns: TODO
            """
            x_axis = np.median(self.x)
            y_axis = x_axis * self.m + self.c
            self._updateM(sm.val)
            self.c = y_axis - self.m * x_axis
            # fl.set_ydata(self.m*xseq + self.c)
            fig.canvas.draw_idle()

        def update_c(val):
            """
            :returns: TODO
            """
            self._updateC(sc.val)
            # fl.set_ydata(self.m*xseq + self.c)
            fig.canvas.draw_idle()

        def update_alpha(val):
            """
            :returns: TODO
            """
            self.alpha = salpha.val
            pts.set_alpha(self.alpha)
            fig.canvas.draw_idle()

        def update_area(val):
            """
            :returns: TODO
            """
            self.area = sarea.val
            pts.set_sizes(np.ones(len(self.x)) * self.area)
            fig.canvas.draw_idle()

        sm.on_changed(update_m)
        sc.on_changed(update_c)
        salpha.on_changed(update_alpha)
        sarea.on_changed(update_area)

        select_col = lambda array, colNum: map(lambda x: x[colNum], array)

        def onselect(verts):
            path = Path(verts)
            xys = pts.get_offsets()
            ind = np.nonzero([path.contains_point(xy) for xy in xys])[0]
            self.x = np.array(select_col(xys[ind], 0))
            self.y = np.array(select_col(xys[ind], 1))

            # self.colorin = colorConverter.to_rgba('red', self.alpha0)
            # self.colorout = colorConverter.to_rgba('blue', self.alpha0)
            # facecolors = np.array([self.colorout for xy in xys])
            # facecolors[ind] = self.colorin
            # pts.set_facecolors(facecolors)

            A = np.vstack([self.x, np.ones(len(self.x))]).T
            self.m, self.c = np.linalg.lstsq(A, self.y)[0]
            self.m1, self.c1 = self.m, self.c
            # fl.set_ydata(self.m*xseq + self.c)
            # hl.set_ydata(np.median(self.y))

            sm.reset()
            sc.reset()

            fig.canvas.draw_idle()

        lasso = LassoSelector(ax, onselect)

        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

        def reset(event):
            sm.reset()
            sc.reset()
            # self.colorout = colorConverter.to_rgba('blue', self.alpha0)
            # facecolors = np.array([self.colorout for y in y0])
            # pts.set_facecolors(facecolors)
            self.y = y0
            self.x = x0
            A = np.vstack([self.x, np.ones(len(self.x))]).T
            self.m, self.c = np.linalg.lstsq(A, self.y)[0]
            self.m1, self.c1 = self.m, self.c
            # fl.set_ydata(self.m*xseq + self.c)
            # hl.set_ydata(np.median(self.y))

        button.on_clicked(reset)

        resetax = plt.axes([0.5, 0.025, 0.1, 0.04])
        button_exit = Button(resetax, 'Ok', color=axcolor, hovercolor='0.975')

        def ok_exit(event):
            plt.close()

        button_exit.on_clicked(ok_exit)

        plt.show()


class BAFPlot(object):

    """Plot BAF hist and density curve"""

    def __init__(self, BAF):
        self._BAF = BAF

    def plot_density(self, bandwidth):
        print self._BAF
        sns.set_style('whitegrid')
        sns.kdeplot(self._BAF, bw=bandwidth)

    def plot_hist(self, bin=0.01):
        print self._BAF
        plt.hist(self._BAF, bins='auto')
        plt.show()

def plotMeanShift(sid, X, labels, cluster_centers):

    # The function to write the data used in
    # plotting to a file.
    # sid : the id for the stripe
    # X: the X axis used in plotting
    # k_List: k value associated with the X axis
    def writeToFile(sid, X_list, k_List):
        if not os.path.exists("plot_data/BAF_plot_data"):
            os.makedirs("plot_data/BAF_plot_data")
        if not os.path.isfile("plot_data/BAF_plot_data/BAFPlot_data.txt"):
            outputFile = open("plot_data/BAF_plot_data/BAFPlot_data.txt", "wr")
            outputFile.write("#{0}\t{1}\t{2}\n".format("BAF", "class", "stripe"))
        else:
            outputFile = open("plot_data/BAF_plot_data/BAFPlot_data.txt", "a")
        for i in range(0, len(k_List)):
            for j in range(0, len(X_list[i])):
                outputFile.write("{0}\t{1}\t{2}\n".format(X_list[i][j], k_List[i], sid))
        outputFile.close()

    def writeCenterToFile(sid,cluster_center_list):
        if not os.path.exists("plot_data/BAF_plot_data"):
            os.makedirs("plot_data/BAF_plot_data")
        if not os.path.isfile("plot_data/BAF_plot_data/BAFPlot_centers_data.txt"):
            outputFile = open("plot_data/BAF_plot_data/BAFPlot_centers_data.txt","wr")
            outputFile.write("#cluster_centers\tstripe\n")
        else:
            outputFile = open("plot_data/BAF_plot_data/BAFPlot_centers_data.txt", "a")
        for i in range(0,len(cluster_center_list)):
            outputFile.write("{0}\t{1}\n".format(cluster_center_list[i],sid))

        outputFile.close()

    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)

    print("number of estimated clusters : %d" % n_clusters_)

    plt.figure(1)
    plt.clf()

    #the following three parameters are for write to file only
    X_list = [] # the x coordinate that get's plotted
    k_list = [] # the corresponding k for every X get's plotted
    cluster_center_list = []

    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    for k, col in zip(range(n_clusters_), colors):
        my_members = labels == k
        cluster_center = cluster_centers[k]
        plt.plot(X[my_members, 0], X[my_members, 1], col + '.')
        plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=14)
        X_list.append(X[my_members, 0])
        k_list.append(k)
        cluster_center_list.append(cluster_center[0])

    writeToFile(sid,X_list,k_list)
    writeCenterToFile(sid,cluster_center_list)
    plt.title(str(sid) + ' estimated number of clusters: %d' % n_clusters_)
    plt.xlabel("BAF")
    plt.show()




def main():
    """
    :returns: TODO

    """
    from data import Segment

    segments = []

    num = 200000

    for i in range(num):
        seg = Segment()
        seg.tReadNum = np.random.randint(200, 300)
        seg.nReadNum = np.random.randint(300, 400)
        seg.gc = np.random.rand()
        segments.append(seg)

    gsp = GCStripePlot(segments, 20000)
    gsp.plot()
    print gsp.output()


if __name__ == "__main__":
    main()

'''
# =============================================================================
#      FileName: stripe.py
#          Desc: Stripe aggregation and decomposition
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-11-24 13:53:39
#       History:
# =============================================================================
'''

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import heapq
import sys
from collections import Counter

import numpy as np
from scipy.cluster import hierarchy
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
from scipy.stats.mstats import gmean
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.datasets.samples_generator import make_blobs

import constants


class Stripe:
    def __init__(self):
        # 记录对应原始seg的索引
        self.index_seg = None

        self.paired_counts = None

        self.rdr = -1.0

    def init_seg(self, seg_list, seg_idx):
        self.index_seg = seg_idx

        self._getRD(seg_list)
        self._getBAF(seg_list)

    def _getRD(self, seg_list):
        # 获取几何平均值
        ratios = [
            seg.tumor_reads_num * 1.0 / seg.normal_reads_num
            for seg in seg_list
        ]
        self.rdr = gmean(ratios)

    def _getBAF(self, seg_list):
        self.paired_counts = np.array(
            [[], [], [], [], [], []], dtype=int).transpose()

        for seg in seg_list:
            self.paired_counts = np.vstack((self.paired_counts,
                                            seg.paired_counts))

    def log_likelihood(self, phi, tp, update_tree=True, new_state=0):
        if update_tree:
            ##################################################
            # some useful info about the tree,
            # used by CNV related computations,
            u.set_node_height(self.tssb)
            u.set_path_from_root_to_node(self.tssb)
            u.map_datum_to_node(self.tssb)
            ##################################################
        return self.__log_complete_likelihood__(
            phi, self.mu_r, self.mu_v, tp, new_state)


    def __log_complete_likelihood__(self, phi, mu_r, mu_v, tp, new_state=0):

        if self.cnv:
            ll = []
            poss_n_genomes = self.compute_n_genomes(tp, new_state)
            poss_n_genomes = [x for x in poss_n_genomes if x[1] > 0]
            for (nr, nv) in poss_n_genomes:
                mu = (nr * mu_r + nv*(1-mu_r)) / (nr + nv)
                ll.append(
                    u.log_binomial_likelihood(self.a[tp],
                                              self.d[tp],
                                              mu) +
                    log(1.0 / len(poss_n_genomes)) + self._log_bin_norm_const
                    [tp])
            if len(poss_n_genomes) == 0:
                ll.append(log(1e-99))  # to handle cases with zero likelihood
            llh = u.logsumexp(ll)
        else:  # CNV datum
            mu = (1 - phi) * mu_r + phi*mu_v  # (mu_r=0.999, mu_v=0.5)
            llh = u.log_binomial_likelihood(
                self.a[tp], self.d[tp], mu) + self._log_bin_norm_const[tp]
        return llh


class MergeSeg(object):
    """The stripe objects, including load, property operations"""

    def __init__(self, data):
        """import data object

        :data: TODO

        """
        self._data = data
        self.stripes = []  # stripes

    def get(self):
        """TODO: Docstring for get.
        :returns: TODO

        """
        self._aggregation(y_down, y_up, stripe_num, noise_stripe_num=2)

    def _aggregation(self, y_down, y_up, stripe_num, noise_stripe_num=2):
        """The aggregation operations for segments in data

        :returns: stripes data structure

        """
        assert stripe_num > 0

        reads_depth_ratio_log = []

        # here should keep index
        yc_v = np.array([
            np.log(seg.tumor_reads_num + 1) - np.log(seg.normal_reads_num + 1)
            for seg in self._data.segments
        ])

        # 记录是否是outlier
        status_yc_v = np.logical_and(yc_v > y_min, yc_v < y_max)

        y_fcd = yc_v.reshape(yc_v.shape[0], 1)
        clusters = hierarchy.fclusterdata(
            y_fcd, stripe_num + noise_stripe_num, criterion="distance")

        # 此处应该只获取最大和最小值之间的条带，且要保留原始位置，以方便索引
        # 此处获取最小和最大值之间的条带的方法是：直接去除这些位置不列入计算范围

        # 此处应该是去除了outlier之后的Counter

        mccs = Counter(
            clusters[status_yc_v]).most_common(stripe_num + noise_stripe_num)

        for c_id, _ in mccs:
            # 对每一个条带进行裂解操作，生成子条带, return
            self._decomposition(c_id, clusters, status_yc_v)

    def _decomposition(self, c_id, clusters, status_yc_v):
        """The decomposition operations for segments in data

        :parameters: TODO
        :returns: TODO

        """
        # 获得该类别的所有结点idx：
        # 即，clusters 中与c_id相等且，在status_yc_v中的位置
        ca = np.argwhere(clusters == c_id).flatten()
        sa = np.argwhere(status_yc_v).flatten()
        mstrip_seg_idx = np.intersectid(ca, sa)

        # 这里的基于BAF的归类处理分为3个步骤

        # 首先进行所有seg的BAF的密度估计，然后获得峰值    类别定位
        # 然后对每一个seg进行归类，按照内部投票的方式     Seg归类
        # 然后返回

        # 这里需要有一个记录原始向量中位置的向量
        seg_list = [self._data.segments[idx] for idx in mstrip_seg_idx]

        paired_counts_all = np.array(
            [[], [], [], [], [], []], dtype=int).transpose()
        for seg in seg_list:
            paired_counts_all = np.vstack((paired_counts_all,
                                           seg.paired_counts))

        a_T = paired_counts_all[:, 2]
        b_T = paired_counts_all[:, 3]
        d_T = a_T + b_T
        l_T = np.min(paired_counts_all[:, 2:4], axis=1)
        p_T = l_T * 1.0 / d_T

        # status_p_T_v = np.logical_and(p_T > p_T_min, p_T < p_T_max).flatten()

        y = np.ones(p_T.shape)
        p_T_y = np.hstack((p_T, y))
        bandwidth = estimate_bandwidth(p_T_y, quantile=0.2, n_samples=500)
        ms.fit(X)
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_
        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)

        seg_label = [
            self._getSegLabl(seg, cluster_centers) for seg in seg_list
        ]

        for label in set(seg_label):
            if label == -1:
                continue
            sub_seg_list = [
                seg for seg, idx in enumerate(seg_list)
                if seg_label[idx] == label
            ]
            sub_seg_idx = [
                mstrip_seg_idx[idx] for seg, idx in enumerate(seg_list)
                if seg_label[idx] == label
            ]
            temp_stripe = Stripe()
            temp_stripe.init_seg(sub_seg_list, sub_seg_idx)
            self.stripes.append(temp_stripe)

    def _getSegLabl(self, seg, cluster_centers):
        if seg.paired_counts is None:
            return -1

        a_T_seg = seg.paired_counts[:, 2]
        b_T_seg = seg.paired_counts[:, 3]
        d_T_seg = a_T_seg + b_T_seg
        l_T_seg = np.min(seg.paired_counts[:, 2:4], axis=1)
        p_T_seg = l_T_seg * 1.0 / d_T_seg

        dis_seg = np.abs(p_T_seg[:, None] - cluster_centers[:, 0])
        labels_seg = np.argmin(dis_seg, axis=1)

        return Counter(labels_seg).most_common(1)[0][0]

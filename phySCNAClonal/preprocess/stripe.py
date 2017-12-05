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
from random import randint

import numpy as np
from scipy.cluster import hierarchy
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
from scipy.stats.mstats import gmean
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.datasets.samples_generator import make_blobs

import constants
from pydp.densities import Density, log_poisson_pdf
from utils import (get_cn_allele_config, get_loga, get_mu_E_joint,
                   log_binomial_likelihood, mad_based_outlier)


class Stripe:
    def __init__(self):
        # 记录对应原始seg的索引
        self.index_seg = None

        self.paired_counts = None

        self.tumor_reads_num = -1.0
        self.normal_reads_num = -1.0
        self.rdr = -1.0

        self.baseline_label = False
        self._baseline = -1

    def init_seg(self, seg_list, seg_idx):
        self.index_seg = seg_idx

        self._init_RD(seg_list)
        self._init_BAF(seg_list)

    def _init_RD(self, seg_list):
        # 获取几何平均值
        tumor_reads_num = [seg.tumor_reads_num for seg in seg_list]
        normal_reads_num = [seg.normal_reads_num for seg in seg_list]

        self.tumor_reads_num = gmean(tumor_reads_num)
        self.normal_reads_num = gmean(normal_reads_num)

        ratios = [
            seg.tumor_reads_num * 1.0 / seg.normal_reads_num
            for seg in seg_list
        ]
        self.rdr = gmean(ratios)

    def _init_BAF(self, seg_list):
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


    def __log_likelihood_CN_BAF(self, phi):

        # if seg.phi_last is not None and Phi.phi in seg.phi_last.keys(): # return seg.phi_last[Phi.phi].likelihood
        # else:
            # if seg.phi_last is None:
                # seg.phi_last = {}

            # if len(seg.phi_last) >= constants.SEG_RECORD_SIZE:
                # idx = randint(0, constants.SEG_RECORD_SIZE - 1)
                # seg.phi_last.pop(seg.phi_last.keys()[idx])

        # 此处是否添加记录
        ll, cn, pi = self._getStripeCNBAF(phi)
        # seg.phi_last[Phi.phi] = SegmentResultData(ll, cn, pi)
        return ll


    def _getStripeCNBAF(self, phi):
        max_copy_number = constants.MAX_COPY_NUMBER

        cns = None
        if self.baseline_label == "True":
            cns = [2]
        elif self.rdr > self._baseline:
            cns = range(2, max_copy_number + 1)
        else:
            cns = range(0, 2 + 1)

        ll_pi_s = [self._getLLStripe(cn, phi) for cn in cns]
        (ll, pi) = max(ll_pi_s, key=lambda x: x[0])
        cn = cns[ll_pi_s.index((ll, pi))]
        return ll, cn, pi

    def _getLLStripe(self, cn, phi):
        allele_config = get_cn_allele_config(contants.MAX_COPY_NUMBER)
        rd_weight = constants.RD_WEIGHT

        ll_seg = 0
        ll_rd = self._getRD(cn, phi)

        allele_types = allele_config[cn]
        if 0 == seg.paired_counts.shape[0]:
            ll_baf = 0
            pi = "*"
        else:
            ll_baf, pi = self._getBAF(cn, allele_types, phi)
        ll_seg = ll_baf + ll_rd
        return ll_seg, pi

    def _getRD(self, cn, phi):
        c_N = constants.COPY_NUMBER_NORMAL

        bar_c = phi * cn + (1.0 - phi) * c_N

        lambda_possion = (
            bar_c / c_N) * self._baseline * (self.normal_reads_num + 1) #not minus 1 ? better
        if lambda_possion < 0:
            lambda_possion = 0

        ll_RD = log_poisson_pdf(self.tumor_reads_num, lambda_possion)
        return ll_RD

    def _getBAF(self, seg, cn, allele_types, phi):
        c_N = constants.COPY_NUMBER_NORMAL
        mu_N = constants.MU_N

        # keys, ppmm values 0.5
        mu_G = np.array(allele_types.values())

        mu_E = get_mu_E_joint(mu_N, mu_G, c_N, cn, phi)

        if self.paired_counts.shape[0] > 1:
            b_T_j = np.min(self.paired_counts[:, 2:4], axis=1)
            d_T_j = np.sum(self.paired_counts[:, 2:4], axis=1)
            baf = b_T_j * 1.0 / d_T_j
            outlier = mad_based_outlier(baf)
            BAF = np.delete(self.paired_counts, list(outlier.astype(int)), axis=0)
            b_T_j = np.min(BAF[:, 2:4], axis=1)
            d_T_j = np.sum(BAF[:, 2:4], axis=1)

        else:
            b_T_j = np.min(self.paired_counts[:, 2:4], axis=1)
            d_T_j = np.sum(self.paired_counts[:, 2:4], axis=1)
            pass

        # add prior or not?
        ll = log_binomial_likelihood(b_T_j, d_T_j, mu_E)
        ll_bafs = ll.sum(axis=0)
        idx_max = ll_bafs.argmax(axis=0)
        ll_baf = ll_bafs[idx_max]
        pi = allele_types[allele_types.keys()[idx_max]]
        return ll_baf, pi


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

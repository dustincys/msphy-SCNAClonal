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

# import json

class Stripe:
    def __init__(self):
        # 记录对应原始seg的索引
        self.stripe_name = ""
        self.stripe_id = ""
        self.segs_idx = None

        self.paired_counts = None

        self.tumor_reads_num = -1.0
        self.normal_reads_num = -1.0

        # 似乎用不到，实际上就是:
        # self.tumor_reads_num /self.normal_reads_num
        # self.rdr = -1.0

        self.baseline_label = False
        # 似乎不应该放在这里
        self._baseline = -1

        # 这两个应该放在这里
        self.copy_number = -1
        self.genotype = ""

        # phi应该放在node结点中
        # self.phi = 0.0

        self.tssb = None
        self.node = None  # this is the node where the datum resides

    def init_segs(self, segs_list, segs_idx):
        self.segs_idx = segs_idx

        self._init_RD(segs_list)
        self._init_BAF(segs_list)

    def _init_RD(self, segs_list):
        # 获取几何平均值
        tumor_reads_num = [seg.tumor_reads_num for seg in segs_list]
        normal_reads_num = [seg.normal_reads_num for seg in segs_list]

        self.tumor_reads_num = gmean(tumor_reads_num)
        self.normal_reads_num = gmean(normal_reads_num)

    def _init_BAF(self, segs_list):
        self.paired_counts = np.array(
            [[], [], [], [], [], []], dtype=int).transpose()

        for seg in segs_list:
            self.paired_counts = np.vstack((self.paired_counts,
                                            seg.paired_counts))

    def _log_likelihood(self, phi, update_tree=True):
        if update_tree:
            ##################################################
            # some useful info about the tree,
            # used by CNV related computations,
            u.set_node_height(self.tssb)
            u.set_path_from_root_to_node(self.tssb)
            u.map_datum_to_node(self.tssb)
            ##################################################

        # 注意：此处应该不受CN\genotype的限制，但不记录目标的CN和genotype
        # 因为此处的parameter不是准确的phi,所以无法达到最优，但为了能够抽样得到
        # 最佳结构。此处要使CN和genotype自由发挥

        # 在时间上的先后顺序能够明确影响测序数据的位置，只有重叠位置的时候才会
        # 发生

        return self.__log_likelihood_RD_BAF(phi)


    def __log_likelihood_RD_BAF(self, phi):
        # 此处是否添加记录
        copy_numbers = None
        if seg.baseline_label == "True":
            copy_numbers = [2]
        elif get_loga(seg) > self._baseline:
            copy_numbers = range(2, self._max_copy_number + 1)
        else:
            copy_numbers = range(0, 2 + 1)

        ll_pi_s = [self._getLLStripe(copy_number, phi) for copy_number in
                   copy_numbers]
        ll, pi = max(ll_pi_s, key=lambda x: x[0])
        cn = ll_pi_s.index((ll, pi))

        self.copy_number = cn
        self.genotype = pi

        return ll


    def _getLLStripe(self, copy_number, phi):
        rd_weight = constants.RD_WEIGHT_TSSB

        ll_stripe = 0
        ll_rd = self._getRD(copy_number, phi)
        allele_types = self._allele_config[copy_number]

        # remove the weak baf point
        self._augBAF(copy_number)

        if 0 == self.paired_counts.shape[0]:
            ll_baf = 0
            pi = "*"
        else:
            ll_baf, pi = self._getBAF(self, copy_number, allele_types, phi)

        ll_stripe = ll_rd * rd_weight + ll_baf * (1 - rd_weight)

        return ll_stripe, pi


    def _augBAF(self, copy_number):
        # todo: remove the baf point is not a good idea
        threshold = constants.BAF_THRESHOLD * self._coverage

        if copy_number > 2:
            d_T_j = np.sum(self.paired_counts[:, 2:4], axis=1)
            idx_rm = tuple(np.where(d_T_j < threshold)[0])
            self.paired_counts = np.delete(self.paired_counts, idx_rm, axis=0)
        else:
            pass

    def _getRD(self, copy_number, phi):
        c_N = constants.COPY_NUMBER_NORMAL

        bar_c = phi * copy_number + (1.0 - phi) * c_N

        lambda_possion = (
            bar_c / c_N) * self._baseline * (seg.normal_reads_num + 1)
        if lambda_possion < 0:
            lambda_possion = 0

        ll_RD = log_poisson_pdf(seg.tumor_reads_num, lambda_possion)
        return ll_RD

    def _getBAF(self, copy_number, allele_types, phi):
        c_N = constants.COPY_NUMBER_NORMAL
        mu_N = constants.MU_N

        mu_G = np.array(allele_types.values())
        mu_E = get_mu_E_joint(mu_N, mu_G, c_N, copy_number, phi)

        if seg.paired_counts.shape[0] > 1:
            b_T_j = np.min(seg.paired_counts[:, 2:4], axis=1)
            d_T_j = np.sum(seg.paired_counts[:, 2:4], axis=1)
            baf = b_T_j * 1.0 / d_T_j
            outlier = mad_based_outlier(baf)
            BAF = np.delete(seg.paired_counts, list(outlier.astype(int)), axis=0)
            b_T_j = np.min(BAF[:, 2:4], axis=1)
            d_T_j = np.sum(BAF[:, 2:4], axis=1)

        else:
            b_T_j = np.min(seg.paired_counts[:, 2:4], axis=1)
            d_T_j = np.sum(seg.paired_counts[:, 2:4], axis=1)
            pass

        ll = log_binomial_likelihood(b_T_j, d_T_j, mu_E)
        ll_bafs = ll.sum(axis=0)
        idx_max = ll_bafs.argmax(axis=0)
        ll_baf = ll_bafs[idx_max]
        pi = allele_types[allele_types.keys()[idx_max]]
        return ll_baf, pi


class DataStripes(object):
    """The stripe objects, including load, property operations"""

    def __init__(self, data):
        """import data object

        :data: TODO

        """
        self._data = data
        self.stripes = []  # stripes

        self.baseline = -1

    def get(self):
        """TODO: Docstring for get.
        :returns: TODO

        """
        self._aggregation(y_down, y_up, stripe_num, noise_stripe_num=2)

    def output_txt(self, outFileName):
        with open(outFileName, 'w') as outFile:
            outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                "id", "segs_idx", "paired_counts", "tumor_reads_num",
                "normal_reads_num", "baseline_label"))

            for s in self.stripes:
                a_T = s.paired_counts[:,2]
                b_T = s.paired_counts[:,3]
                a_T_strl = np.array_str(a_T).strip("[]").split()
                b_T_strl = np.array_str(b_T).strip("[]").split()

                outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    s.stripe_id,
                    ",".join(s.segs_idx),
                    "{0}|{1}".format(",".join(a_T_strl), ",".join(b_T_strl)),
                    s.tumor_reads_num,
                    s.normal_reads_num,
                    s.baseline_label)
            pass

    def _aggregation(self, y_down, y_up, stripe_num, noise_stripe_num=2):
        """The aggregation operations for segments in data

        :returns: stripes data structure

        """
        assert stripe_num > 0

        reads_depth_ratio_log = []

        # here should keep idx
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
        segs_list = [self._data.segments[idx] for idx in mstrip_seg_idx]

        paired_counts_all = np.array(
            [[], [], [], [], [], []], dtype=int).transpose()
        for seg in segs_list:
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
            self._getSegLabl(seg, cluster_centers) for seg in segs_list
        ]

        for label in set(seg_label):
            if label == -1:
                continue
            sub_seg_list = [
                seg for seg, idx in enumerate(segs_list)
                if seg_label[idx] == label
            ]
            sub_seg_idx = [
                mstrip_seg_idx[idx] for seg, idx in enumerate(segs_list)
                if seg_label[idx] == label
            ]
            temp_stripe = Stripe()
            temp_stripe.stripe_id = "{0}_{1}".format(str(c_id), str(idx))
            temp_stripe.init_segs(sub_seg_list, sub_seg_idx)
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

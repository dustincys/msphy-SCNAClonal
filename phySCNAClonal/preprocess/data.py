'''
# =============================================================================
#      FileName: data.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-11-24 13:53:26
#       History:
# =============================================================================
'''
import sys
import numpy as np
import scipy.cluster.hierarchy as hcluster

from segment import Segment
from scipy.stats.mstats import gmean
from collections import Counter


class Data:

    self.seg_num = 0
    self.baseline = -1

    self.segments = []

#       peak reange
    self.pr = 0

    def __init__(self, max_copynumber = 6, coverage = 30):

        self.max_copynumber = max_copynumber
        self.coverage = coverage

        self.allele_config = get_cn_allele_config(max_copy_number)


    def load_segments_bed(self, bed_file_normalized_name):
        """

        :bed_file_normalized_name: TODO
        :returns: TODO

        """

        bed_chroms, bed_starts, bed_ends, tumor_reads, normal_reads, gcs =\
            BEDnParser(bed_file_normalized_name)
        get_chrom_format(bed_chroms)
        bed_num = len(bed_chroms)

        for i in range(0, bed_num):
            chrom_idx = chrom_name_to_idx(bed_chroms[i])
            seg_name = get_segment_name(
                bed_chroms[i], bed_starts[i], bed_ends[i])

            normal_reads_num = normal_reads[i]
            tumor_reads_num = tumor_reads[i]

            segment_i = Segment()
            segment_i.name = seg_name
            segment_i.chrom_idx = chrom_idx
            segment_i.chrom_name = bed_chroms[i]
            segment_i.start = bed_starts[i]
            segment_i.end = bed_ends[i]
            segment_i.normal_reads_num = normal_reads_num
            segment_i.tumor_reads_num = tumor_reads_num

            if 0 == normal_reads_num:
                segment_i.log2_ratio = -float('Inf')
            else:
                segment_i.log2_ratio = np.log2(1.0 *
                                               tumor_reads_num/normal_reads_num)

            segment_i.log_ratio = np.log(1.0 * (tumor_reads_num + 1.0) /
                                         (normal_reads_num + 1.0))
            segment_i.gc = gcs[i]

            self.segments.append(segment_i)
            self.seg_num += 1

    def load_segments(self, normal_bam, tumor_bam, bed_file_name):
        chrom_idx_list = constants.CHROM_IDX_LIST
        chrom_start = constants.CHROM_START

        sam_SQ = normal_bam.header['SQ']
        sam_chrom_format = get_chrom_format(map(lambda x: x['SN'], sam_SQ))
        chrom_lens, chrom_idxs = get_chrom_lens_idxs(chrom_idx_list, sam_SQ)

        bed_chroms, bed_starts, bed_ends, gcs = BEDParser(bed_file_name)
        get_chrom_format(bed_chroms)
        bed_num = len(bed_chroms)

        for i in range(0, bed_num):
            chrom_idx = chrom_name_to_idx(bed_chroms[i])
            chrom_name = chrom_idx_to_name(chrom_idx, sam_chrom_format)
            seg_name = get_segment_name(chrom_name, bed_starts[i], bed_ends[i])

            if chrom_idx not in chrom_idx_list:
                print 'Chromsome {0} not found, segment {1} excluded...'.format(
                    bed_chroms[i], seg_name)
                sys.stdout.flush()
                continue

            chrom_lst_idx = chrom_idxs.index(chrom_idx)

            if bed_starts[i] < chrom_start or bed_ends[
                    i] > chrom_lens[chrom_lst_idx]:
                print 'Out of range chromsome {0}, segment {1} excluded...'.\
                    format(bed_chroms[i], seg_name)
                sys.stdout.flush()
                continue

            normal_reads_num = normal_bam.count(
                chrom_name, bed_starts[i], bed_ends[i])
            tumor_reads_num = tumor_bam.count(
                chrom_name, bed_starts[i], bed_ends[i])

            segment_i = Segment()
            segment_i.name = seg_name
            segment_i.chrom_idx = chrom_idx
            segment_i.chrom_name = chrom_name
            segment_i.start = bed_starts[i]
            segment_i.end = bed_ends[i]
            segment_i.normal_reads_num = normal_reads_num
            segment_i.tumor_reads_num = tumor_reads_num
            segment_i.log2_ratio = np.log2(1.0 *
                                           tumor_reads_num/normal_reads_num)

            segment_i.log_ratio = np.log(1.0 * (tumor_reads_num + 1.0) /
                                         (normal_reads_num + 1.0))
            self.segments.append(segment_i)
            self.seg_num += 1

    def get_LOH_frac(self):
        for j in range(0, self.seg_num):
            self.segments[j].LOH_frac = get_LOH_frac(
                self.segments[j].paired_counts)

    def get_APM_frac(self):
        """
        :returns: TODO

        """
        for j in range(0, self.seg_num):
            self.segments[j].APM_frac = get_APM_frac_MAXMIN(
                self.segments[j].paired_counts)

    def get_LOH_status(self, baseline_thred, flag_runpreprocess=False):

        if flag_runpreprocess:
            LOH_num = 0
            FLOH_num = 0
            for j in range(0, self.seg_num):
                self.segments[j].LOH_status = get_LOH_status(
                    self.segments[j].LOH_frac, baseline_thred)
                if self.segments[j].LOH_status == "TRUE":
                    LOH_num = LOH_num + 1
                elif self.segments[j].LOH_status == "FALSE":
                    FLOH_num = FLOH_num + 1

            print "LOH_num/seg_num = {0}/{1}".format(LOH_num, self.seg_num)
            print "FLOH_num/seg_num = {0}/{1}".format(FLOH_num, self.seg_num)
        else:
            print "get_LOH_status function called from model."

    def get_APM_status(self, baseline_thred_APM):
        APM_num = 0
        for j in range(0, self.seg_num):
            self.segments[j].APM_status = get_APM_status(
                self.segments[j].APM_frac, baseline_thred_APM)
            if self.segments[j].APM_status == "TRUE":
                APM_num = APM_num + 1

        print "APM_num/seg_num = {0}/{1}".format(APM_num, self.seg_num)

    def compute_Lambda_S_LOH(self, max_copynumber, subclone_num,
                             flag_runpreprocess=False):
        """ compute the Lambda S, through hierarchy clustering
        """
        if not flag_runpreprocess:
            print "compute_Lambda_S function called from model"
            sys.stdout.flush()
            return

        thresh = constants.HC_THRESH

        reads_depth_ratio_log = []
        reads_depth_ratio = []
        for j in range(0, self.seg_num):
            if self.segments[j].APM_status == 'TRUE' and\
                    self.segments[j].LOH_status == 'FALSE':
                ratio = self.segments[
                    j].tumor_reads_num*1.0/self.segments[j].normal_reads_num
                reads_depth_ratio_log.append(np.log(ratio))
                reads_depth_ratio.append(ratio)

        reads_depth_ratio = np.array(reads_depth_ratio)
        reads_depth_ratio_log = np.array(reads_depth_ratio_log)
        if reads_depth_ratio_log.shape[0] == 0:
            print 'Error: no APM-LOH position found, existing...'
            print 'Either the baseline_thred_APM is too large, or the constants\
            APM_N_MIN is too large; Or, the baseline_thred_LOH is too small'
            sys.exit(1)

        reads_depth_ratio_log = reads_depth_ratio_log.reshape(
            reads_depth_ratio_log.shape[0], 1)
        y = np.ones(reads_depth_ratio_log.shape)
        reads_depth_ratio_log = np.hstack((reads_depth_ratio_log, y))
        clusters = hcluster.fclusterdata(
            reads_depth_ratio_log, thresh, criterion="distance")
        mccs = Counter(clusters).most_common(max_copynumber * subclone_num)

        rdr_min = float('Inf')
        cluster_min = -1
        for i in range(0, len(mccs)):
            cluster_temp = mccs[i][0]
            print "cluster temp : {}".format(cluster_temp)
            rdr_temp = gmean(reads_depth_ratio[clusters == cluster_temp])
            print "rdr_temp"
            print "log: {}".format(np.log(rdr_temp))
            if rdr_min > rdr_temp:
                rdr_min = rdr_temp
                cluster_min = cluster_temp

        print mccs
        print "log baseline: {}".format(np.log(rdr_min))
        sys.stdout.flush()

        cluster_flag = (clusters == cluster_min)
        baseline_num = 0
        rdr_i = 0
        for j in range(0, self.seg_num):
            if self.segments[j].APM_status == 'TRUE' and\
                    self.segments[j].LOH_status == 'FALSE':
                if cluster_flag[rdr_i]:
                    self.segments[j].baseline_label = 'TRUE'
                    baseline_num = baseline_num + 1
                else:
                    self.segments[j].baseline_label = 'FALSE'
                rdr_i = rdr_i + 1
            else:
                self.segments[j].baseline_label = 'FALSE'

        print "baseline_num: {}".format(baseline_num)

        if baseline_num == 0:
            print 'Error: No diploid segments found, existing...'
            sys.exit(1)

        self.baseline = rdr_min

    def _log_likelihood(self, id, phi, update_tree=True, new_state=0):

        if update_tree:
            ##################################################
            # some useful info about the tree,
            # used by CNV related computations,
            u.set_node_height(self.tssb)
            u.set_path_from_root_to_node(self.tssb)
            u.map_datum_to_node(self.tssb)
            ##################################################

        seg = self.segments[id]
        ll, cn, pi = self._getSegResData(seg, phi)

        return ll

    def _getSegResData(self, seg, phi):
        copy_numbers = None
        if seg.baseline_label == "True":
            copy_numbers = [2]
        elif get_loga(seg) > self.baseline:
            copy_numbers = range(2, self.max_copy_number + 1)
        else:
            copy_numbers = range(0, 2 + 1)

        ll_pi_s = [self._getLLSeg(seg, copy_number, phi) for copy_number in
                   copy_numbers]
        (ll, pi) = max(ll_pi_s, key=lambda x: x[0])
        cn = ll_pi_s.index((ll, pi))
        return ll, cn, pi

    def _getLLSeg(self, seg, copy_number, phi):
        ll_seg = 0
        ll_rd = self._getRD(seg, copy_number, phi)
        allele_types = self._allele_config[copy_number]
        self._augBAF(seg, copy_number)
        if 0 == seg.paired_counts.shape[0]:
            ll_baf = 0
            pi = "*"
        else:
            ll_baf, pi = self._getBAF(seg, copy_number, allele_types, phi)
        ll_seg = ll_baf + ll_rd
        return ll_seg, pi

    def _augBAF(self, seg, copy_number):
        if copy_number > 2:
            threshold = constants.BAF_THRESHOLD * self._coverage
            d_T_j = np.sum(seg.BAF[:, 2:4], axis=1)
            idx_rm = tuple(np.where(d_T_j < threshold)[0])
            seg.BAF = np.delete(seg.BAF, idx_rm, axis=0)
        else:
            pass

    def _getRD(self, seg, copy_number, phi):
        c_N = constants.COPY_NUMBER_NORMAL
        bar_c = phi * copy_number + (1.0 - phi) * c_N
        print "____>>> _getRD: bar_c, c_N, self._baseline, seg.normal_reads_num____"
        print bar_c, c_N, self._baseline, seg.normal_reads_num
        print "_________end _getRD:bar_c, c_N, self._baseline, seg.normal_reads_num______________"

        lambda_possion = (
            bar_c / c_N) * self._baseline * (seg.normal_reads_num + 1) #not minus 1 ? better
        if lambda_possion < 0:
            lambda_possion = 0
        print "____>>> _getRD: seg.tumor_reads_num, lambda_possion____"
        print seg.tumor_reads_num, lambda_possion
        print "_________end _getRD:seg.tumor_reads_num, lambda_possion______________"

        ll_RD = log_poisson_pdf(seg.tumor_reads_num, lambda_possion)
        return ll_RD

    def _getBAF(self, seg, copy_number, allele_types, phi):
        c_N = constants.COPY_NUMBER_NORMAL
        mu_N = constants.MU_N
        # keys, ppmm values 0.5
        mu_G = np.array(allele_types.values())

        print "____>>> _getBAF: mu_N, mu_G, c_N, copy_number, phi____"
        print mu_N, mu_G, c_N, copy_number, phi
        print "_________end _getBAF:mu_N, mu_G, c_N, copy_number, phi______________"

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

        # add prior or not?
        ll = log_binomial_likelihood(b_T_j, d_T_j, mu_E)
        ll_bafs = ll.sum(axis=0)
        idx_max = ll_bafs.argmax(axis=0)
        ll_baf = ll_bafs[idx_max]
        pi = allele_types[allele_types.keys()[idx_max]]
        return ll_baf, pi

    # computes the binomial parameter
    def compute_n_genomes(self, tp, new_state=0):
        def descend(nd, new_state):
            # this is needed for Metropolis-Hastings likelihood computations
            pi = nd.pi1[tp] if new_state else nd.pi[tp]
            ssm_node = self.node.path[-1]
            mr_cnv = self.find_most_recent_cnv(nd)
            ancestors = nd.get_ancestors()
            if (ssm_node not in ancestors) and (not mr_cnv):
                self.nr1 += pi * 2
                self.nr2 += pi * 2
                self.nr3 += pi * 2
                self.nr4 += pi * 2
            elif ssm_node in ancestors and (not mr_cnv):
                self.nr1 += pi
                self.nv1 += pi
                self.nr2 += pi
                self.nv2 += pi
                self.nr3 += pi
                self.nv3 += pi
                self.nr4 += pi
                self.nv4 += pi
            elif (ssm_node not in ancestors) and mr_cnv:
                self.nr1 += pi * (mr_cnv[1] + mr_cnv[2])
                self.nr2 += pi * (mr_cnv[1] + mr_cnv[2])
                self.nr3 += pi * (mr_cnv[1] + mr_cnv[2])
                self.nr4 += pi * (mr_cnv[1] + mr_cnv[2])
            elif ssm_node in ancestors and mr_cnv:
                self.nr3 += pi * max(0, (mr_cnv[1]+mr_cnv[2] - 1))
                self.nv3 += pi * min(1, mr_cnv[1]+mr_cnv[2])
                self.nr4 += pi * max(0, (mr_cnv[1] + mr_cnv[2] - 1))
                self.nv4 += pi * min(1, mr_cnv[1]+mr_cnv[2])

                if ssm_node in mr_cnv[0].node.get_ancestors():
                    self.nr1 = self.nr1 + pi * mr_cnv[1]
                    self.nv1 = self.nv1 + pi * mr_cnv[2]
                    self.nr2 = self.nr2 + pi * mr_cnv[2]
                    self.nv2 = self.nv2 + pi * mr_cnv[1]
                else:
                    self.nr1 = self.nr1 + pi * max(0, (mr_cnv[1]+mr_cnv[2] - 1))
                    self.nv1 = self.nv1 + pi * min(1, mr_cnv[1]+mr_cnv[2])
                    self.nr2 = self.nr2 + pi * max(0,
                                                   (mr_cnv[1] + mr_cnv[2] - 1))
                    self.nv2 = self.nv2 + pi * min(1, mr_cnv[1]+mr_cnv[2])
            else:
                print "PANIC"

        nodes = self.tssb.root['node'].tssb.get_nodes()
        self.nr1 = 0
        self.nv1 = 0
        self.nr2 = 0
        self.nv2 = 0
        self.nr3 = 0
        self.nv3 = 0
        self.nr4 = 0
        self.nv4 = 0
        for nd in nodes:
            descend(nd, new_state)
        if len(self.cnv) == 1 and self.node == self.cnv[0][0].node:
            out = [
                (self.nr1, self.nv1),
                (self.nr2, self.nv2),
                (self.nr3, self.nv3),
                (self.nr4, self.nv4)]
        else:
            out = [(self.nr1, self.nv1), (self.nr2, self.nv2)]
        return out

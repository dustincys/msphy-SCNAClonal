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
from stripe import Stripe, StripePool


class SegmentPool:

    self.baseline = -1 // 此处baseline是对数值
    self.segments = []

    def __init__(self, idx = 0, maxCopyNumber = 6, coverage = 30):
        self.idx = idx

        self.maxCopyNumber = maxCopyNumber
        self.coverage = coverage

        self._alleleConfig = get_cn_allele_config(max_copy_number)


    def load_seg_bed(self, bedName):
        """

        :bedName: TODO
        :returns: TODO

        """

        bedChroms, bedStarts, bedEnds, tReadNums, nReadNums, gcs =\
            BEDnParser(bedName)
        get_chrom_format(bedChroms)
        bedNum = len(bedChroms)

        for i in range(0, bedNum):
            chromIdx = chrom_name_to_idx(bedChroms[i])
            segName = get_segment_name(
                bedChroms[i], bedStarts[i], bedEnds[i])

            nReadNum = nReadNums[i]
            tReadNum = tReadNums[i]

            tempSeg = Segment()
            tempSeg.name = segName
            tempSeg.chromIdx = chromIdx
            tempSeg.chromName = bedChroms[i]
            tempSeg.start = bedStarts[i]
            tempSeg.end = bedEnds[i]
            tempSeg.nReadNum = nReadNum
            tempSeg.tReadNum = tReadNum

            # if 0 == nReadNum:
                # tempSeg.log2_ratio = -float('Inf')
            # else:
                # tempSeg.log2_ratio = np.log2(1.0 *
                                               # tReadNum/nReadNum)

            # tempSeg.log_ratio = np.log(1.0 * (tReadNum + 1.0) /
                                         # (nReadNum + 1.0))
            tempSeg.gc = gcs[i]

            self.segments.append(tempSeg)

    def load_seg_bam(self, normal_bam, tumor_bam, bed_file_name):
        chromIdxL = constants.CHROM_IDX_LIST
        chromStart = constants.CHROM_START

        samSQ = normal_bam.header['SQ']
        samChromFormat = get_chrom_format(map(lambda x: x['SN'], samSQ))
        chromLens, chromIdxs = get_chrom_lens_idxs(chromIdxL, samSQ)

        bedChroms, bedStarts, bedEnds, gcs = BEDParser(bed_file_name)
        get_chrom_format(bedChroms)
        bedNum = len(bedChroms)

        for i in range(0, bedNum):
            chromIdx = chrom_name_to_idx(bedChroms[i])
            chromName = chrom_idx_to_name(chromIdx, samChromFormat)
            segName = get_segment_name(chromName, bedStarts[i], bedEnds[i])

            if chromIdx not in chromIdxL:
                print 'Chromsome {0} not found, segment {1} excluded...'.format(
                    bedChroms[i], segName)
                sys.stdout.flush()
                continue

            chromListIdx = chromIdxs.index(chromIdx)

            if bedStarts[i] < chromStart or bedEnds[
                    i] > chromLens[chromListIdx]:
                print 'Out of range chromsome {0}, segment {1} excluded...'.\
                    format(bedChroms[i], segName)
                sys.stdout.flush()
                continue

            nReadNum = normal_bam.count(
                chromName, bedStarts[i], bedEnds[i])
            tReadNum = tumor_bam.count(
                chromName, bedStarts[i], bedEnds[i])

            tempSeg = Segment()
            tempSeg.name = segName
            tempSeg.chromIdx = chromIdx
            tempSeg.chromName = chromName
            tempSeg.start = bedStarts[i]
            tempSeg.end = bedEnds[i]
            tempSeg.nReadNum = nReadNum
            tempSeg.tReadNum = tReadNum
            # tempSeg.log2_ratio = np.log2(1.0 *
                                           # tReadNum/normal_reads_num)

            # tempSeg.log_ratio = np.log(1.0 * (tReadNum + 1.0) /
                                         # (normal_reads_num + 1.0))
            self.segments.append(tempSeg)

    def get_baseline(self, maxCopyNumber, subcloneNum, isPreprocess=False):
        self._get_LOH_frac()
        self._get_LOH_status()
        self._get_APM_frac()
        self._get_APM_status()
        self._calc_baseline()

        self._get_baseline_stripe()

        # 此处转换成了stripe之后再返回
        return self.get_seg_by_tag()

    def get_seg_by_tag(self, tag="BASELINE"):
        """
        return seg list
        """
        return filter(lambda item: item.tag==tag, self.segments)


    def _get_baseline_stripe(self):
        # 此处应该获取基线的条带，使用StripePool中的功能获取条带
        # def __init__(self, segPool, baseline=0.0, yDown, yUp, stripeNum, noiseStripeNum=2):

        # 获取yDown 和yUp 等相应的参数值，此处应该使用手动输入
        yDown = constants.YDOWNL[self.idx]
        yUp = constants.YUPL[self.idx]
        # 此处应该近似为最大拷贝数与亚克隆数的乘积，作为手工输入也可以
        stripeNum = constants.STRIPENUML[self.idx]
        noiseStripeNum = constants.NOISESTRIPENUML[self.idx]

        tempSP = StripePool(self, self.baseline, yDown, yUp, stripeNum, noiseStripeNum)
        tempSP.get()

        # 从条带中找到与基线最接近的条带作为基线条带

        self.nReadNum = -1.0
        self.tReadNum = -1.0

        tempDis = float("Inf")
        targetSP = None
        for sp in tempSP.stripes:
            spr = np.abs(self.baseline - np.log((sp.tReadNum + 1.0) /
                                                (sp.tReadNum + 1.0)))
            if tempDis > spr:
                tempDis = spr
                targetSP = sp

        # 将该条带之内的所有片段标注标签
        assert targetSP is not None

        for idx in targetSP.segsIdxL:
            self.segments[idx].tag = "BASELINE"


    def _calc_baseline(self, maxCopyNumber, subcloneNum, isPreprocess=False):
        """
        compute the Lambda S, through hierarchy clustering
        """
        if not isPreprocess:
            print >> sys.stdout, "compute_Lambda_S function called from model"
            return

        thresh = constants.HC_THRESH

        rdRatioLog = []
        for j in range(0, len(self.segments)):
            if self.segments[j].APMStatus == 'TRUE' and\
                    self.segments[j].LOHStatus == 'FALSE':
                ratio = self.segments[j].tReadNum*1.0/\
                    self.segments[j].normal_reads_num
                rdRatioLog.append(np.log(ratio))

        rdRatioLog = np.array(rdRatioLog)
        if rdRatioLog.shape[0] == 0:
            print >> sys.stderr, 'Error: no APM-LOH position found, existing...'
            print >> sys.stderr, 'Either the baselineThredAPM is too large, or\
            the constants APM_N_MIN is too large; Or, the baseline_thred_LOH is\
            too small'
            sys.exit(1)

        rdRatioLog = rdRatioLog.reshape(rdRatioLog.shape[0], 1)
        y = np.ones(rdRatioLog.shape)
        rdRatioLog = np.hstack((rdRatioLog, y))
        clusters = hcluster.fclusterdata(rdRatioLog,
                                         thresh,
                                         criterion="distance")
        mccs = Counter(clusters).most_common(maxCopyNumber * subcloneNum)

        rdrMinLog = float('Inf')
        clusterMin = -1
        for i in range(0, len(mccs)):
            clusterTemp = mccs[i][0]
            print >> sys.stdout, "cluster temp : {}".format(clusterTemp)
            tempRdrLog = mean(rdRatioLog[clusters == clusterTemp])
            print >> sys.stdout, "tempRdrLog"
            print >> sys.stdout, "log: {}".format(tempRdrLog)
            if rdrMinLog > tempRdrLog:
                rdrMinLog = tempRdrLog
                clusterMin = clusterTemp

        print >> sys.stdout, mccs
        print >> sys.stdout, "log baseline: {}".format(rdrMinLog)
        sys.stdout.flush()

        clusterFlag = (clusters == clusterMin)
        baselineNum = 0
        rdrIdx = 0
        for j in range(0, len(self.segments)):
            if self.segments[j].APMStatus == 'TRUE' and\
                    self.segments[j].LOHStatus == 'FALSE':
                if clusterFlag[rdrIdx]:
                    self.segments[j].baselineLabel = 'TRUE'
                    baselineNum = baselineNum + 1
                else:
                    self.segments[j].baselineLabel = 'FALSE'
                rdrIdx = rdrIdx + 1
            else:
                self.segments[j].baselineLabel = 'FALSE'

        print >> sys.stdout, "baselineNum: {}".format(baselineNum)

        if baselineNum == 0:
            print >> sys.stderr, 'Error: No diploid segments found, existing...'
            sys.exit(1)

        self.baseline = rdrMinLog

    def _get_LOH_frac(self):
        for j in range(0, len(self.segments)):
            self.segments[j].LOHFrac = get_LOH_frac(
                self.segments[j].pairedCounts)

    def _get_APM_frac(self):
        for j in range(0, len(self.segments)):
            self.segments[j].APMFrac = get_APM_frac_MAXMIN(
                self.segments[j].pairedCounts)

    def _get_LOH_status(self, baseThred, isPreprocess=False):
        if isPreprocess:
            LOHNum = 0
            FLOHNum = 0
            for j in range(0, len(self.segments)):
                self.segments[j].LOHStatus = get_LOH_status(
                    self.segments[j].LOHFrac, baseThred)
                if self.segments[j].LOHStatus == "TRUE":
                    LOHNum = LOHNum + 1
                elif self.segments[j].LOHStatus == "FALSE":
                    FLOHNum = FLOHNum + 1

            print >> sys.stdout, "LOHNum/segNum = {0}/{1}".format(LOHNum, len(self.segments))
            print >> sys.stdout, "FLOHNum/segNum = {0}/{1}".format(FLOHNum, len(self.segments))
        else:
            print >> sys.stdout, "get_LOH_status function called from model."

    def _get_APM_status(self, baselineThredAPM):
        APMNum = 0
        for j in range(0, len(self.segments)):
            self.segments[j].APMStatus = get_APM_status(
                self.segments[j].APMFrac, baselineThredAPM)
            if self.segments[j].APMStatus == "TRUE":
                APMNum = APMNum + 1

        print "APMNum/segNum = {0}/{1}".format(APMNum, len(self.segments))

    # def _log_likelihood(self, id, phi, update_tree=True, new_state=0):

        # if update_tree:
            # ##################################################
            # # some useful info about the tree,
            # # used by CNV related computations,
            # u.set_node_height(self.tssb)
            # u.set_path_from_root_to_node(self.tssb)
            # u.map_datum_to_node(self.tssb)
            # ##################################################

        # seg = self.segments[id]
        # ll, cn, pi = self._getSegResData(seg, phi)

        # return ll

    # def _getSegResData(self, seg, phi):
        # copyNumbers = None
        # if seg.baselineLabel == "True":
            # copyNumbers = [2]
        # elif get_loga(seg) > self.baseline:
            # copyNumbers = range(2, self.max_copy_number + 1)
        # else:
            # copyNumbers = range(0, 2 + 1)

        # llPiS = [self._getLLSeg(seg, copyNumber, phi) for copyNumber in
                   # copyNumbers]
        # (ll, pi) = max(llPiS, key=lambda x: x[0])
        # cn = llPiS.index((ll, pi))
        # return ll, cn, pi

    # def _getLLSeg(self, seg, copyNumber, phi):
        # llSeg = 0
        # llRd = self._getRD(seg, copyNumber, phi)
        # alleleType = self._alleleConfig[copyNumber]
        # self._augBAF(seg, copyNumber)
        # if 0 == seg.pairedCounts.shape[0]:
            # llBAFs = 0
            # pi = "*"
        # else:
            # llBAFs, pi = self._getBAF(seg, copyNumber, alleleType, phi)
        # llSeg = llBAFs + llRd
        # return llSeg, pi

    # def _augBAF(self, seg, copyNumber):
        # if copyNumber > 2:
            # threshold = constants.BAF_THRESHOLD * self._coverage
            # dTj = np.sum(seg.BAF[:, 2:4], axis=1)
            # idxRm = tuple(np.where(dTj < threshold)[0])
            # seg.BAF = np.delete(seg.BAF, idxRm, axis=0)
        # else:
            # pass

    # def _getRD(self, seg, copyNumber, phi):
        # cN = constants.COPY_NUMBER_NORMAL
        # bar_c = phi * copyNumber + (1.0 - phi) * cN
        # print "____>>> _getRD: bar_c, cN, self._baseline, seg.normal_reads_num____"
        # print bar_c, cN, self._baseline, seg.normal_reads_num
        # print "_________end _getRD:bar_c, cN, self._baseline, seg.normal_reads_num______________"

        # lambda_possion = (
            # bar_c / cN) * self._baseline * (seg.normal_reads_num + 1) #not minus 1 ? better
        # if lambda_possion < 0:
            # lambda_possion = 0
        # print "____>>> _getRD: seg.tReadNum, lambda_possion____"
        # print seg.tReadNum, lambda_possion
        # print "_________end _getRD:seg.tReadNum, lambda_possion______________"

        # ll_RD = log_poisson_pdf(seg.tReadNum, lambda_possion)
        # return ll_RD

    # def _getBAF(self, seg, copyNumber, alleleType, phi):
        # cN = constants.COPY_NUMBER_NORMAL
        # mu_N = constants.MU_N
        # # keys, ppmm values 0.5
        # mu_G = np.array(alleleType.values())

        # print "____>>> _getBAF: mu_N, mu_G, cN, copyNumber, phi____"
        # print mu_N, mu_G, cN, copyNumber, phi
        # print "_________end _getBAF:mu_N, mu_G, cN, copyNumber, phi______________"

        # mu_E = get_mu_E_joint(mu_N, mu_G, cN, copyNumber, phi)

        # if seg.pairedCounts.shape[0] > 1:
            # b_T_j = np.min(seg.pairedCounts[:, 2:4], axis=1)
            # dTj = np.sum(seg.pairedCounts[:, 2:4], axis=1)
            # baf = b_T_j * 1.0 / dTj
            # outlier = mad_based_outlier(baf)
            # BAF = np.delete(seg.pairedCounts, list(outlier.astype(int)), axis=0)
            # b_T_j = np.min(BAF[:, 2:4], axis=1)
            # dTj = np.sum(BAF[:, 2:4], axis=1)

        # else:
            # b_T_j = np.min(seg.pairedCounts[:, 2:4], axis=1)
            # dTj = np.sum(seg.pairedCounts[:, 2:4], axis=1)
            # pass

        # # add prior or not?
        # ll = log_binomial_likelihood(b_T_j, dTj, mu_E)
        # ll_bafs = ll.sum(axis=0)
        # idx_max = ll_bafs.argmax(axis=0)
        # llBAFs = ll_bafs[idx_max]
        # pi = alleleType[alleleType.keys()[idx_max]]
        # return llBAFs, pi

    # # computes the binomial parameter
    # def compute_n_genomes(self, tp, new_state=0):
        # def descend(nd, new_state):
            # # this is needed for Metropolis-Hastings likelihood computations
            # pi = nd.pi1[tp] if new_state else nd.pi[tp]
            # ssm_node = self.node.path[-1]
            # mr_cnv = self.find_most_recent_cnv(nd)
            # ancestors = nd.get_ancestors()
            # if (ssm_node not in ancestors) and (not mr_cnv):
                # self.nr1 += pi * 2
                # self.nr2 += pi * 2
                # self.nr3 += pi * 2
                # self.nr4 += pi * 2
            # elif ssm_node in ancestors and (not mr_cnv):
                # self.nr1 += pi
                # self.nv1 += pi
                # self.nr2 += pi
                # self.nv2 += pi
                # self.nr3 += pi
                # self.nv3 += pi
                # self.nr4 += pi
                # self.nv4 += pi
            # elif (ssm_node not in ancestors) and mr_cnv:
                # self.nr1 += pi * (mr_cnv[1] + mr_cnv[2])
                # self.nr2 += pi * (mr_cnv[1] + mr_cnv[2])
                # self.nr3 += pi * (mr_cnv[1] + mr_cnv[2])
                # self.nr4 += pi * (mr_cnv[1] + mr_cnv[2])
            # elif ssm_node in ancestors and mr_cnv:
                # self.nr3 += pi * max(0, (mr_cnv[1]+mr_cnv[2] - 1))
                # self.nv3 += pi * min(1, mr_cnv[1]+mr_cnv[2])
                # self.nr4 += pi * max(0, (mr_cnv[1] + mr_cnv[2] - 1))
                # self.nv4 += pi * min(1, mr_cnv[1]+mr_cnv[2])

                # if ssm_node in mr_cnv[0].node.get_ancestors():
                    # self.nr1 = self.nr1 + pi * mr_cnv[1]
                    # self.nv1 = self.nv1 + pi * mr_cnv[2]
                    # self.nr2 = self.nr2 + pi * mr_cnv[2]
                    # self.nv2 = self.nv2 + pi * mr_cnv[1]
                # else:
                    # self.nr1 = self.nr1 + pi * max(0, (mr_cnv[1]+mr_cnv[2] - 1))
                    # self.nv1 = self.nv1 + pi * min(1, mr_cnv[1]+mr_cnv[2])
                    # self.nr2 = self.nr2 + pi * max(0,
                                                   # (mr_cnv[1] + mr_cnv[2] - 1))
                    # self.nv2 = self.nv2 + pi * min(1, mr_cnv[1]+mr_cnv[2])
            # else:
                # print "PANIC"

        # nodes = self.tssb.root['node'].tssb.get_nodes()
        # self.nr1 = 0
        # self.nv1 = 0
        # self.nr2 = 0
        # self.nv2 = 0
        # self.nr3 = 0
        # self.nv3 = 0
        # self.nr4 = 0
        # self.nv4 = 0
        # for nd in nodes:
            # descend(nd, new_state)
        # if len(self.cnv) == 1 and self.node == self.cnv[0][0].node:
            # out = [
                # (self.nr1, self.nv1),
                # (self.nr2, self.nv2),
                # (self.nr3, self.nv3),
                # (self.nr4, self.nv4)]
        # else:
            # out = [(self.nr1, self.nv1), (self.nr2, self.nv2)]
        # return out

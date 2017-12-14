'''
# =============================================================================
#      FileName: BamToDataConverter.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2016-10-08 10:54:48
#       History:
# =============================================================================
'''

import sys
import pickle as pkl
from multiprocessing import Pool
import numpy as np
import pysam

from mcmc import MCMCLM

from phySCNAClonal.preprocess.data import Data
from phySCNAClonal.preprocess.stripe import Stripe, DataStripes
from phySCNAClonal.preprocess.iofun import PairedCountsIterator, PairedPileupIterator
from phySCNAClonal.preprocess.utils import get_BAF_counts, normal_heterozygous_filter

from plotGC import GCStripePlot


class phySCNAClonal_Converter:

    def __init__(self, normal_bam_filename, tumor_bam_filename,
                 reference_genome_filename, input_filename_base, segments_bed,
                 BICseq_bed_fileName_corrected, pkl_path="",
                 max_copynumber=6, subclone_num=1, baseline_thred_LOH=0.3,
                 baseline_thred_APM=0.01, min_depth=20, min_bqual=10,
                 min_mqual=10,  process_num=1):
        self.normal_bam_filename = normal_bam_filename
        self.tumor_bam_filename = tumor_bam_filename
        self.reference_genome_filename = reference_genome_filename
        self.input_filename_base = input_filename_base
        self.segments_bed = segments_bed
        self.BICseq_bed_fileName_corrected = BICseq_bed_fileName_corrected
        self.pkl_path = pkl_path

        self.max_copynumber = max_copynumber
        self.subclone_num = subclone_num
        self.baseline_thred_LOH = baseline_thred_LOH
        self.baseline_thred_APM = baseline_thred_APM

        self.min_depth = min_depth
        self.min_bqual = min_bqual
        self.min_mqual = min_mqual
        print "process_num = {}".format(process_num)
        self.process_num = process_num

        self.data = Data()
        self.stripes = None

    def convert(self, method, pkl_flag=False):
        if pkl_flag and self.pkl_path != "":
            print "load pkl from"
            print self.pkl_path
            infile = open(self.pkl_path, 'rb')
            self.data = pkl.load(infile)
            infile.close()
        else:
            self._load_segments_bed()
            print "phySCNAClonal converter converting"

            if "auto" == method:
                self._MCMC_gccorrection()
            elif "visual" == method:
                self._visual_gccorrection()
                sys.stdout.flush()
            self._get_counts()

        self.visualize()
        self._baseline_selection()

        data_file_name = self.input_filename_base + '.phySCNAClonal.input.pkl'
        outfile = open(data_file_name, 'wb')
        pkl.dump(self.data, outfile, protocol=2)
        outfile.close()


        self.dataStripes = DataStripes(self.data)
        self.dataStripes.get()

        stripes_file_name = self.input_filename_base + '.phySCNAClonal.stripes.input.pkl'
        outfile = open(stripes_file_name, 'wb')
        pkl.dump(self.dataStripes, outfile, protocol=2)
        outfile.close()

        stripes_file_name = self.input_filename_base + '.phySCNAClonal.stripes.input.txt'
        self.dataStripes.output_txt(stripes_file_name)

    def _MCMC_gccorrection(self):
        """
        The interception is irrelevant for correction, set as median
        MCMCLM only returns the m and c, then correct the data here
        """
        mcmclm = MCMCLM(self.data, 0, self.subclone_num, self.max_copynumber)
        m, c = mcmclm.run()
        print "MCMC slope = {}".format(m)
        self._correct(m, c)

    def _correct(self, slope, intercept):

        x = np.array(map(lambda seg: seg.gc, self.data.segments))
        y = np.array(map(lambda seg: np.log(seg.tumor_reads_num + 1) -
                         np.log(seg.normal_reads_num + 1),
                         self.data.segments))

        K = np.percentile(y, 50)
        A = slope * x + intercept
        y_corrected = y - A + K

        for i in range(len(y_corrected)):
            self.data.segments[i].tumor_reads_num = np.exp(
                y_corrected[i] +
                np.log(self.data.segments[i].normal_reads_num + 1)
            ) - 1
            self.data.segments[i].log_ratio = np.log(
                (y_corrected[i] + 1.0) /
                (self.data.segments[i].normal_reads_num + 1.0)
            )

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(slope, intercept)

    def visualize(self):
        gsp = GCStripePlot(self.data.segments, len(self.data.segments))
        print "total number: {}".format(self.data.seg_num)
        gsp.plot()
        x, y, m, c = gsp.output()
        print "x, y, m, c"
        print x, y, m, c

    def _visual_gccorrection(self):
        gsp = GCStripePlot(self.data.segments, len(self.data.segments))
        print "total number: {}".format(self.data.seg_num)
        gsp.plot()
        x, y, m, c = gsp.output()
        print "x, y, m, c"
        print x, y, m, c
        self._correct(m, c)

    def _baseline_selection(self):
        print "begin baseline selection.."
        self._get_LOH_frac()
        self._get_LOH_status()
        self._get_APM_frac()
        self._get_APM_status()
        self._compute_Lambda_S()

    def _get_APM_status(self):
        self.data.get_APM_status(self.baseline_thred_APM)

    def _get_LOH_status(self):
        self.data.get_LOH_status(self.baseline_thred_LOH,
                                 flag_runpreprocess=True)

    def _compute_Lambda_S(self):
        print "begin compute lambda s .."
        self.data.compute_Lambda_S_LOH(self.max_copynumber, self.subclone_num,
                                       flag_runpreprocess=True)

    def _load_segmentsn(self):
        """
        :returns: TODO

        """
        normal_bam = pysam.Samfile(self.normal_bam_filename, 'rb')
        tumor_bam = pysam.Samfile(self.tumor_bam_filename, 'rb')

        print 'Loading normalized segments by {0}...'.format(self.segments_bed)
        sys.stdout.flush()
        self.data.load_segmentsn(normal_bam, tumor_bam, self.segments_bed)

        normal_bam.close()
        tumor_bam.close()

    def _load_segments(self):
        normal_bam = pysam.Samfile(self.normal_bam_filename, 'rb')
        tumor_bam = pysam.Samfile(self.tumor_bam_filename, 'rb')

        print 'Loading segments with gc by {0}...'.format(self.segments_bed)
        sys.stdout.flush()
        self.data.load_segments(normal_bam, tumor_bam, self.segments_bed)

        normal_bam.close()
        tumor_bam.close()

    def _load_segments_bed(self):
        print 'Loading segments with gc by {0}...'.format(self.segments_bed)
        sys.stdout.flush()
        self.data.load_segments_bed(self.segments_bed)

    def _get_counts(self):
        seg_num = self.data.seg_num
        process_num = self.process_num
        print "process_num = {}".format(process_num)

        if process_num > seg_num:
            process_num = seg_num

        pool = Pool(processes=process_num)

        args_list = []

        for j in range(0, seg_num):
            seg_name = self.data.segments[j].name
            chrom_name = self.data.segments[j].chrom_name
            chrom_idx = self.data.segments[j].chrom_idx
            start = self.data.segments[j].start
            end = self.data.segments[j].end

            args_tuple = (
                seg_name,
                chrom_name,
                chrom_idx,
                start,
                end,
                self.normal_bam_filename,
                self.tumor_bam_filename,
                self.reference_genome_filename,
                self.min_depth,
                self.min_bqual,
                self.min_mqual)

            args_list.append(args_tuple)

        counts_tuple_list = pool.map(process_by_segment, args_list)

        for j in range(0, seg_num):
            paired_counts_j, BAF_counts_j = counts_tuple_list[j]

            self.data.segments[j].paired_counts = paired_counts_j
            self.data.segments[j].BAF_counts = BAF_counts_j

    def _get_LOH_frac(self):
        self.data.get_LOH_frac()

    def _get_APM_frac(self):
        self.data.get_APM_frac()
# ===============================================================================
#  Function
# ===============================================================================


def process_by_segment(args_tuple):
    seg_name, chrom_name, chrom_idx, start, end, normal_bam_filename,\
        tumor_bam_filename, reference_genome_filename, min_depth, min_bqual,\
        min_mqual = args_tuple

    print 'Preprocessing segment {0}...'.format(seg_name)
    sys.stdout.flush()

    normal_bam = pysam.Samfile(normal_bam_filename, 'rb')
    tumor_bam = pysam.Samfile(tumor_bam_filename, 'rb')
    ref_genome_fasta = pysam.Fastafile(reference_genome_filename)

    normal_pileup_iter = normal_bam.pileup(chrom_name, start, end)
    tumor_pileup_iter = tumor_bam.pileup(chrom_name, start, end)

    paired_pileup_iter = PairedPileupIterator(
        normal_pileup_iter, tumor_pileup_iter, start, end)
    paired_counts_iter = PairedCountsIterator(
        paired_pileup_iter,
        ref_genome_fasta,
        chrom_name,
        chrom_idx,
        min_depth,
        min_bqual,
        min_mqual)

    paired_counts_j, BAF_counts_j = iterator_to_counts(paired_counts_iter)
    counts_tuple_j = (paired_counts_j, BAF_counts_j)

    normal_bam.close()
    tumor_bam.close()
    ref_genome_fasta.close()

    return counts_tuple_j


def iterator_to_counts(paired_counts_iter):
    buff = 100000

    paired_counts_j = np.array([[], [], [], [], [], []], dtype=int).transpose()
    BAF_counts_j = np.zeros((100, 100))
    buff_counts = []
    i = 0

    for counts in paired_counts_iter:
        buff_counts.append(counts)
        i = i + 1

        if i < buff:
            continue

        buff_counts = np.array(buff_counts)

        if buff_counts.shape[0] != 0:
            BAF_counts_buff = get_BAF_counts(buff_counts)
            BAF_counts_j += BAF_counts_buff

        buff_counts_filtered = normal_heterozygous_filter(buff_counts)

        if buff_counts_filtered.shape[0] != 0:
            paired_counts_j = np.vstack((paired_counts_j, buff_counts_filtered))

        buff_counts = []
        i = 0

    buff_counts = np.array(buff_counts)

    if buff_counts.shape[0] != 0:
        BAF_counts_buff = get_BAF_counts(buff_counts)
        BAF_counts_j += BAF_counts_buff

    buff_counts_filtered = normal_heterozygous_filter(buff_counts)

    if buff_counts_filtered.shape[0] != 0:
        paired_counts_j = np.vstack((paired_counts_j, buff_counts_filtered))

    return (paired_counts_j, BAF_counts_j)

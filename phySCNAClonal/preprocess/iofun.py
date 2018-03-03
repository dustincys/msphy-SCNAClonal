'''
# =============================================================================
#      FileName: iofun.py
#          Desc: functions for get allele info from bam
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-03-03 08:54:49
#       History: Andrew Roth, Yi Li
# =============================================================================
'''
import sys
from collections import Counter

import numpy as np
import pysam

ascii_offset = 33

#JointSNVMix
class PairedCountsIterator:
    def __init__(self, pairedPileupIter, refGenomeFasta, chromName, chromIdx,
                 minDepth=20, minBqual=10, minMqual=10):
        self.pairedPileupIter = pairedPileupIter
        self.refGenomeFasta = refGenomeFasta
        self.chromName = chromName
        self.chromIdx = chromIdx
        self.minDepth = minDepth
        self.minBqual = minBqual
        self.minMqual = minMqual

    def __iter__(self):
        return self

    def next(self):
        normalColumn, tumorColumn = self.pairedPileupIter.next()

        while True:
            if normalColumn.n < self.minDepth or tumorColumn.n < self.minDepth:
                normalColumn, tumorColumn = self.pairedPileupIter.next()
                continue

            pos = normalColumn.pos
            chromIdx = self.chromIdx
            refBase = self.refGenomeFasta.fetch(self.chromName, pos, pos + 1).upper()

            if refBase == '':
                print 'Error: %s does not match the reference of the bam files' \
                % self.refGenomeFasta.filename
                sys.exit(-1)

            pairedCounts = self._get_paired_counts(normalColumn, tumorColumn, chromIdx, pos, refBase)

            if pairedCounts == None:
                normalColumn, tumorColumn = self.pairedPileupIter.next()
                continue
            else:
                return pairedCounts

    def _get_paired_counts(self, normalColumn, tumorColumn, chromIdx, pos, refBase):
        normalBases = self._parse_pileup_column(normalColumn)
        tumorBases = self._parse_pileup_column(tumorColumn)

        normalNonRefBase, normalCounts = self._get_counts(refBase, normalBases)
        tumorNonRefBase, tumorCounts = self._get_counts(refBase, tumorBases)

        # Check again for lines below read depth. The first check above speeds things up, though redundant.
        normalDepth = normalCounts[0] + normalCounts[1]
        tumorDepth = tumorCounts[0] + tumorCounts[1]

        if normalDepth < self.minDepth or tumorDepth < self.minDepth:
            return None

        # Shift index to one based position.
        oneBasedPos = pos + 1

        pairedCounts = []
        pairedCounts.extend(normalCounts)
        pairedCounts.extend(tumorCounts)
        pairedCounts.append(chromIdx)
        pairedCounts.append(pos)

        return pairedCounts

    def _parse_pileup_column(self, pileupColumn):
        bases = []

        for read in pileupColumn.pileups:
            if read.is_del:
                continue

            if hasattr(read, 'query_position') == True:
                qpos = read.query_position
                bqual = read.alignment.query_qualities[qpos]
            elif hasattr(read, 'qpos') == True:
                qpos = read.qpos
                bqual = ord(read.alignment.qual[qpos]) - ascii_offset
            else:
                raise Exception("Error in pysam qpos/query_position.")

            mqual = read.alignment.mapq

            if mqual < self.minMqual:
                continue

            if bqual < self.minBqual:
                continue

            base = read.alignment.seq[qpos].upper()
            bases.append(base)

        return bases

    def _get_counts(self, refBase, bases, nonRefBase=None):
        counter = Counter(bases)

        nonRefBase, counts = self._parse_counts(refBase, counter, nonRefBase)

        return nonRefBase, counts

    def _parse_counts(self, refBase, counter, nonRefBase=None):
        refCounts = counter[refBase]

        del counter[refBase]
        del counter['N']

        # Check if there is any non-ref bases.
        if nonRefBase is not None:
            nonRefCounts = counter[nonRefBase]
        else:
            if len(counter) > 0:
                nonRefBase, nonRefCounts = counter.most_common(1)[0]
            else:
                nonRefBase = 'N'
                nonRefCounts = 0

        counts = (refCounts, nonRefCounts)

        return nonRefBase, counts

#JointSNVMix
class PairedPileupIterator:
    def __init__(self, normalIter, tumorIter, segmentStart, segmentEnd):
        self.normalIter = normalIter
        self.tumorIter = tumorIter
        self.segmentStart = segmentStart
        self.segmentEnd = segmentEnd

    def __iter__(self):
        return self

    def next(self):
        normalColumn = self.normalIter.next()
        tumorColumn = self.tumorIter.next()

        while True:
            normalPos = normalColumn.pos
            tumorPos = tumorColumn.pos

            if normalPos == tumorPos:
                if normalPos >= self.segmentStart and normalPos <= self.segmentEnd:
                    return normalColumn, tumorColumn
                else:
                    normalColumn = self.normalIter.next()
                    tumorColumn = self.tumorIter.next()
            elif normalPos < tumorPos:
                normalColumn = self.normalIter.next()
            elif normalPos > tumorPos:
                tumorColumn = self.tumorIter.next()
            else:
                raise Exception("Error in paired pileup iterator.")

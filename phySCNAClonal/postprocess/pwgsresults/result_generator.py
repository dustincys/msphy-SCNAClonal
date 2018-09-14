import numpy as np
from collections import defaultdict

import sys
import os
sys.path.append(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import phySCNAClonal.model.util2
import json


class ResultGenerator(object):
    def generate(self, treeFile, SCNALFile):
        """
        Load SCNA data from SCNALFile, write population and data parameters
        into SCNA data file
        """
        reader = util2.TreeReader(treeFile)
        try:
            # to record params, manually
            params = json.loads(reader.read_extra_file('params.json'))
        except KeyError:
            # File not present in archive, likely because it originates from an
            # older run.
            params = {}
        reader.close()

        SCNAPool, isStripe = self.__load_SCNAL_from_pkl(SCNALFile)

        summaries = {}
        allMutAss = {}
        allMutDataParam = {}

        for idx, llh, dp, pops, mutAss, structure, mutPops in\
                self._summarize_all_pops(treeFile):
            summaries[idx] = {
                'llh': llh,
                'structure': structure,
                'populations': pops
            }
            allMutDataParam[idx] = dp
            allMutAss[idx] = mutAss

            currentSCNAPool = copy.deepcopy(SCNAPool)
            self._update_SCNAPool(mutPops, dp, currentSCNAPool, isStripe)
            summaries[idx]['SCNAPool'] = currentSCNAPool

        return summaries, allMutAss, params, isStripe

    def _update_SCNAPool(self, mutPops, dp, SCNAPool, isStripe):

        if isStripe:
            for stripe, sdp, idx in zip(
                SCNAPool.stripes, dp, range(SCNAPool.stripes)):
                assert stripe.sid == dp[1]
                assert idx == int(dp[0])

                pop = mutPops[dp[1]]

                for segIdx in stripe.segsIdxL:
                    SCNAPool.segPool.segments[segIdx].copyNumber = in(dp[2])
                    SCNAPool.segPool.segments[segIdx].genotype = str(dp[2])
                    SCNAPool.segPool.segments[segIdx].phi = pop
                    pass
                pass

        else:
            for segment, sdp, idx in zip(
                SCNAPool.segments, dp, range(SCNAPool.segments)):

                segmentId = "{0}_{1}_{2}".format(
                    segment.chromName, str(segment.start), str(segment.end))
                assert segmentId == dp[1]
                assert idx == int(dp[0])

                pop = mutPops[dp[1]]

                segments.copyNumber = in(dp[2])
                segments.genotype = str(dp[2])
                segments.phi = pop
                pass

    def _summarize_all_pops(self, treeFile, isStripe):
        # 此处根据SCNAPool的类型进行生成DataParameter
        reader = util2.TreeReader(treeFile)
        for idx, llh, tree, dp in reader.load_trees_and_metadata(
                removeEmptyVertices=True):
            yield (idx, llh, dp) + self._summarize_pops(tree)
        reader.close()

    def _summarize_pops(self, tree):
        pops = {}
        mutPops = {}
        structure = defaultdict(list)
        # Note that there will be an entry in mutAssignments for a given
        # subclone only if it has at least one SSM or CNV. This assumption
        # holds true for all PhyloWGS trees, so one can ascertain the number of
        # cancerous populations via len(mutAssignments).

        mutAssignments = defaultdict(lambda: [])
        idx = [0]

        def _traverse_r(vertex, parent):
            mutations = vertex.get_data()
            # vertex.params represents phis (i.e., population freqs) associated
            # with each sample.
            cellPrev = vertex.param
            currentIdx = idx[0]

            SCNANum = 0
            for mut in mutations:
                # 此处最好要加入segIdx，直接进行后续的分析
                # 需要记录位置
                mutPops[mut.sid] = cellPrev
                mutAssignments[currentIdx].append(mut.sid)
                SCNANum += 1

            # Preorder traversal is consistent with printo_latex.py, meaning index
            # values should correspond to same vertices.
            pops[currentIdx] = {
                'cellular_prevalence': cellPrev,
                'num_SCNAs': SCNANum,
            }

            # Visit children in order of decreasing phi.
            children = sorted(
                vertex.children(),
                key=lambda v: v.param,
                reverse=True)
            for child in children:
                idx[0] += 1
                structure[currentIdx].append(idx[0])
                _traverse_r(child, currentIdx)

        _traverse_r(tree.root['node'], None)
        return (pops, mutAssignments, structure, mutPops)

    def __load_SCNAL_from_pkl(self, inputFilePath):
        inputFile = open(inputFilePath, 'rb')
        SCNAPoolL = pkl.load(inputFile)
        lastSCNAPool = SCNAPoolL[-1]
        if type(SCNAPool).__name__ == "StripePool":
            return lastSCNAPool, True
        elif type(SCNAPool).__name__ == "SegmentPool":
            return lastSCNAPool, False

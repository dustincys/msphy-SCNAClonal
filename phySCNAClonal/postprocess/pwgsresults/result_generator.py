#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: result_generator.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-09-16 13:48:50
#       History:
# =============================================================================
'''

import copy
import json
import os
import pickle as pkl
import sys
from collections import defaultdict

import numpy as np

from phySCNAClonal.model.util2 import TreeReader

sys.path.append(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))


class ResultGenerator(object):
    def generate(self, treeFile, SCNALFile):
        """
        Load SCNA data from SCNALFile, write population and data parameters
        into SCNA data file
        """
        treeNum = -1
        reader = TreeReader(treeFile)
        try:
            # to record params, manually
            params = json.loads(reader.read_extra_file('params.json'))
            treeNum = reader.num_trees()
        except KeyError:
            # File not present in archive, likely because it originates from an
            # older run.
            params = {}
        reader.close()

        SCNAPool, isStripe = self.__load_SCNAL_from_pkl(SCNALFile)

        summaries = {}
        allMutAss = {}
        allMutDataParam = {}
        partialDict = {}

        for idx, llh, dp, tree, pops, mutAss, structure, mutPops in\
                self._summarize_all_pops(treeFile, isStripe):

            self._generate_partialData(tree, partialDict)

            summaries[idx] = {
                'llh': llh,
                'structure': structure,
                'populations': pops,
                'tree': tree
            }
            allMutDataParam[idx] = dp
            allMutAss[idx] = mutAss

            currentSCNAPool = copy.deepcopy(SCNAPool)
            self._update_SCNAPool(mutPops, dp, currentSCNAPool, isStripe)
            summaries[idx]['SCNAPool'] = currentSCNAPool

        return treeNum, summaries, allMutAss, params, partialDict, isStripe


    def _generate_partialData(self, tree, pd):

        def descend(root):
            childrenIds = set()
            for node in root['children']:
                childrenIds = childrenIds.union(descend(node))
            if len(root['node'].data) > 0:
                for dataId in root['node'].data:
                    for childId in childrenIds:
                        self._update_partialDict(pd, dataId, childId)
            else:
                return childrenIds

            if len(root['node'].data) > 0:
                return root['node'].data

        descend(tree.root)


    def _update_partialDict(self, pd, fromNode, toNode):
        if fromNode in pd.keys():
            if toNode in pd[fromNode].keys():
                pd[fromNode][toNode] += 1.0
            else:
                pd[fromNode][toNode] = 1.0
        else:
            pd[fromNode] = {}
            pd[fromNode][toNode] = 1.0


    def _update_SCNAPool(self, mutPops, dp, SCNAPool, isStripe):

        blSegsL = None

        if isStripe:
            for stripe, sdp, idx in zip(
                SCNAPool.stripes, dp, range(len(SCNAPool.stripes))):
                assert stripe.name == sdp[1]
                assert idx == int(sdp[0])

                pop = mutPops[sdp[1]]

                for segIdx in stripe.segsIdxL:
                    targetSeg = SCNAPool.segPool.segments[segIdx]
                    targetSeg.copyNumber = int(sdp[2])
                    targetSeg.genotype = str(sdp[3])
                    targetSeg.phi = pop
                    pass
                pass

            blSegsL = filter(lambda item: item.tag == 'BASELINE',
                             SCNAPool.segPool.segments)
        else:
            for segment, sdp, idx in zip(
                SCNAPool.segments, dp, range(len(SCNAPool.segments))):

                segmentId = "{0}_{1}_{2}".format(
                    segment.chromName, str(segment.start), str(segment.end))
                assert segmentId == sdp[1]
                assert idx == int(sdp[0])

                pop = mutPops[sdp[1]]

                segment.copyNumber = int(sdp[2])
                segment.genotype = str(sdp[2])
                segment.phi = pop
                if segment.tag == 'BASELINE':
                    segment.copyNumber = 2
                    segment.genotype = 'PM'
                    segment.phi = 1.0
                    pass
                pass

            blSegsL = filter(lambda item: item.tag == 'BASELINE',
                             SCNAPool.segments)

        for blSeg in blSegsL:
            blSeg.copyNumber = 2
            blSeg.genotype = 'PM'
            blSeg.phi = 1.0
            pass

    def _summarize_all_pops(self, treeFile, isStripe):
        # 此处根据SCNAPool的类型进行生成DataParameter
        reader = TreeReader(treeFile)
        for idx, llh, tree, dp in reader.load_trees_and_metadata(
                removeEmptyVertices=True):
            yield (idx, llh, dp, tree) + self._summarize_pops(tree)
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
        pklFile = open(inputFilePath, 'rb')
        SCNAPool = pkl.load(pklFile)
        pklFile.close()
        if not hasattr(SCNAPool, 'segments'):
            return SCNAPool, True
        else:
            return SCNAPool, False

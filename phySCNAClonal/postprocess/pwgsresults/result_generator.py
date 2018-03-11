import numpy as np
from collections import defaultdict

import sys
import os
sys.path.append(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import phySCNAClonal.model.util2
import json


class ResultGenerator(object):
    def generate(self, treeFile, segPoolLFile, includeStripeNames):
        """
        segPoolL: load the segment information
        """
        reader = util2.TreeReader(treeFile)
        firstTree = next(reader.load_trees())
        # cnv_logical_physical_mapping = json.loads(reader.read_extra_file('cnv_logical_physical_mapping.json'))
        try:
            params = json.loads(reader.read_extra_file('params.json'))
        except KeyError:
            # File not present in archive, likely because it originates from an
            # older run.
            params = {}
        reader.close()

        segPool = self.__load_segPoolL_from_pkl(segPoolLFile)

        # 注意此处是否添加add stripe name
        mutList = self._list_mutations(firstTree, segPool, includeStripeNames)
        summaries = {}
        allMutAss = {}
        for idx, llh, pops, mutAss, structure in self._summarize_all_pops(
                treeFile):
            summaries[idx] = {
                'llh': llh,
                'structure': structure,
                'populations': pops,
            }
            allMutAss[idx] = mutAss

        return summaries, mutList, allMutAss, params

    def _summarize_all_pops(self, treeFile):
        reader = util2.TreeReader(treeFile)
        for idx, llh, tree in reader.load_trees_and_metadata(
                removeEmptyVertices=True):
            yield (idx, llh) + self._summarize_pops(tree)
        reader.close()

    def _summarize_pops(self, tree):
        pops = {}
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

            stripeNum = 0
            for mut in mutations:
                # 此处最好要加入segIdx，直接进行后续的分析
                # 需要记录位置
                mutAssignments[currentIdx].append(mut.sid)
                stripeNum += 1

            # Preorder traversal is consistent with printo_latex.py, meaning index
            # values should correspond to same vertices.
            pops[currentIdx] = {
                'cellular_prevalence': cellPrev,
                'stripe_number': stripeNum,
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
        return (pops, mutAssignments, structure)

    def _list_mutations(self, tree, segPool, includeStripeNames,
                        includeSegmentList=True):
        stripes = {}

        def _getSegInfo(segIdxL, segL):
            segs = []
            for idx in segIdxL:
                info = {
                    'name':segL[idx].name,
                    'chrom_name':segL[idx].chrom_name,
                    'start':segL[idx].start,
                    'end':segL[idx].end,
                    'tumor_read_number': segL[idx].tReadNum,
                    'normal_read_number': segL[idx].nNeadNum,
                    'gc': segL[idx].gc,
                }
                segs.append(info)

            return segs

        def _traverse(node):
            for mut in node['node'].get_data():
                stripes[mut.sid] = {
                    'tumor_read_number': mut.tReadNum,
                    'normal_read_number': mut.nNeadNum,
                    'copy_number': mut.copyNumber,
                    'genotype': mut.genotype,
                    'tag': mut.tag,
                }
                if includeStripeNames:
                    stripes[mut.sid]['name'] = mut.name
                if includeSegmentList:
                    stripes[mut.sid]['segment_index'] = mut.segIdxL.split(",")
                    stripes[mut.sid]['segment_info'] = _getSegInfo(mut.segIdxL.split(","), segL)


            for child in node['children']:
                _traverse(child)

        _traverse(tree.root)

        return { 'stripes': stripes }

    def __load_segPoolL_from_pkl(self, inputFilePath):

        inputFile = open(inputFilePath, 'rb')
        segPoolL = pkl.load(inputFile)

        return segPoolL[-1]

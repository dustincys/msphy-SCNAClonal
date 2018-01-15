import numpy as np
from collections import defaultdict

import sys
import os
sys.path.append(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import util2
import json


class ResultGenerator(object):
    def generate(self, tree_file):
        reader = util2.TreeReader(tree_file)
        first_tree = next(reader.load_trees())
        # cnv_logical_physical_mapping = json.loads(reader.read_extra_file('cnv_logical_physical_mapping.json'))
        try:
            params = json.loads(reader.read_extra_file('params.json'))
        except KeyError:
            # File not present in archive, likely because it originates from an older
            # run.
            params = {}
        reader.close()

        # 注意此处是否添加add stripe name
        mutlist = self._list_mutations(first_tree, include_stripe_names,
                                       include_segment_list)
        summaries = {}
        all_mutass = {}
        for idx, llh, pops, mutass, structure in self._summarize_all_pops(
                tree_file):
            summaries[idx] = {
                'llh': llh,
                'structure': structure,
                'populations': pops,
            }
            all_mutass[idx] = mutass

        return summaries, mutlist, all_mutass, params

    def _summarize_all_pops(self, tree_file):
        reader = util2.TreeReader(tree_file)
        for idx, llh, tree in reader.load_trees_and_metadata(
                remove_empty_vertices=True):
            yield (idx, llh) + self._summarize_pops(tree)
        reader.close()

    def _summarize_pops(self, tree):
        pops = {}
        structure = defaultdict(list)
        # Note that there will be an entry in mut_assignments for a given subclone
        # only if it has at least one SSM or CNV. This assumption holds true for all
        # PhyloWGS trees, so one can ascertain the number of cancerous populations
        # via len(mut_assignments).

        mut_assignments = defaultdict(lambda: [])
        idx = [0]

        def _traverse_r(vertex, parent):
            mutations = vertex.get_data()
            # vertex.params represents phis (i.e., population freqs) associated with
            # each sample.
            cell_prev = vertex.params
            current_idx = idx[0]

            num_stripes = 0
            for mut in mutations:
                mut_assignments[current_idx].append(mut.id)
                num_stripes += 1

            # Preorder traversal is consistent with printo_latex.py, meaning index
            # values should correspond to same vertices.
            pops[current_idx] = {
                'cellular_prevalence': cell_prev,
                'num_stripes': num_stripes,
            }

            # Visit children in order of decreasing phi.
            children = sorted(
                vertex.children(),
                key=lambda v: v.param,
                reverse=True)
            for child in children:
                idx[0] += 1
                structure[current_idx].append(idx[0])
                _traverse_r(child, current_idx)

        _traverse_r(tree.root['node'], None)
        return (pops, mut_assignments, structure)

    def _list_mutations(self, tree, include_stripe_names,
                        include_segment_list):
        stripes = {}

        def _traverse(node):
            for mut in node['node'].get_data():
                stripes[mut.stripe_id] = {
                    'tumor_reads_num': mut.tumor_reads_num,
                    'normal_reads_num': mut.normal_reads_num,
                    'baseline_label': mut.baseline_label,
                }
                if include_stripe_names:
                    stripes[mut.stripe_id]['name'] = mut.stripe_name
                if include_segment_list:
                    stripes[mut.stripe_id]['seg_idx'] = mut.segs_idx.split(",")

            for child in node['children']:
                _traverse(child)

        _traverse(tree.root)


        return { 'stripes': stripes }

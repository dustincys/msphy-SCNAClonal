from __future__ import print_function
import numpy as np


class ResultMunger(object):
    def __init__(self, tree_summaries, mutass, isStripe):
        self._tree_summaries = tree_summaries
        self._mutass = mutass
        self._isStripe = isStripe

    def _convert_keys_to_ints(self, dic):
        keys = dic.keys()
        for key in dic.keys():
            dic[int(key)] = dic[key]
            del dic[key]

    def remove_small_nodes(self, minSCNAs):
        for treeIdx, treeFeatures in self._tree_summaries.items():
            smallNodes = self._find_small_nodes(
                treeIdx, treeFeatures['populations'], minSCNAs)
            self._remove_nodes(smallNodes, treeIdx, mutDestination='best')

        return (self._tree_summaries, self._mutass)

    def remove_superclones(self):
        for treeIdx in self._tree_summaries.keys():
            pops = self._tree_summaries[treeIdx]['populations']
            structure = self._tree_summaries[treeIdx]['structure']

            rootIdx = 0
            assert len(structure[rootIdx]) > 0
            if len(structure[rootIdx]) != 1:
                # Ignore polyclonal trees, since it's not clear how to handle them.
                continue
            clonalIdx = structure[rootIdx][0]
            # Check assumption about node numbering (which isn't really important).
            assert clonalIdx == 1

            if clonalIdx not in structure or len(structure[clonalIdx]) != 1:
                # Either tree doesn't have subclones, or it has multiple subclones.
                # Ignore trees with >1 child coming off clonal node, since they don't
                # fit our definition of superclonal.
                continue
            childIdx = structure[clonalIdx][0]
            # Check assumption about node numbering (which isn't really important).
            assert childIdx == 2

            # At this point, "clone" is the superclone to remove, and "child" will be
            # the new clonal node.
            clone, child = pops[clonalIdx], pops[childIdx]
            if child['num_SCNAs'] == 0:
                # Prevent division by zero.
                continue
            # 注意，此处不将data个数作为考虑
            # if not (float(clone['num_ssms']) / child['num_ssms'] <= 0.33):
                # Child must have at least 3 times the number of SSMs as parent to be
                # considered superclonal tree.
                #print(treeIdx, 'differential in num_ssms too small', clone['num_ssms'], child['num_ssms'])
                continue
            if not (np.mean(
                    np.abs(
                        np.array(clone['cellular_prevalence']) - np.array(
                            child['cellular_prevalence']))) <= 0.1):
                # Cellular prevalences of clusters must be within 10%.
                #print(treeIdx, 'CPs too dissimilar', clone['num_ssms'], child['num_ssms'])
                continue

            print(treeIdx, 'Superclone:', clone, 'actual clone:', child)
            # Revise cellular prevalence to be weighted mean of both nodes.
            clonal_SCNAs, clonal_cp = clone['num_SCNAs'], np.array(
                clone['cellular_prevalence'])
            child_SCNAs, child_cp = child['num_SCNAs'], np.array(
                child['cellular_prevalence'])
            total_SCNAs = clonal_SCNAs + child_SCNAs
            # Cellular prevalence is a vector, so we must brought it into NumPy so we
            # can do scalar division on it. We must, however, convert it back to a
            # list so we can dump it to JSON.
            clone['cellular_prevalence'] = list(
                ((clonal_cp * clonal_SCNAs) +
                 (child_cp * child_SCNAs)) / total_SCNAs)

            # Remove the "true" clonal node and move its mutations to the superclonal
            # node. This is easier than the reverse, as it means that idx=1 always
            # points to the same node, and this function doesn't refer to any other
            # nodes. Otherwise, I'd be referring to the clonal idx=2 node, which
            # becomes idx=1 partway through.
            self._remove_nodes([childIdx], treeIdx, mutDestination='clonal')
            print(treeIdx, 'New clonal node:', clone)

    def remove_polyclonal_trees(self):
        polyidxs = set()

        for tidx in self._tree_summaries.keys():
            structure = self._tree_summaries[tidx]['structure']
            assert len(structure[0]) > 0
            if len(structure[0]) == 1:
                # Not polyclonal.
                continue
            polyidxs.add(tidx)

        polyclonal_frac = len(polyidxs) / float(len(self._tree_summaries))
        if polyclonal_frac >= 0.8:
            raise Exception(
                '%d%% of trees are polyclonal (%s of %s), so not enough to report good posterior.'
                % (100 * polyclonal_frac, len(polyidxs),
                   len(self._tree_summaries)))

        for pidx in sorted(polyidxs):
            print(pidx, 'polyclonal tree at idx=%s' % pidx)
            del self._tree_summaries[pidx]
            del self._mutass[pidx]

        assert set(self._tree_summaries.keys()) == set(self._mutass.keys())
        num_preceding_poly = 0
        for tidx in sorted(self._tree_summaries.keys()):
            if tidx in polyidxs:
                num_preceding_poly += 1
            else:
                if num_preceding_poly == 0:
                    continue
                newtidx = tidx - num_preceding_poly
                assert newtidx < tidx
                assert newtidx not in self._tree_summaries.keys(
                ) and newtidx not in self._mutass.keys()
                self._tree_summaries[newtidx] = self._tree_summaries[tidx]
                self._mutass[newtidx] = self._mutass[tidx]
                del self._tree_summaries[tidx]
                del self._mutass[tidx]

    def _renumber_nodes(self, treeIdx, subclone_idx_map):
        subclone_idxs = sorted(
            self._tree_summaries[treeIdx]['populations'].keys())

        num_removed = 0
        # We may have removed populations beyond max(subclone_idxs), but as these
        # occurred *after* the highest-indexed of the remaining populations,
        # renumbering is not necessary for them -- i.e., nothing in the tree is
        # affected by their removal.
        for subclone_idx in range(1, max(subclone_idxs) + 1):
            if subclone_idx not in self._tree_summaries[treeIdx][
                    'populations']:
                # Node was removed.
                num_removed += 1
            elif num_removed > 0:
                # Node not removed, but something before it was, so renumber.
                subclone_idx_map[subclone_idx] = subclone_idx - num_removed

        # By proceeding in sorted order, we guarantee we're not overwriting a
        # single element twice, which would give the wrong values. Why? Since the
        # new_idx of a node is always less than its original index, we guarantee
        # that we only move it "down". Since we proceed in sorted order, we have
        # already moved any nodes that preceded it, so we don't overwrite them.
        for subclone_idx in subclone_idxs:
            if subclone_idx not in subclone_idx_map:
                # No preceding nodes were removed, so do nothing.
                continue
            # Node remains, so must renumber it.
            new_idx = subclone_idx_map[subclone_idx]

            self._tree_summaries[treeIdx]['populations'][
                new_idx] = self._tree_summaries[treeIdx]['populations'][
                    subclone_idx]
            del self._tree_summaries[treeIdx]['populations'][subclone_idx]

            if subclone_idx in self._tree_summaries[treeIdx]['structure']:
                self._tree_summaries[treeIdx]['structure'][
                    new_idx] = self._tree_summaries[treeIdx]['structure'][
                        subclone_idx]
                del self._tree_summaries[treeIdx]['structure'][subclone_idx]

        # We must also renumber children in the structure -- just renumbering
        # parents isn't enough.
        for subclone_idx, children in self._tree_summaries[treeIdx][
                'structure'].items():
            self._tree_summaries[treeIdx]['structure'][subclone_idx] = [
                subclone_idx_map[c] if c in subclone_idx_map else c
                for c in children
            ]

    def _correct_mut_counts(self, populations, treeIdx):
        for sidx, subclone in populations.items():
            # Note that only mutass entries for subclones with (> 0 SSMs or > 0
            # CNVs) will exist. Thus, no mutass entry will exist for node 0, as it
            # never has SSMs or CNVs.
            if sidx == 0:
                continue
            subclone["num_SCNAs"] = len(self._mutass[treeIdx][sidx])

    def _remove_nodes(self, nodes, treeIdx, mutDestination):
        subclone_idx_map = {}
        treeFeatures = self._tree_summaries[treeIdx]

        # Remove summary stats about population.
        for node_idx in nodes:
            del treeFeatures['populations'][node_idx]
            # Mark node as removed. Use subclone_idx_map to track both node removals
            # and renumberings.
            subclone_idx_map[node_idx] = None

        self._remove_nodes_from_tree_structure(subclone_idx_map,
                                               treeFeatures['structure'])
        self._renumber_nodes(treeIdx, subclone_idx_map)
        self._reassign_muts(treeIdx, subclone_idx_map, mutDestination)
        self._correct_mut_counts(treeFeatures['populations'], treeIdx)

    def _find_small_nodes(self, treeIdx, populations, minSCNAs):
        smallNodes = set()

        subclone_idxs = sorted(populations.keys())
        last_phi = None
        last_idx = None

        if minSCNAs >= 1:
            # This is a count of SSMs, so use it without adjustment (but ensure it's an int).
            minSCNAs = int(minSCNAs)
        else:
            # This is a fraction of total SSMs.
            numData = 0
            currentSCNAPool = self._tree_summaries[treeIdx]['SCNAPool']
            if self._isStripe:
                numData = len(currentSCNAPool.stripes)
            else:
                numData = len(currentSCNAPool.segments)

            minSCNAs = int(round(float(minSCNAs) * numData))

        for subclone_idx in subclone_idxs:
            for p, children in self._tree_summaries[treeIdx][
                    'structure'].items():
                if subclone_idx in children:
                    parent = p
                    break
            subclone = populations[subclone_idx]
            # Ensure this node's phi is <= the phi of its preceding sibling node, if any exists.
            if subclone_idx > 0 and last_idx in self._tree_summaries[treeIdx]['structure'][parent]:
                assert subclone['cellular_prevalence'] <= last_phi
            last_phi = subclone['cellular_prevalence']
            last_idx = subclone_idx

            if subclone_idx == 0 or subclone['num_SCNAs'] >= minSCNAs:
                continue
            smallNodes.add(subclone_idx)

        return smallNodes

    def _remove_nodes_from_tree_structure(self, subclonal_idx_map,
                                          tree_structure):
        def _find_parent(struct, idx):
            for parent, children in struct.items():
                if idx in children:
                    return parent
            raise Exception('Could not find parent of %s in %s' % (idx,
                                                                   struct))

        removed = set([
            N for N in subclonal_idx_map.keys() if subclonal_idx_map[N] is None
        ])

        for rem in removed:
            parent = _find_parent(tree_structure, rem)
            # Remove node from parent
            tree_structure[parent] = [
                c for c in tree_structure[parent] if c != rem
            ]
            # Assign removed node's children to their grandparent
            if rem in tree_structure:
                tree_structure[parent] += tree_structure[rem]
                del tree_structure[rem]
            # Sort, since order may not be preserved.
            tree_structure[parent].sort()
            # If no children remain after deletion, remove child list from tree.
            if len(tree_structure[parent]) == 0:
                del tree_structure[parent]

    def _move_muts_to_best_node(self, muts, mutass, populations):
        for mut_type in ('ssms', 'cnvs'):
            for mut_id in muts[mut_type]:
                mut_stats = self._mutlist[mut_type][mut_id]
                ref_reads = np.mean(mut_stats['ref_reads'])
                total_reads = np.mean(mut_stats['total_reads'])
                # Note this doesn't take into account CNVs that may skew the
                # relationship between VAF and phi -- we assume that there is one
                # maternal and one paternal copy, and that only one of these is
                # mutated.
                implied_phi = 2 * (total_reads - ref_reads) / total_reads
                implied_phi = min(implied_phi, 1.0)

                lowest_phi_delta = 1
                best_node = None
                for pidx, pop in populations.items():
                    phi_delta = abs(
                        np.mean(pop['cellular_prevalence']) - implied_phi)
                    # Don't allow assignments to the non-cancerous root node.
                    if phi_delta < lowest_phi_delta and pidx != 0:
                        lowest_phi_delta = phi_delta
                        best_node = pidx

                mutass[best_node][mut_type].append(mut_id)


    def _move_muts_to_clonal_node(self, muts, mutass, populations, structure):
        rootIdx = 0
        clonalIdx = 1
        assert clonalIdx in structure[rootIdx]
        for mut_type in ('ssms', 'cnvs'):
            mutass[clonalIdx][mut_type] += muts[mut_type]

    def _reassign_muts(self, treeIdx, subclone_idx_map, destination='best'):
        deleted_muts = []
        mutass = self._mutass[treeIdx]

        for sidx in sorted(subclone_idx_map.keys()):
            new_idx = subclone_idx_map[sidx]
            # This ensures we're not improperly overwriting assignments.
            assert new_idx < sidx

            if new_idx is None:  # Node was removed.
                deleted_muts.append(mutass[sidx])
            else:
                mutass[new_idx] = mutass[sidx]
            del mutass[sidx]

        for dm in deleted_muts:
            if destination == 'best':
                self._move_muts_to_best_node(
                    dm, mutass, self._tree_summaries[treeIdx]['populations'])
            elif destination == 'clonal':
                self._move_muts_to_clonal_node(
                    dm, mutass, self._tree_summaries[treeIdx]['populations'],
                    self._tree_summaries[treeIdx]['structure'])
            else:
                raise Exception('Unknown destination: %s' % destination)

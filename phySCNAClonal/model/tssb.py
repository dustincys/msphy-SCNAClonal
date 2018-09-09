#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: tssb.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-12-06 15:39:03
#       History: phylowgs
# =============================================================================
'''
import sys
from time import *

import scipy.stats
from gwpy.segments import Segment, SegmentList
from numpy import *
from numpy.random import *

from phySCNAClonal.preprocess.utils import get_cn_allele_config
from phySCNAClonal.model.usegsampler.segsupportive import MultiRangeSampler
from util import betapdfln, boundbeta, logsumexp, sticks_to_edges

from copy import deepcopy


class TSSB(object):
    minDpAlpha = 1.0
    maxDpAlpha = 50.0
    minDpGamma = 1.0
    maxDpGamma = 10.0
    minAlphaDecay = 0.05
    maxAlphaDecay = 0.80

    def __init__(self,
                 dpAlpha=1.0,
                 dpGamma=1.0,
                 rootNode=None,
                 data=None,
                 minDepth=0,
                 maxDepth=15,
                 alphaDecay=1.0,
                 baseline=0,
                 maxCopyNumber=6):

        self.maxCopyNumber = maxCopyNumber
        # 此处的baseline是指什么
        self.baseline = baseline
        self.alleleConfig = get_cn_allele_config(maxCopyNumber)

        if rootNode is None:
            raise Exception("Root node must be specified.")

        self.minDepth = minDepth
        self.maxDepth = maxDepth
        self.dpAlpha = dpAlpha
        self.dpGamma = dpGamma
        self.alphaDecay = alphaDecay
        self.data = data
        self.dataNum = 0 if data is None else len(
            data)  # data.shape[0] #shankar
        self.root = {
            'node': rootNode,
            'main': boundbeta(1.0, dpAlpha) if self.minDepth == 0 else 0.0,
            'sticks': empty((0, 1)),
            'children': [],
            'tag': False
        }
        rootNode.tssb = self

        if False:
            dataU = rand(self.dataNum)
            self.assignments = []
            for n in range(self.dataNum):
                (c, path) = self.find_node(dataU[n])
                c.add_datum(n)
                self.assignments.append(c)
        else:
            self.assignments = []
            for n in range(self.dataNum):
                self.root['node'].add_datum(n)
                self.assignments.append(self.root['node'])

    def add_data(self, data):
        (weights, nodes) = self.get_mixture()
        newDataNum = len(data)  # data.shape[0] #shankar
        for n in range(newDataNum):
            logprobs = []
            for k, node in enumerate(nodes):
                logprobs.append(log(weights[k]) + node.logprob(data[n],
                                                               self.alleleConfig,
                                                               self.baseline,
                                                               self.maxCopyNumber))
            logprobs = array(logprobs)
            probs = exp(logprobs - logsumexp(logprobs))
            best_k = sum(rand() > cumsum(probs))
            nodes[best_k].add_datum(n + self.dataNum)
            self.assignments.append(nodes[best_k])
        self.data = vstack([self.data, data])
        self.dataNum += newDataNum

        # shankar

    #    def clear_data(self):
    #        dims = self.data.shape[1]
    #        for n in range(self.dataNum):
    #            self.assignments[n].remove_datum(n)
    #        self.assignments = []
    #        self.data        = empty((0,dims))
    #        self.dataNum    = 0

    def resample_node_params(self, iters=1):
        for iter in range(iters):

            def descend(root):
                for index, child in enumerate(root['children']):
                    descend(child)
                root['node'].resample_params()

            descend(self.root)

    def resample_assignments(self, tag):
        def path_lt(path1, path2):
            if len(path1) == 0 and len(path2) == 0:
                return 0
            elif len(path1) == 0:
                return 1
            elif len(path2) == 0:
                return -1
            s1 = "".join(map(lambda i: "%03d" % (i), path1))
            s2 = "".join(map(lambda i: "%03d" % (i), path2))

            return cmp(s2, s1)

        epsilon = finfo(float64).eps
        # this is not useful
        lengths = []

        # change data range
        # 需要保存索引号码

        tagL = array([int(item.tag) for item in self.data])
        timeTagSeq = sort(unique(tagL))
        currentTimeTagIdx = where(timeTagSeq == tag)[0][0]
        uSegL = MultiRangeSampler(0,1)
        if currentTimeTagIdx == 0:
            self.reset_time_tag()
        else:
            self.mark_time_tag(timeTagSeq[currentTimeTagIdx - 1])
            uNegtive = self.get_u_segL()
            uSegL.remove(uNegtive)


        for n in where(tagL == tag)[0]:
            # change to segment operation
            tempUSegL = deepcopy(uSegL)
            minU = tempUSegL.lowerBoundary
            maxU = tempUSegL.upperBoundary

            llhMapD = {}
            # Get an initial uniform variate.
            ancestors = self.assignments[n].get_ancestors()
            current = self.root
            indices = []
            for anc in ancestors[1:]:
                index = map(lambda c: c['node'],
                            current['children']).index(anc)
                current = current['children'][index]
                indices.append(index)

            # Here indices could be used as path to locate v stick

            currentVStick = self._locate_v_stick(indices)
            isInNegativeSpace = currentVStick['tag']

            llhS = -float('Inf')
            if isInNegativeSpace:
                oldLlh = -float('Inf')
                llhS = -float('Inf')
            else:
                oldLlh = self.assignments[n].logprob(self.data[n:n + 1],
                                                     self.alleleConfig,
                                                     self.baseline,
                                                     self.maxCopyNumber)
                llhMapD[self.assignments[n]] = oldLlh
                llhS = log(rand()) + oldLlh

            while True:
                newU = tempUSegL.sample()
                (newNode, newPath) = self.find_node(newU)
                if newNode.parent() is None:
                    # shankar: to make root node empty
                    newNode = newNode.children()[0]
                    newPath = [0]
                oldNode = self.assignments[n]
                oldNode.remove_datum(n)
                newNode.add_datum(n)
                self.assignments[n] = newNode
                if newNode in llhMapD:
                    newLlh = llhMapD[newNode]
                else:
                    ####################################
                    #  Record most likely copy number  #
                    ####################################
                    newLlh = newNode.logprob(self.data[n:n + 1],
                                             self.alleleConfig, self.baseline,
                                             self.maxCopyNumber)
                    if -float('Inf') == newLlh:
                        print >> sys.stderr, "Slice sampler weirdness"
                    llhMapD[newNode] = newLlh
                if newLlh > llhS:
                    break
                elif abs(maxU - minU) < epsilon:
                    if -float('Inf') == newLlh:
                        raise Exception("Slice sampler weirdness.")

                    newNode.remove_datum(n)
                    oldNode.add_datum(n)
                    self.assignments[n] = oldNode
                    print >> sys.stderr, "Slice sampler shrank down.  Keep current state."
                    break
                else:
                    newNode.remove_datum(n)
                    oldNode.add_datum(n)
                    self.assignments[n] = oldNode
                    pathComp = path_lt(indices, newPath)
                    if pathComp < 0:
                        tempUSegL.removeLeft(newU)
                    elif pathComp >= 0:  # temporary fix only!!!!!!
                        tempUSegL.removeRight(newU)
                    else:
                        raise Exception("Slice sampler weirdness.")
                    minU = tempUSegL.lowerBoundary
                    maxU = tempUSegL.upperBoundary
            lengths.append(len(newPath))
        lengths = array(lengths)

    def _locate_v_stick(self, indices):
        """locate and return v stick by indices

        :indices: The indices to current data node
        :returns: v stick

        """
        target = self.root
        for index in indices:
            target = target['children'][index]

        return target

    def cull_tree(self):
        def descend(root):
            counts = array(map(lambda child: descend(child), root['children']))
            keep = len(trim_zeros(counts, 'b'))

            for child in root['children'][keep:]:
                child['node'].kill()
                del child['node']

            root['sticks'] = root['sticks'][:keep]
            root['children'] = root['children'][:keep]

            return sum(counts) + root['node'].num_local_data()

        descend(self.root)

    def resample_sticks(self):
        def descend(root, depth=0):

            dataDown = 0
            indices = range(len(root['children']))
            indices.reverse()
            for i in indices:
                child = root['children'][i]
                childData = descend(child, depth + 1)
                postAlpha = 1.0 + childData
                postBeta = self.dpGamma + dataDown
                root['sticks'][i] = boundbeta(
                    postAlpha,
                    postBeta) if depth != 0 else .999999  # shankar
                dataDown += childData

            # Resample the main break.
            dataHere = root['node'].num_local_data()
            postAlpha = 1.0 + dataHere
            postBeta = (self.alphaDecay**depth) * self.dpAlpha + dataDown
            root['main'] = boundbeta(
                postAlpha, postBeta) if self.minDepth <= depth else 0.0
            if depth == 0:
                root['main'] = 1e-30  # to make root node empty (shankar)

            return dataHere + dataDown

        descend(self.root)

    def resample_stick_orders(self):
        def descend(root, depth=0):
            if not root['children']:
                return

            newOrder = []
            represented = set(
                filter(lambda i: root['children'][i]['node'].has_data(),
                       range(len(root['children']))))
            # 每一phi stick 的实际长度
            allWeights = diff(hstack([0.0, sticks_to_edges(root['sticks'])]))
            while True:
                if not represented:
                    break

                u = rand()
                while True:
                    subIndices = filter(lambda i: i not in newOrder,
                                         range(root['sticks'].shape[0]))
                    # 此处添加了剩余空间的长度
                    subWeights = hstack(
                        [allWeights[subIndices], 1.0 - sum(allWeights)])
                    # 每一个空间所占用的比重
                    subWeights = subWeights / sum(subWeights)
                    # 随机获得一个位置，此位置之前完整空间的个数
                    index = sum(u > cumsum(subWeights))

                    if index == len(subIndices):
                        root['sticks'] = vstack(
                            [root['sticks'],
                             boundbeta(1, self.dpGamma)])
                        root['children'].append({
                            'node': root['node'].spawn(),
                            'main': boundbeta(1.0, (self.alphaDecay ** (depth + 1)) * self.dpAlpha) if self.minDepth <= (depth + 1) else 0.0, # 此处minDepth 应该是手动控制的
                            'sticks': empty((0, 1)),
                            'children': [],
                            'tag': False
                        })
                        allWeights = diff(
                            hstack([0.0, sticks_to_edges(root['sticks'])]))
                    else:
                        index = subIndices[index]
                        break
                newOrder.append(index)
                represented.discard(index)

            newChildren = []
            for k in newOrder:
                child = root['children'][k]
                newChildren.append(child)
                descend(child, depth + 1)

            for k in filter(lambda k: k not in newOrder,
                            range(root['sticks'].shape[0])):
                root['children'][k]['node'].kill()
                del root['children'][k]['node']

            root['children'] = newChildren
            root['sticks'] = zeros((len(root['children']), 1))

        descend(self.root)

        # Immediately resample sticks.
        self.resample_sticks()

    def resample_hypers(self, dpAlpha=True, alphaDecay=True, dpGamma=True):
        def dp_alpha_llh(dpAlpha, alphaDecay):
            def descend(dpAlpha, root, depth=0):
                llh = betapdfln(root['main'], 1.0, (alphaDecay**depth) *
                                dpAlpha) if self.minDepth <= depth else 0.0
                for child in root['children']:
                    llh += descend(dpAlpha, child, depth + 1)
                return llh

            return descend(dpAlpha, self.root)

        if dpAlpha:
            upper = self.maxDpAlpha
            lower = self.minDpAlpha
            llhS = log(rand()) + dp_alpha_llh(self.dpAlpha, self.alphaDecay)
            while True:
                newDpAlpha = (upper - lower) * rand() + lower
                newLlh = dp_alpha_llh(newDpAlpha, self.alphaDecay)
                if newLlh > llhS:
                    break
                elif newDpAlpha < self.dpAlpha:
                    lower = newDpAlpha
                elif newDpAlpha > self.dpAlpha:
                    upper = newDpAlpha
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.dpAlpha = newDpAlpha

        if alphaDecay:
            upper = self.maxAlphaDecay
            lower = self.minAlphaDecay
            llhS = log(rand()) + dp_alpha_llh(self.dpAlpha, self.alphaDecay)
            while True:
                newAlphaDecay = (upper - lower) * rand() + lower
                newLlh = dp_alpha_llh(self.dpAlpha, newAlphaDecay)
                if newLlh > llhS:
                    break
                elif newAlphaDecay < self.alphaDecay:
                    lower = newAlphaDecay
                elif newAlphaDecay > self.alphaDecay:
                    upper = newAlphaDecay
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.alphaDecay = newAlphaDecay

        def dp_gamma_llh(dpGamma):
            def descend(dpGamma, root):
                llh = 0
                for i, child in enumerate(root['children']):
                    llh += betapdfln(root['sticks'][i], 1.0, dpGamma)
                    llh += descend(dpGamma, child)
                return llh

            return descend(dpGamma, self.root)

        if dpGamma:
            upper = self.maxDpGamma
            lower = self.minDpGamma
            llhS = log(rand()) + dp_gamma_llh(self.dpGamma)
            while True:
                newDpGamma = (upper - lower) * rand() + lower
                newLlh = dp_gamma_llh(newDpGamma)
                if newLlh > llhS:
                    break
                elif newDpGamma < self.dpGamma:
                    lower = newDpGamma
                elif newDpGamma > self.dpGamma:
                    upper = newDpGamma
                else:
                    raise Exception("Slice sampler shrank to zero!")
            self.dpGamma = newDpGamma

    def draw_data(self, dataNum=1, **args):
        self.data = []
        self.assignments = []
        for n in range(dataNum):
            u = rand()
            (node, path) = self.find_node(u)
            self.data.append(node.sample(args))
            self.assignments.append(node)
            node.add_datum(n)
            self.dataNum += 1
        self.data = concatenate(self.data)
        return self.data

    def resample_data(self, **args):
        for n in range(self.dataNum):
            u = rand()
            (node, path) = self.find_node(u)
            self.assignments[n].remove_datum(n)
            node.add_datum(n)
            self.assignments[n] = node
            self.data[n] = node.sample(args)[0]

    def mark_time_tag(self, tag):
        """
        mark each node's time tag status.
        All the nodes in the tssb tree, that contains offspings data with tag <=
        given tag.
        """
        def descend(root, tag):
            if 0 == len(root['children']):
                root['tag'] = False

                timeTags = [int(item.tag) for item in root['node'].get_data() if
                            int(item.tag) <= tag]
                if 0 < len(timeTags):
                    return True
                else:
                    return False
            else:
                if 0 < sum([descend(child, tag) for child in root['children']]):
                    root['tag'] = True
                else:
                    root['tag'] = False

                if root['tag']:
                    return True
                else:
                    timeTags = [int(item.tag) for item in root['node'].get_data() if
                                int(item.tag) <= tag]
                    if 0 < len(timeTags):
                        return True
                    else:
                        return False

        descend(self.root, tag)

    def reset_time_tag(self):

        def descend(root):
            root['tag'] = False
            for child in root['children']:
                descend(child)

        descend(self.root)

    def find_node(self, u):
        def descend(root, u, depth=0):
            if depth >= self.maxDepth:
                # print >>sys.stderr, "WARNING: Reached maximum depth."
                return (root['node'], [])
            elif u < root['main']:
                if root['tag']:
                     print >>sys.stderr, "Negative space located!!."

                return (root['node'], [])
            else:
                # Rescale the uniform variate to the remaining interval.
                u = (u - root['main']) / (1.0 - root['main'])

                # Perhaps break sticks out appropriately.
                if depth > 0:
                    while not root['children'] or (
                            1.0 - prod(1.0 - root['sticks'])) < u:
                        root['sticks'] = vstack([
                            root['sticks'],
                            # 注意此处为右边界
                            boundbeta(1, self.dpGamma) if depth != 0 else .999
                        ])  # shankar
                        root['children'].append({
                            'node':
                            root['node'].spawn(),
                            'main':
                            boundbeta(1.0, (self.alphaDecay**
                                            (depth + 1)) * self.dpAlpha)
                            if self.minDepth <= (depth + 1) else 0.0,
                            'sticks':
                            empty((0, 1)),
                            'children': [],
                            'tag': False
                        })

                    edges = 1.0 - cumprod(1.0 - root['sticks'])
                    index = sum(u > edges)
                    edges = hstack([0.0, edges])
                    u = (u - edges[index]) / (edges[index + 1] - edges[index])

                    (node, path) = descend(root['children'][index], u,
                                           depth + 1)
                else:
                    index = 0
                    (node, path) = descend(root['children'][index], u,
                                           depth + 1)

                path.insert(0, index)

                return (node, path)

        return descend(self.root, u)

    def get_u_segL(self):
        """
        return the non-supportive range of u in the find_node function.
        here, all the target node is pre-marked in the mark_time_tag function.
        """
        def descend(root):
            if not root['tag']:
                return None, None
            elif 0 == len(root['children']):
                return array([0]), array(root['main'])
            else:
                # get edges
                starts = array([])
                ends = array([])

                edges = 1.0 - cumprod(1.0 - root['sticks'])
                edges = hstack([0.0, edges])
                for index, child in zip(range(len(root['children'])),
                                        root['children']):
                    startsChild, endsChild = descend(child)
                    if startsChild is None or endsChild is None:
                        continue
                    else:
                        startsChild = startsChild *\
                            (edges[index + 1] - edges[index]) + edges[index]
                        endsChild = endsChild *\
                            (edges[index + 1] - edges[index]) + edges[index]
                        starts = r_[starts, startsChild]
                        ends = r_[ends, endsChild]

                starts = starts * (1 - root['main']) + root['main']
                ends = ends * (1 - root['main']) + root['main']

                starts = r_[starts, 0]
                ends = r_[ends, root['main']]
                return starts, ends

        starts, ends = descend(self.root)
        assert(len(starts) == len(ends))

        sl = SegmentList([])
        for i in range(len(starts)):
            sl.append(Segment(starts[i], ends[i]))
        sl.coalesce()

        return sl


    def get_nodes(self):
        def descend(root):
            node = [root['node']]
            for child in root['children']:
                child_nodes = descend(child)
                node.extend(child_nodes)
            return node

        return descend(self.root)

    def get_mixture(self):
        def descend(root, mass):
            weight = [mass * root['main']]
            node = [root['node']]
            edges = sticks_to_edges(root['sticks'])
            weights = diff(hstack([0.0, edges]))

            for i, child in enumerate(root['children']):
                (child_weights, child_nodes) = descend(
                    child, mass * (1.0 - root['main']) * weights[i])
                weight.extend(child_weights)
                node.extend(child_nodes)
            return (weight, node)

        # 返回两个向量
        return descend(self.root, 1.0)

    def complete_data_log_likelihood(self):
        weights, nodes = self.get_mixture()
        llhs = []
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.num_local_data() * log(weights[i]) +
                            node.data_log_likelihood(self.alleleConfig,
                                                     self.baseline,
                                                     self.maxCopyNumber))
        return sum(array(llhs))

    def complete_log_likelihood(self):
        weights, nodes = self.get_mixture()
        llhs = [
            self.dp_alpha_llh(self.dpAlpha, self.alphaDecay),
            self.dp_gamma_llh(self.dpGamma)
        ]
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.data_log_likelihood(self.alleleConfig,
                                                     self.baseline,
                                                     self.maxCopyNumber))
        return sum(array(llhs))

    def dp_alpha_llh(self, dpAlpha, alphaDecay):
        def descend(dpAlpha, root, depth=0):
            llh = betapdfln(root['main'], 1.0, (alphaDecay**depth) *
                            dpAlpha) if self.minDepth <= depth else 0.0
            for child in root['children']:
                llh += descend(dpAlpha, child, depth + 1)
            return llh

        return descend(dpAlpha, self.root)

    def dp_gamma_llh(self, dpGamma):
        def descend(dpGamma, root):
            llh = 0
            for i, child in enumerate(root['children']):
                llh += betapdfln(root['sticks'][i], 1.0, dpGamma)
                llh += descend(dpGamma, child)
            return llh

        return descend(dpGamma, self.root)

    def print_graph(self, fh, base_width=5000, min_width=5):
        print >> fh, """graph: { title:            "TSSB Graph"  \
                                portsharing:      no            \
                                smanhattanedges:  yes           \
                                equalydist:       yes           \
                                layout_algorithm: tree          \
                                node.fontname:    "helvR8"      \
                                node.height:      25 """


        print >> fh, """node: { label:"%0.5f" title:"%s" width:%d}""" % (
            self.root['main'], "X",
            max(int(self.root['main'] * base_width), min_width))

        def descend(root, name, mass):
            total = 0.0
            edges = sticks_to_edges(root['sticks'])
            weights = diff(hstack([0.0, edges]))
            for i, child in enumerate(root['children']):
                childName = "%s-%d" % (name, i)
                childMass = mass * weights[i] * child['main']
                print >> fh, """node: {  label:"%0.5f" title:"%s" width:%d}""" % (
                    childMass, childName,
                    max(int(childMass * base_width), min_width))
                print >> fh, """edge: { source:"%s" target:"%s" anchor:1}""" % (
                    name, childName)
                total += childMass + descend(child, childName,
                                              mass * weights[i] *
                                              (1.0 - child['main']))
            return total

        print >> fh, """}"""

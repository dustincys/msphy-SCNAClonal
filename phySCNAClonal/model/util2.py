#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
from numpy import *
import cPickle as pickle
import zipfile
import shutil

import scipy.stats as stat
from scipy.stats import beta, binom
from scipy.special import gammaln
from math import exp, log

import csv
# Allow long lines in .cnv files, which can potentially list thousands of SSMs
# in one CNV. According to http://stackoverflow.com/a/15063941, this value can
# be as much as a C long. C longs are guaranteed to accommodate at least this
# value, which is a signed 32-bit int.
csv.field_size_limit(2147483647)

from tssb import *


def log_factorial(n):
    return gammaln(n + 1)


def log_bin_coeff(n, k):
    return log_factorial(n) - log_factorial(k) - log_factorial(n - k)


def log_binomial_likelihood(x, n, mu):
    return x * log(mu) + (n - x) * log(1 - mu)


def log_beta(a, b):
    return gammaln(a) + gammaln(b) - gammaln(a + b)


def logsumexp(X, axis=None):
    maxes = numpy.max(X, axis=axis)
    return numpy.log(numpy.sum(numpy.exp(X - maxes), axis=axis)) + maxes


def parse_physical_cnvs(pcnvs):
    physicalCnvs = []

    for physical_cnv in pcnvs.split(';'):
        fields = physical_cnv.split(',')
        cnv = dict([F.split('=', 1) for F in fields])
        for key in ('start', 'end', 'major_cn', 'minor_cn'):
            cnv[key] = int(cnv[key])
        cnv['cell_prev'] = [float(C) for C in cnv['cell_prev'].split('|')]
        physicalCnvs.append(cnv)

    return physicalCnvs


def load_data(inputFilePath):
    # load stripes data

    inputFile = open(inputFilePath, 'rb')
    dataStripes = pkl.load(inputFile)

    return (dataStripes.stripes, dataStripes.baseline)

#################################################
# some useful functions to get some info about,
# the tree, used by CNV related computations


def set_node_height(tssb):
    tssb.root['node'].ht = 0

    def descend(root, ht):
        for child in root.children():
            child.ht = ht
            descend(child, ht+1)
    descend(tssb.root['node'], 1)


def set_path_from_root_to_node(tssb):
    nodes = tssb.get_nodes()
    for node in nodes:
        node.path = node.get_ancestors()


def map_datum_to_node(tssb):
    nodes = tssb.get_nodes()
    for node in nodes:
        for datum in node.get_data():
            datum.node = node
#################################################


def check_bounds(p, l=0.0001, u=.9999):
    if p < l:
        p = l
    if p > u:
        p = u
    return p

# removes the empty nodes from the tssb tree
# Does not removes root as it is not required
# root: root of the current tree
# parent: parent of the root
# Note this funciton modifies the sticks so they remain valid.


def remove_empty_nodes(root, parent=None):
    for child in list(root['children']):
        remove_empty_nodes(child, root)
    if (root['node'].get_data() == []):
        if (root['children'] == []):  # leaf
            if (parent is not None):
                ind = parent['children'].index(root)
                parent['children'].remove(root)
                root['node'].kill()
                parent['sticks'] = delete(parent['sticks'], ind, 0)
            return
        else:
            if (parent is not None):
                parent_ = root['node'].parent()
                ind = parent['children'].index(root)
                for i, child in enumerate(list(root['children'])):
                    parent['children'].append(child)
                    toappend = zeros((1, 1))
                    toappend[0] = root['sticks'][i]
                    parent['sticks'] = append(parent['sticks'], toappend, 0)
                    root['children'].remove(child)
                for child in list(root['node'].children()):
                    child._parent = parent_
                    parent_.add_child(child)
                    root['node'].remove_child(child)
                parent['children'].remove(root)
                parent['sticks'] = delete(parent['sticks'], ind, 0)
                root['node'].kill()


def rm_safely(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno == 2:  # Ignore "no such file" errors
            pass
        else:
            raise e


class CorruptZipFileError(Exception):
    pass


class BackupManager(object):
    def __init__(self, filenames):
        self._filenames = filenames
        self._backupFileNameL = [os.path.realpath(fn) +
            '.backup' for fn in self._filenames]

    def save_backup(self):
        for fn, backupFn in zip(self._filenames, self._backupFileNameL):
            # 文件复制
            shutil.copy2(fn, backupFn)

    def restore_backup(self):
        for fn, backupFn in zip(self._filenames, self._backupFileNameL):
            shutil.copy2(backupFn, fn)

    def remove_backup(self):
        for backupFn in self._backupFileNameL:
            try:
                os.remove(backupFn)
            except OSError:
                pass


class StateManager(object):
    defaultLastStateFn = 'state.last.pickle'
    defaultInitialStateFn = 'state.initial.pickle'

    def __init__(self):
        self._initialStateFn = StateManager.defaultInitialStateFn
        self._lastStateFn = StateManager.defaultLastStateFn

    def _write_state(self, state, stateFn):
        with open(stateFn, 'w') as stateFile:
            pickle.dump(state, stateFile, protocol=pickle.HIGHEST_PROTOCOL)

    def write_state(self, state):
        self._write_state(state, self._lastStateFn)

    def load_state(self):
        with open(self._lastStateFn) as stateFile:
            return pickle.load(stateFile)

    def load_initial_state(self):
        with open(self._initialStateFn) as stateFile:
            return pickle.load(stateFile)

    def write_initial_state(self, state):
        self._write_state(state, self._initialStateFn)

    def delete_state_file(self):
        rm_safely(self._lastStateFn)

    def state_exists(self):
        return os.path.isfile(self._lastStateFn)


class TreeWriter(object):
    defaultArchiveFn = 'trees.zip'

    def __init__(self, resumeRun=False):
        self._archiveFn = TreeWriter.defaultArchiveFn
        if resumeRun:
            self._ensure_archive_is_valid()
        else:
            # Remove file to avoid unwanted behaviour. By the zipfile module's
            # behaviour, given that we open the file with the "a" flag, if a
            # non-zip file exists at this path, a zip file will be appended to
            # the file; otherwise, if the file is already a zip, additional
            # files will be written into the zip. On a new run, neither case is
            # something we want.
            rm_safely(self._archiveFn)

    def add_extra_file(self, filename, data):
        self._open_archive()
        self._archive.writestr(filename, data)
        self._close_archive()

    def _ensure_archive_is_valid(self):
        with zipfile.ZipFile(self._archiveFn) as zipf:
            if zipf.testzip() is not None:
                raise CorruptZipFileError(
                    'Corrupt zip file: %s' %
                    self._archiveFn)

    def _open_archive(self):
        self._archive = zipfile.ZipFile(
            self._archiveFn,
            'a',
            compression=zipfile.ZIP_DEFLATED,
            allowZip64=True)

    def _close_archive(self):
        self._archive.close()

    def write_trees(self, serializedTrees):
        self._open_archive()
        for st, idx, llh in serializedTrees:
            isBurnin = idx < 0
            prefix = isBurnin and 'burnin' or 'tree'
            treefn = '%s_%s_%s' % (prefix, idx, llh)
            self._archive.writestr(treefn, st)
        self._close_archive()


class TreeReader(object):
    def __init__(self, archive_fn):
        self._archive = zipfile.ZipFile(archive_fn)
        infolist = self._archive.infolist()
        treeInfoL = [t for t in infolist if t.filename.startswith('tree_')]
        burninInfoL = [t for t in infolist if t.filename.startswith('burnin_')]

        # Sort by index
        treeInfoL.sort(key=lambda tinfo: self._extract_metadata(tinfo)[0])
        burninInfoL.sort(key=lambda tinfo: self._extract_burnin_idx(tinfo))

        self._treeL = []
        self._burninTreeL = []

        for info in treeInfoL:
            idx, llh = self._extract_metadata(info)
            assert idx == len(self._treeL)
            self._treeL.append((idx, llh, info))
        for info in burninInfoL:
            idx = self._extract_burnin_idx(info)
            assert len(burninInfoL) + idx == len(self._burninTreeL)
            self._burninTreeL.append((idx, info))

    def read_extra_file(self, filename):
        return self._archive.read(filename)

    def num_trees(self):
        return len(self._treeL)

    def close(self):
        self._archive.close()

    def _extract_metadata(self, zinfo):
        tokens = zinfo.filename.split('_')
        idx = int(tokens[1])
        llh = float(tokens[2])
        return (idx, llh)

    def _extract_burnin_idx(self, zinfo):
        idx = int(zinfo.filename.split('_')[1])
        return idx

    def _parse_tree(self, zinfo, removeEmptyVertices=False):
        pickled = self._archive.read(zinfo)
        tree = pickle.loads(pickled)
        if removeEmptyVertices:
            remove_empty_nodes(tree.root)
        return tree

    def load_tree(self, idx, removeEmptyVertices=False):
        tidx, llh, zinfo = self._treeL[idx]
        assert tidx == idx
        return self._parse_tree(zinfo, removeEmptyVertices)

    def load_trees(self, numTrees=None, removeEmptyVertices=False):
        for idx, llh, tree in self.load_trees_and_metadata(
                numTrees, removeEmptyVertices):
            yield tree

    def load_trees_and_burnin(self, removeEmptyVertices=False):
        for tidx, zinfo in self._burninTreeL:
            tree = self._parse_tree(zinfo, removeEmptyVertices)
            yield (tidx, tree)
        for tidx, llh, zinfo in self._treeL:
            tree = self._parse_tree(zinfo, removeEmptyVertices)
            yield (tidx, tree)

    def load_trees_and_metadata(
            self,
            numTrees=None,
            removeEmptyVertices=False):
        # Sort by LLH
        trees = sorted(
            self._treeL,
            key=lambda tidx_llh_zinfo: tidx_llh_zinfo[1],
            reverse=True)

        if numTrees is not None:
            numTrees = min(numTrees, len(trees))
            trees = trees[:numTrees]

        for tidx, llh, zinfo in trees:
            tree = self._parse_tree(zinfo, removeEmptyVertices)
            yield (tidx, llh, tree)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: params.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-12-08 16:43:47
#       History: phylowgs
# =============================================================================
'''
###### code to sample from the paramater posterior p(\phi | data) ########

import os
import subprocess as sp

import numpy
from numpy import *
from numpy.random import dirichlet

import util2 as u2
from tssb import *
from util import dirichletpdfln


def get_c_fnames(tmp_dir):
    def _make_c_fname(name):
        fname = 'c_%s.txt' % (name)
        return os.path.join(tmp_dir, fname)

    FNAME_C_TREE = _make_c_fname('tree')
    FNAME_C_PARAMS = _make_c_fname('params')
    FNAME_C_MH_ARATIO = _make_c_fname('mh_ar')

    return (FNAME_C_TREE, FNAME_C_PARAMS, FNAME_C_MH_ARATIO)


# done for multi-sample


def metropolis(tssb,
               iters=1000,
               std=0.01,
               burnin=0,
               n_stripes=0,
               fin='',
               rseed=1,
               tmp_dir='.'):
    wts, nodes = tssb.get_mixture()

    # file names
    FNAME_STRIPE_DATA = fin
    FNAME_C_TREE, FNAME_C_PARAMS, FNAME_C_MH_ARATIO = get_c_fnames(tmp_dir)

    ## initialize the MH sampler###########
    # sample_cons_params(tssb)
    # update_params(tssb)
    ######################################

    ## prepare to call the c++ code ###########
    u2.set_node_height(tssb)
    write_tree(tssb, FNAME_C_TREE)  # write the current tree to the disk
    u2.map_datum_to_node(tssb)
    ###########################################

    MH_ITR = str(iters)
    MH_STD = str(std)
    N_STRIPE_DATA = str(n_stripes)
    NNODES = str(len(nodes))
    TREE_HEIGHT = str(max([node.ht for node in nodes]) + 1)

    script_dir = os.path.dirname(os.path.realpath(__file__))
    sp.check_call([
        '%s/mh.o' % script_dir,
        MH_ITR,
        MH_STD,
        N_STRIPE_DATA,
        NNODES,
        TREE_HEIGHT,
        FNAME_STRIPE_DATA,
        FNAME_C_TREE,
        FNAME_C_PARAMS,
        FNAME_C_MH_ARATIO])
    ar = str(loadtxt(FNAME_C_MH_ARATIO, dtype='string'))
    # update the tree with the new parameters sampled using the c++ code
    update_tree_params(tssb, FNAME_C_PARAMS)

    return ar


# done for multi-sample

def write_tree(tssb, fname):
    fh = open(fname, 'w')
    wts, nodes = tssb.get_mixture()
    didIntD = dict()

    # 此处stripe的ID应为数值
    # 此处ID 最好重新计算
    # 与sid分开处理
    # sid 记录类别信息，类别中的tag信息，是否是基线信息
    for dat in tssb.data:
        didIntD[dat.id] = int(dat.id)

    def descend(root):
        for child in root.children():
            descend(child)

        # write data#
        cids = ''
        for child in root.children():
            cids += str(child.id) + ','
        cids = cids.strip(',')
        if cids == '':
            cids = str(-1)

        dids = ''
        for dat in root.get_data():
            dids += str(didIntD[dat.id]) + ','
        dids = dids.strip(',')
        if dids == '':
            dids = str(-1)

        line = str(root.id) + '\t' + str( root.param) + '\t' +\
            str(root.pi) + '\t' + str(len(root.children())) + '\t' +\
            cids + '\t' + str(len(root.get_data())) + '\t' + dids +\
            '\t' + str(root.ht)

        fh.write(line)
        fh.write('\n')
        fh.flush()
        ###############

    descend(tssb.root['node'])
    fh.flush()
    fh.close()

# no changes for multi-sample
# data/node state format (parameter independent dot-product weights)
# datum_id	node_id_1,pi,nr,nv;node_id_2,pi,nr,nv;....
# these weights are used to compute data log-likelihood


# 此处不需要传递拷贝数和基因型参数
# def write_data_state(tssb, fname):
    # fh = open(fname, 'w')
    # wts, nodes = tssb.get_mixture()

    # # 此处dat为stripe，记录拷贝数、基因型、以及population frequency
    # # 注意此处基因型无法区分呈现互补数的基因型
    # for dat in tssb.data:
        # if not dat.node:
            # continue  # todo: this won't happen
        # for node in nodes:
            # fh.write(str(dat.id) + '\t' + str(dat.copy_number) + '\t' +
                # str(dat.genotype) + '\t' + str(dat.phi))
            # fh.write('\n')

    # fh.flush()
    # fh.close()


# done for multi-sample

def find_most_recent_cnv(dat, nd):
    out = None
    for n in nd.get_ancestors()[::-1]:
        if n in [x[0].node for x in dat.cnv]:
            out = [x for x in dat.cnv if x[0].node == n][0]
            break
    return out


# done for multi sample


def update_tree_params(tssb, fname):
    wts, nodes = tssb.get_mixture()
    ndict = dict()
    for node in nodes:
        ndict[node.id] = node

    fh = open(fname)
    params = [line.split() for line in fh.readlines()]
    fh.close()

    for p in params:
        ndict[int(p[0])].param = float(p[1])
        ndict[int(p[0])].pi = float(p[2])
    # params=loadtxt('c_params.txt')
    # for p in params:
    #	ndict[p[0]].params = p[1]
    #	ndict[p[0]].pi = p[2]


# def string_to_list(p):
    # p = p.strip(',')
    # return array([float(pp) for pp in p.split(',')])


# done for multi-sample
# tree-structured finite-dimensional stick breaking


def sample_cons_params(tssb):
    def descend(root):

        if root.parent() is None:
            root.param1 = 1
            root.pi1 = root.param1 * rand(1)  # break nu stick
        r = root.param1 - root.pi1  # mass assigned to children
        p = rand(len(root.children()))
        p = r * p * 1. / sum(p)
        index = 0
        for child in root.children():
            child.param1 = p[index]  # break psi sticks
            # break nu stick
            child.pi1 = child.param1 * (rand(1)**
                                                 (len(child.children()) > 0))
            index += 1
        for child in root.children():
            descend(child)

    descend(tssb.root['node'])



def update_params(tssb):
    def descend(root):
        for child in root.children():
            descend(child)
        root.param = root.param1
        root.pi = root.pi1

    descend(tssb.root['node'])


###### old code, not in use #############
# data/node state format (parameter independent dot-product weights)
# datum_id	node_id_1,pi,nr,nv;node_id_2,pi,nr,nv;....
# these weights are used to compute data log-likelihood
def write_data_state1111(tssb):
    fh = open('c_data_states.txt', 'w')
    wts, nodes = tssb.get_mixture()

    for dat in tssb.data:
        if not dat.cnv:
            continue  # nothing to do for CNVs

        if not dat.node:
            continue  # todo: this won't happen

        ancestors = dat.node.get_ancestors()  # path from root to ssm node

        mr_cnv = dat.cnv[0]  # CNV corresponding to the SSM

        dat.state1 = ''  # maternal
        dat.state2 = ''  # paternal

        # do this until we encounter the SSM node,
        # i.e., along the path from root to the SSM node
        visited_cnv = False
        for node in ancestors:

            if node != mr_cnv[0].node and visited_cnv == False:  # until CNV is encountered
                dat.state1 += str(node.id) + ',' + str(2) + ',' + str(0) + ';'
            else:
                visited_cnv = True
                dat.state1 += str(node.id) + ',' + str(
                    mr_cnv[1] + mr_cnv[2]) + ',' + str(0) + ';'
            dat.state2 = dat.state1

        # do this after the SSM node, i.e, for all nodes in the subtree below the SSM node
        # [node_id, nr, nv] format
        def descend(n, d):
            if n == mr_cnv[0].node:
                # maternal
                d.state1 += str(n.id) + ',' + str(mr_cnv[1]) + ',' + str(
                    mr_cnv[2]) + ';'
                # paternal
                d.state2 += str(n.id) + ',' + str(mr_cnv[2]) + ',' + str(
                    mr_cnv[1]) + ';'
            else:
                d.state1 += str(n.id) + ',' + str(mr_cnv[1] + mr_cnv[2] -
                                                  1) + ',' + str(1) + ';'
                d.state2 = d.state1
            for child in n.children():
                descend(child, d)

        # traverse the tree below the ssm node
        for child in node.children():
            descend(child, dat)

        fh.write(
            str(dat.id[1:]) + '\t' + dat.state1.strip(';') + '\t' +
            dat.state2.strip(';'))
        fh.write('\n')

    fh.flush()
    fh.close()

#!/usr/bin/env python2
import cPickle

from numpy import array,around,mean,std
from numpy.random import *

from phySCNAClonal.model.tssb import TSSB
from phySCNAClonal.model.util2 import TreeReader, TreeWriter, remove_empty_nodes

from ete2 import *
import heapq
from subprocess import call
import argparse
import sys


def compute_lineages(archiveFn, treeNumPrint, stripePoolFilePath):
    stripeL, baseline = load_data(stripePoolFilePath)

    stripeLNum = len(stripeL)

    treeReader = TreeReader(archiveFn)
    treeNum = treeReader.num_trees()  # number of MCMC samples

    G = dict()
    stripeNameIdxD = dict([(stripe.name, str(i)) for i, stripe in enumerate(stripeL)])

    for idx, llh, tssb in treeReader.load_trees_and_metadata():
        wts, nodes = tssb.get_mixture()
        nodeDataLD = dict([(id(node), '') for i, node in enumerate(nodes)])
        dataIdxSL = []

        # process descendants

        def descend(root):
            data = root.get_data()
            ndata = len(data)
            dataIdxS = ''
            if ndata > 0:
                if root.parent() is not None:
                    dataIdxS = nodeDataLD[id(root.parent())]  # ??
                for stripe in data:
                    # idx_idx_..
                    dataIdxS += stripeNameIdxD[stripe.name] + '_'
                nodeDataLD[id(root)] = dataIdxS
                dataIdxSL.append(_sort(dataIdxS))
            for child in root.children():
                descend(child)

        descend(tssb.root['node'])

        dataIdxSS = sort_and_merge(dataIdxSL)
        if dataIdxSS in G:
            G[dataIdxSS].append(idx)
        else:
            G[dataIdxSS] = [idx]

    postTrees = []
    for ptree in G.keys():
        tssbIdxL = G[ptree]
        prob = round(len(tssbIdxL) * 1.0 / treeNum, 4)
        # print 'posterior probability: ' + repr(prob)
        idx = tssbIdxL[0]
        # Note that trees aren't ordered by likelihood -- only posterior
        # probabaility.
        heapq.heappush(postTrees, (-1.0 * prob, tssbIdxL))

    # print the trees in latex format
    try:
        os.mkdir('posterior_trees')
    except OSError as e:
        if e.errno == 17:  # Directory exists
            pass
        else:
            raise e
    fIdx = 0

    if treeNumPrint is None:
        treeNumPrint = len(postTrees)

    while len(postTrees) > 0 and fIdx < treeNumPrint:
        score, idxL = heapq.heappop(postTrees)
        score = -score

        # aggregate frequencies
        freqs = dict()
        for idx in idxL:
            tssb = treeReader.load_tree(idx)
            remove_empty_nodes(tssb.root, None)

            def descend(root):
                for child in root.children():
                    descend(child)
                names = ''
                for dat in root.get_data():
                    names += dat.name + ';'
                names = names.strip(';')
                if names in freqs:
                    freqs[names].append(root.param)
                else:
                    if len(names) > 0 or root.parent() is None:
                        freqs[names] = []
                        freqs[names].append(root.param)

            descend(tssb.root['node'])

        texFileName = 'posterior_trees/tree_%s_%s.tex' % (fIdx, score)
        print_best_tree(treeReader.load_tree(idxL[0]), texFileName, score, freqs)

        # Call pdflatex. To permit it to find standalone.* files,
        # change into PhyloWGS directory to run the command, then
        # change back to previous directory.
        script_dir = os.path.dirname(os.path.realpath(__file__))
        old_wd = os.getcwd()
        os.chdir(script_dir)
        try:
            call([
                'pdflatex', '-interaction=nonstopmode',
                '-output-directory=%s/posterior_trees/' % old_wd,
                '%s/%s' % (old_wd, texFileName)
            ])
        except OSError:  # pdflatex not available, do not die
            print >> sys.stderr, 'pdflatex not available'
        os.chdir(old_wd)

        fIdx += 1
    treeReader.close()


def _sort(string):
    string = sorted(string.split('_'))
    sstr = ''
    for s in string:
        sstr += s + '_'
    return sstr


def sort_and_merge(dataIdxSL):
    sglist = sorted(dataIdxSL)
    sstr = ''
    for s in sglist:
        sstr += s + ';'
    return sstr


### printing stuff #################


def print_best_tree(tssb, fout, score, freqs):
    remove_empty_nodes(tssb.root, None)  # removes empty leaves

    #wts, nodes = tssb.get_mixture()
    #w = dict([(n[1], n[0]) for n in zip(wts,nodes)])

    print_tree_latex(tssb, fout, score, freqs)


################ LATEX PRINTING ######################
global count

# writes code for tree
# root: root of the tree
# treeFile: string with latex code


def write_tree(root, treeFile):
    global count
    count += 1
    treeFile += 'node {{{0}}}'.format(count)
    for child in root.children():
        treeFile += 'child {'
        treeFile = write_tree(child, treeFile)
        treeFile += '}'
    return treeFile


# writes code for index
# root: root of the tree
# treeFile: string with latex code


def print_index(root, treeFile, freqs):
    global count
    count += 1
    treeFile += '{0} & '.format(count)
    treeFile += '%s & ' % (len(root.get_data()))

    # print param
    names = ''
    for dat in root.get_data():
        names += dat.name + ';'
    names = names.strip(';')
    freq = array(freqs[names])
    treeFile += '{0}'.format(
        str(around(mean(freq[:, i]), 3)) + ' $\pm$ ' + str(
            around(std(freq[:, i]), 3)))
    treeFile = treeFile[:-2]
    treeFile += '\\\\\n'

    for child in root.children():
        treeFile = print_index(child, treeFile, freqs)
    return treeFile


# writes the latex code
# tssb: tssb structure of the tree
# fout: output file for latex


def print_tree_latex(tssb, fout, score, freqs):
    global count
    #remove_empty_nodes(tssb.root, None)

    fout = open(fout, 'w')
    count = -1
    # treeFile='\documentclass{article}\n'
    treeFile = '\documentclass{standalone}\n'
    treeFile += '\usepackage{tikz}\n'
    treeFile += '\usepackage{multicol}\n'
    treeFile += '\usetikzlibrary{fit,positioning}\n'
    treeFile += '\\begin{document}\n'
    treeFile += '\\begin{tikzpicture}\n'
    treeFile += '\\node (a) at (0,0){\n'
    treeFile += '\\begin{tikzpicture}\n'
    treeFile += '[grow=east, ->, level distance=20mm,\
	every node/.style={circle, minimum size = 8mm, thick, draw =black,inner sep=2mm},\
	every label/.append style={shape=rectangle, yshift=-1mm},\
	level 2/.style={sibling distance=50mm},\
	level 3/.style={sibling distance=20mm},\
	level 4/.style={sibling distance=20mm},\
	every edge/.style={-latex, thick}]\n'

    treeFile += '\n\\'
    treeFile = write_tree(tssb.root['node'], treeFile)
    treeFile += ';\n'
    treeFile += '\\end{tikzpicture}\n'
    treeFile += '};\n'
    count = -1
    treeFile += '\\node (b) at (a.south)[anchor=north,yshift=-.5cm]{\n'
    treeFile += '\\begin{tikzpicture}\n'
    treeFile += '\\node (table){\n'
    treeFile += '\\begin{tabular}{|c|l|'
    for i in range(len(tssb.root['node'].param)):
        treeFile += 'l|'
    treeFile += '}\n'
    treeFile += '\\hline\n'
    treeFile += 'Node & \multicolumn{{1}}{{|c|}}{{Stripes}} & \multicolumn{{{0}}}{{|c|}}{{Clonal frequencies}}\\\\\n'.format(
        len(tssb.root['node'].param))
    treeFile += '\\hline\n'
    treeFile = print_index(tssb.root['node'], treeFile, freqs)
    treeFile += '\\hline\n'
    treeFile += '\\end{tabular}\n'
    treeFile += '};\n'
    treeFile += '\\end{tikzpicture}\n'
    treeFile += '};\n'
    treeFile += '\\node at (b.south) [anchor=north,yshift=-.5cm]{Posterior probability: ' + str(
        score) + '};\n'
    treeFile += '\\end{tikzpicture}\n'
    treeFile += '\end{document}\n'
    fout.write(treeFile)
    fout.close()


def process(args):
    compute_lineages(TreeWriter.defaultArchiveFn, args.treeNumPrint,
                     args.stripePoolFilePath)


def main():
    process()


if __name__ == "__main__":
    main()

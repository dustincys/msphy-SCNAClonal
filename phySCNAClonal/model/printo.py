import sys
from subprocess import call

from ete2 import *
from numpy import *
from numpy.random import *

from phySCNAClonal.model.tssb import *
from phySCNAClonal.model.util import *
# for tssb show tree, comment
# from phySCNAClonal.model.util2 import remove_empty_nodes, TreeReader
from collections import Counter

import os


def show_tree_structure2(tssb,
                        texFolder,
                        pdfFolder,
                        filePrefix,
                        currentTimeTag,
                        allTimeTags,
                        toCompile=False,
                        toShow=False):

    texFileFullPath = os.path.join(texFolder, filePrefix + ".tex")
    print_tree_latex2(tssb, texFileFullPath, currentTimeTag, allTimeTags)


    if toCompile:
        try:
            callList = ["/usr/local/texlive/2017/bin/x86_64-linux/pdflatex",
                        "-interaction=nonstopmode",
                        "-output-directory={}".format(pdfFolder),
                        texFileFullPath]
            call(" ".join(callList), shell=True)
            call("rm {0}/{1}.aux {0}/{1}.log".format(pdfFolder, filePrefix),
                 shell=True)
        except OSError as oser:  # pdflatex not available, do not die
            print >> sys.stderr, 'pdflatex not available'
        if toShow:
            try:
                callList = ["/usr/bin/okular",
                            os.path.basename(texFileFullPath)]
                call(" ".join(callList), shell=True)
            except OSError:  # pdflatex not available, do not die
                print >> sys.stderr, 'okular not available'

def show_tree_structure3(tssb,
                        texFolder,
                        pdfFolder,
                        filePrefix,
                        toCompile=False,
                        toShow=False):

    texFileFullPath = os.path.join(texFolder, filePrefix + ".tex")
    print_tree_latex3(tssb, texFileFullPath)


    if toCompile:
        try:
            callList = ["/usr/local/texlive/2017/bin/x86_64-linux/pdflatex",
                        "-interaction=nonstopmode",
                        "-output-directory={}".format(pdfFolder),
                        texFileFullPath]
            call(" ".join(callList), shell=True)
            call("rm {0}/{1}.aux {0}/{1}.log".format(pdfFolder, filePrefix),
                 shell=True)
        except OSError as oser:  # pdflatex not available, do not die
            print >> sys.stderr, 'pdflatex not available'
        if toShow:
            try:
                callList = ["/usr/bin/okular",
                            os.path.basename(texFileFullPath)]
                call(" ".join(callList), shell=True)
            except OSError:  # pdflatex not available, do not die
                print >> sys.stderr, 'okular not available'

ctr=0
def print_top_trees(treeArchive, fout, k=5):
    global ctr;
    foutF = open(fout,'w')
    treeReader = TreeReader(treeArchive)
    i = 0
    for idx, (tidx, llh, tree, dp) in enumerate(
        treeReader.load_trees_and_metadata(k)):
            ctr=0
            remove_empty_nodes(tree.root, None)
            # print top K trees in ascii
            if i < 5 :
                print_best_tree_pdf(tree, fout+"temp{}.tex".format(i))
            print_best_tree(tree, foutF)
            i = i + 1

    treeReader.close()
    foutF.close()

def print_best_tree(tssb,fout):
    nodes = tssb.get_nodes()
    nnodes = sum( [ 1 for node in nodes if len(node.get_data()) ] )

    #fout=open(fout,'w')
    t = Tree();t.name='0'
    fout.write('id, \t phi, \t nChildren, \t nGenes, \t stripes \n')
    print_node2(tssb.root,None,t,fout)
    fout.write('\n\n')
    fout.write(t.get_ascii(show_internal=True))
    fout.write('\n\n')
    fout.write('Number of non-empty nodes in the tree: ' +repr(nnodes))
    fout.write('\n\n\n')
    #fout.close()

def print_node2(node, parent,tree,fout):
    global ctr;
    num_data = node['node'].num_data()
    node_name  = ctr ; ctr+=1;

    stripes = node['node'].get_data()
    stripeId = ''
    if len(stripes)>0:

        stripeId = stripes[0].sid #name
        for spIdx in arange(1,len(stripes)):
            stripeId = stripeId + '; ' + stripes[spIdx].sid #name

    out_str = str(node_name) + ',\t' + str(around(node['node'].param,3)) +\
        ',\t' + str(len(node['children'])) + ',\t' + str(len(stripes)) + ',\t' +\
        stripeId +  '\n'

    fout.write(out_str)
    for child in node['children']:
        name_string = str(ctr)#+'('+str(len(child['node'].get_data()))+')'
        print_node2(child, node_name,tree.add_child(name=name_string),fout)

### printing stuff #################
def print_best_tree_pdf(tssb,fout,score=0):
    #wts, nodes = tssb.get_mixture()
    #w = dict([(n[1], n[0]) for n in zip(wts,nodes)])
    print_tree_latex(tssb,fout,score)


################ LATEX PRINTING ######################
global count
# writes code for tree
# root: root of the tree
# tree_file: string with latex code
def write_tree(root):
    global count
    count+=1

    reStr = "{0}".format(count)

    if 'tag' not in root.keys():
        print "not tag"
    if root['tag']:
        reStr += ",draw=red"
    else:
        reStr += ",draw=black"

    for child in root['children']:
        reStr += write_tree(child)

    return "[{}]".format(reStr)


def write_tree2(root, timeTag, timeTags):
    global count
    count+=1

    def colorize(num, t, timeTag):
        if t == timeTag and num > 0:
            return "\\textcolor{red}{" + str(num) + "}"
        else:
            return str(num)

    reStr = ""
    tempTimeTags = Counter([int(item.tag) for item in root['node'].get_data()])
    if len(tempTimeTags) > 0:
        for t in timeTags:
            reStr += "{0}$|$".format(colorize(tempTimeTags[t], t, timeTag))
        reStr = reStr.strip("$|$")
    else:
        pass

    if 'tag' not in root.keys():
        print "not tag"
    if root['tag']:
        reStr += ",draw=black,fill=black!20!white"
    else:
        reStr += ",draw=black"

    for child in root['children']:
        reStr += write_tree2(child, timeTag, timeTags)

    return "[{}]".format(reStr)


def write_tree3(root):
    # for debugging crossing & summing rule
    global count
    count+=1

    reStr = "[{"
    for d in  root['node'].data:
        reStr += "{0},".format(d)
    reStr = reStr.strip(",")
    reStr += "}"

    for child in root['children']:
        reStr += write_tree3(child)

    return reStr + "]"


# writes code for index
# root: root of the tree
# tree_file: string with latex code
def print_index(root, tree_file):
    global count
    count+=1
    tree_file+='{0} & '.format(count)

    stripes=''
    for sp in root['node'].get_data():
        # stripes+='{0}, '.format(sp.name)
        stripes += "{0}, ".format(sp.tag)
    stripes += str(root['tag'])
    stripes = stripes.replace("_","-")

    tree_file+=stripes.strip().strip(',')

    if root['node'].get_data()==[]:
        tree_file+='-- '

    tree_file+=' & '
    tree_file+=str(around(root['node'].param,3))
    tree_file+='\\\\\n'

    for child in root['children']:
        tree_file=print_index(child, tree_file)
    return tree_file

# writes the latex code
# tssb: tssb structure of the tree
# fout: output file for latex
def print_tree_latex(tssb,fout,score):
    global count
    fout = open(fout,'w')
    count=-1
    #tree_file='\documentclass{article}\n'
    tree_file='\\documentclass{standalone}\n'
    tree_file+='\\usepackage{forest}\n'
    tree_file+='\\usepackage{tikz}\n'
    tree_file+='\\usepackage{multicol}\n'
    tree_file+='\\usetikzlibrary{fit,positioning}\n'
    tree_file+='\\begin{document}\n'
    tree_file+='\\begin{tikzpicture}\n'
    tree_file+='\\node (a) at (0,0){\n'
    tree_file+='\\begin{forest}\n'
    tree_file+='before typesetting nodes={for descendants={edge=->}}\n'
    tree_file+=write_tree(tssb.root)
    tree_file+='\\end{forest}\n'
    tree_file+='};\n'
    count=-1
    tree_file+='\\node (b) at (a.south)[anchor=north,yshift=-.5cm]{\n'
    tree_file+='\\begin{tikzpicture}\n'
    tree_file+='\\node (table){\n'
    tree_file+='\\begin{tabular}{|c|p{5cm}|p{5cm}|'
    tree_file+='l|'
    tree_file+='}\n'
    tree_file+='\\hline\n'
    tree_file+='Node & \multicolumn{{1}}{{|c|}}{{Mutations}} & \multicolumn{{1}}{{|c|}}{{Frequencies}} \\\\\n'.format("1")
    tree_file+='\\hline\n'
    tree_file=print_index(tssb.root, tree_file)
    tree_file+='\\hline\n'
    tree_file+='\\end{tabular}\n'
    tree_file+='};\n'
    tree_file+='\\end{tikzpicture}\n'
    tree_file+='};\n'
    #tree_file+='\\node at (b.south) [anchor=north,yshift=-.5cm]{Posterior probability: ' + str(score) + '};\n'
    tree_file+='\\end{tikzpicture}\n'
    tree_file+='\end{document}\n'
    fout.write(tree_file)
    fout.close()

def print_tree_latex2(tssb,fout,timeTag,timeTags):
    global count
    fout = open(fout,'w')
    count=-1
    #tree_file='\documentclass{article}\n'
    tree_file='\\documentclass{standalone}\n'
    tree_file+='\\usepackage{forest}\n'
    tree_file+='\\usepackage{xcolor}\n'
    tree_file+='\\usepackage{tikz}\n'
    tree_file+='\\usepackage{multicol}\n'
    tree_file+='\\usetikzlibrary{fit,positioning}\n'
    tree_file+='\\begin{document}\n'
    tree_file+='\\begin{forest}\n'
    tree_file+='before typesetting nodes={for descendants={edge=->}}\n'
    tree_file+=write_tree2(tssb.root, timeTag, timeTags)
    tree_file+='\\end{forest}\n'
    tree_file+='\end{document}\n'
    fout.write(tree_file)
    fout.close()


def print_tree_latex3(tssb,fout):
    global count
    fout = open(fout,'w')
    count=-1
    #tree_file='\documentclass{article}\n'
    tree_file='\\documentclass{standalone}\n'
    tree_file+='\\usepackage{forest}\n'
    tree_file+='\\usepackage{xcolor}\n'
    tree_file+='\\usepackage{tikz}\n'
    tree_file+='\\usepackage{multicol}\n'
    tree_file+='\\usetikzlibrary{fit,positioning}\n'
    tree_file+='\\begin{document}\n'
    tree_file+='\\begin{forest}\n'
    tree_file+='before typesetting nodes={for descendants={edge=->}}\n'
    tree_file+=write_tree3(tssb.root)
    tree_file+='\\end{forest}\n'
    tree_file+='\end{document}\n'
    fout.write(tree_file)
    fout.close()

import sys
from subprocess import call

from ete2 import *
from numpy import *
from numpy.random import *

from tssb import *
from util import *
from util2 import remove_empty_nodes, TreeReader

ctr=0
def print_top_trees(treeArchive, fout, k=5):
    global ctr;
    foutF = open(fout,'w')
    treeReader = TreeReader(treeArchive)
    i = 0
    for idx, (tidx, llh, tree) in enumerate(treeReader.load_trees_and_metadata(k)):
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
def write_tree(root, tree_file):
    global count
    count+=1
    tree_file+='node {{{0}}}'.format(count)
    for child in root.children():
        tree_file+='child {'
        tree_file=write_tree(child, tree_file)
        tree_file+='}'
    return tree_file

# writes code for index
# root: root of the tree
# tree_file: string with latex code
def print_index(root, tree_file):
    global count
    count+=1
    tree_file+='{0} & '.format(count)

    stripes=''
    for sp in root.get_data():
        # stripes+='{0}, '.format(sp.name)
        stripes += "{0}, ".format(sp.tag)
    stripes = stripes.replace("_","-")

    tree_file+=stripes.strip().strip(',')

    if root.get_data()==[]:
        tree_file+='-- '

    tree_file+=' & '
    tree_file+=str(around(root.param,3))
    tree_file+='\\\\\n'

    for child in root.children():
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
    tree_file+='\\usepackage{tikz}\n'
    tree_file+='\\usepackage{multicol}\n'
    tree_file+='\\usetikzlibrary{fit,positioning}\n'
    tree_file+='\\begin{document}\n'
    tree_file+='\\begin{tikzpicture}\n'
    tree_file+='\\node (a) at (0,0){\n'
    tree_file+='\\begin{tikzpicture}\n'
    tree_file+='[grow=east, ->, level distance=20mm,\
        every node/.style={circle, minimum size = 8mm, thick, draw =black,inner sep=2mm},\
        every label/.append style={shape=rectangle, yshift=-1mm},\
        level 2/.style={sibling distance=50mm},\
        level 3/.style={sibling distance=20mm},\
        level 4/.style={sibling distance=20mm},\
        every edge/.style={-latex, thick}]\n'
    tree_file+='\n\\'
    tree_file=write_tree(tssb.root['node'], tree_file)
    tree_file+=';\n'
    tree_file+='\\end{tikzpicture}\n'
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
    tree_file=print_index(tssb.root['node'], tree_file)
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


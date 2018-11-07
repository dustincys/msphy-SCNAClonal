#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: figure_generator.py
#          Desc: generate figure for each tree
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-09-14 19:58:20
#       History:
# =============================================================================
'''

import os
from collections import Counter
from subprocess import call

import numpy as np

from phySCNAClonal.model.util2 import remove_empty_nodes
from phySCNAClonal.preprocess.utils import AnswerIndex


class FigureGenerator(object):

    """Docstring for FigureGenerator. """

    def __init__(self):
        """Initialize input object"""
        self._get_output_prefixes()
        pass

    def __init__(self, outputPrefix = "figureGenerator", answerFilePath = None):
        self._outputPrefix = outputPrefix
        self._answerFilePath = answerFilePath
        self._get_output_prefixes()

    def set_output_prefix(self, outputPrefix):
        self._outputPrefix = outputPrefix
        self._get_output_prefixes()

    def set_answer_file_path(self, answerFilePath = None):
        self._answerFilePath = answerFilePath

    def _get_output_prefixes(self):
        self._treeOutputPrefix = self._outputPrefix + "_tree"
        self._rOutputPrefix = self._outputPrefix + "_r"
        self._copyNumberOutputPrefix = self._outputPrefix + "_copyNumber"
        self._phiOutputPrefix = self._outputPrefix + "_phi"

    def draw_tree(self, tssb, timeTagsL, score=None):
        texFilePath = self._treeOutputPrefix + ".tex"
        remove_empty_nodes(tssb.root, None)  # removes empty leaves
        self._print_tree_latex(tssb, timeTagsL, texFilePath, score)
        self._draw_pdf(texFilePath)

    def _draw_pdf(self, texFilePath):
        pdfFolder = os.path.dirname(texFilePath)
        try:
            callList = ["/usr/local/texlive/2017/bin/x86_64-linux/pdflatex",
                        "-interaction=nonstopmode",
                        "-output-directory={}".format(pdfFolder),
                        texFilePath]
            call(" ".join(callList), shell=True)
            call("rm {0}.aux {0}.log".format(self._treeOutputPrefix),
                shell=True)
        except OSError as oser:  # pdflatex not available, do not die
            print >> sys.stderr, 'pdflatex not available'

    def draw_data_param(self, SCNAPool, isStripe = False):
        rTablePath = self._rOutputPrefix + ".txt"
        copyNumberPath = self._copyNumberOutputPrefix + ".png"

        self._output_r_table(SCNAPool,
                             rTablePath,
                             self._answerFilePath,
                             isStripe)

        self._r_draw_CopyNumber(rTablePath, copyNumberPath)
        self._r_draw_phi(rTablePath, self._phiOutputPrefix)

    def draw_all_result(self, SCNAPool, tssb, isStripe=False, score=None):
        timeTagsL = []
        if isStripe:
            timeTagsL = sorted(list(set([int(item.tag) for item in
                                         SCNAPool.stripes])))
        else:
            timeTagsL = sorted(list(set([int(item.tag) for item in
                                         SCNAPool.segments])))

        self.draw_tree(tssb, timeTagsL, score)
        self.draw_data_param(SCNAPool, isStripe)

    def _print_tree_latex(self, tssb, timeTagsL, outputFilePath, score):
        count = {'nodeId': -1}

        def write_tree(root, timeTagsL, count):
            count['nodeId'] += 1

            reStr = "\\textcolor{blue}{" + str(count['nodeId']) +  "}("
            tempTimeTags = Counter([int(item.tag) for item in root['node'].get_data()])
            if len(tempTimeTags) > 0:
                for t in timeTagsL:
                    reStr += "{0}$|$".format(tempTimeTags[t])
                reStr = reStr.strip("$|$")
            else:
                pass
            reStr += ")"
            reStr += ",draw=black"

            for child in root['children']:
                reStr += write_tree(child, timeTagsL, count)

            return "[{}]".format(reStr)

        def print_index(root, treeFile, count):
            count['nodeId'] += 1

            treeFile += "\\textcolor{blue}{" +str(count['nodeId']) +  "} & "

            SCNAStr=''
            for s in root['node'].get_data():
                SCNAStr += "{0}$\|${1}, ".format(s.tag, s.name)

            #  TODO:  <14-09-18, Chu Yanshuo> #
            SCNAStr = SCNAStr.replace("_","\\_")

            treeFile += SCNAStr.strip().strip(',')

            if root['node'].get_data()==[]:
                treeFile+='-- '

            treeFile+=' & '
            treeFile+=str(np.around(root['node'].param,3))
            treeFile+='\\\\\n'

            for child in root['children']:
                treeFile=print_index(child, treeFile, count)
            return treeFile

        treeFile = '\documentclass{standalone}\n'
        treeFile +='\\usepackage{forest}\n'
        treeFile +='\\usepackage{xcolor}\n'
        treeFile += '\usepackage{tikz}\n'
        treeFile += '\usepackage{multicol}\n'
        treeFile += '\usetikzlibrary{fit,positioning}\n'
        treeFile += '\\begin{document}\n'
        treeFile += '\\begin{tikzpicture}\n'
        treeFile += '\\node (a) at (0,0){\n'
        treeFile += '\\begin{forest}\n'
        treeFile += 'before typesetting nodes={for descendants={edge=->}}\n'
        treeFile += write_tree(tssb.root, timeTagsL, count)
        treeFile += '\\end{forest}\n'
        treeFile += '};\n'
        count['nodeId'] = -1
        treeFile += '\\node (b) at (a.south)[anchor=north,yshift=-.5cm]{\n'
        treeFile += '\\begin{tikzpicture}\n'
        treeFile += '\\node (table){\n'
        treeFile += '\\begin{tabular}{|c|l|l|}\n'
        treeFile += '\\hline\n'
        treeFile += 'Node & SCNAs (stage$\|$ID) & Clonal frequencies\\\\\n'
        treeFile += '\\hline\n'
        treeFile = print_index(tssb.root, treeFile, count)
        treeFile += '\\hline\n'
        treeFile += '\\end{tabular}\n'
        treeFile += '};\n'
        treeFile += '\\end{tikzpicture}\n'
        treeFile += '};\n'
        treeFile += '\\node at (b.south) [anchor=north,yshift=-.5cm]{Posterior probability: ' + str(
            score) + '};\n'
        treeFile += '\\end{tikzpicture}\n'
        treeFile += '\end{document}\n'

        with open(outputFilePath, 'w') as outputFile:
            outputFile.write(treeFile)

    def _output_r_table(self, SCNAPool, outputFilePath, answerFilePath=None,
                        isStripe=False):
        segPool = None
        if isStripe:
            segPool = SCNAPool.segPool
        else:
            segPool = SCNAPool

        if answerFilePath is None:
            with open(outputFilePath, 'w') as outputFile:
                outputFile.write(segPool.segments[0].toName() +"\n")
                for seg in segPool.segments:
                    outputFile.write(seg.toString() +"\n")
        else:
            ansIdx = AnswerIndex(answerFilePath)

            with open(outputFilePath, 'w') as outputFile:
                outputFile.write(segPool.segments[0].toName() +
                                "\tcopyNumberAnswer\tgenotypeAnswer\tphiAnswer\n")
                for seg in segPool.segments:
                    value = ansIdx.getValue(seg.chromName,
                                            seg.start,
                                            seg.end)
                    value = value.strip()
                    listValue = value.split('\t')

                    outputFile.write(seg.toString() +"\t{0}\t{1}\t{2}\n".format(
                        listValue[-3], listValue[-2], listValue[-1]))
                    pass
                pass
            pass
        pass

    def _r_draw_CopyNumber(self, rTableFilePath, outputFilePath):
        """
        call r script to generate copy number figure
        """
        rScriptFilePath = os.path.join(os.path.dirname(
            os.path.realpath(__file__)), 'drawCopyNumber.R')
        try:
            callList = ["/usr/bin/Rscript", rScriptFilePath,
                        rTableFilePath, outputFilePath]
            call(" ".join(callList), shell=True)
        except OSError as oser:
            print >> oser
            print >> sys.stderr, 'Rscript not available'

    def _r_draw_phi(self, rTableFilePath, outputFilePath):
        """
        call r script to generate phi figure
        """
        rScriptFilePath = os.path.join(os.path.dirname(
            os.path.realpath(__file__)), 'drawPhi.R')
        try:
            callList = ["/usr/bin/Rscript", rScriptFilePath,
                        rTableFilePath, outputFilePath]
            call(" ".join(callList), shell=True)
        except OSError as oser:
            print >> oser
            print >> sys.stderr, 'Rscript not available'


#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from phySCNAClonal.postprocess.pwgsresults.figure_generator import \
    FigureGenerator
from phySCNAClonal.postprocess.pwgsresults.result_generator import \
    ResultGenerator


def output_partial_data(outFilePath, pd, treeNum):
    outFileP = open(outFilePath+".linkProb.txt", 'w')
    outFileT = open(outFilePath+".tree.txt", "w")

    for fromNode in pd.keys():
        totalLink = sum([pd[fromNode][toNode] for toNode in pd[fromNode]])
        for toNode in pd[fromNode]:
            outFileP.write("{0}\t{1}\t{2}\n".format(
                fromNode, toNode, pd[fromNode][toNode]*10.0/totalLink))
            outFileT.write("{0}\t{1}\t{2}\n".format(
                fromNode, toNode, pd[fromNode][toNode]*10.0/treeNum))

    outFileT.close()
    outFileP.close()

def output_result_distribution_data(outFilePath, summaries, bestIdx):
    outFileDist = open(outFilePath+".result_distribution.txt", 'w')
    outFileBest = open(outFilePath+".bestTree.txt", 'w')


    headLine = True
    for treeIdx in summaries.keys():
        SCNAL = summaries[treeIdx]['SCNAL']
        for scna in SCNAL:
            if headLine:
                outFileDist.write(scna.toName()+"\n")
                headLine = False
            outFileDist.write(scna.toString()+"\n")

    headLine = True
    SCNAL = summaries[bestIdx]['SCNAL']
    for scna in SCNAL:
        if headLine:
            outFileBest.write(scna.toName()+"\n")
            headLine = False
        outFileBest.write(scna.toString()+"\n")

    outFileDist.close()
    outFileBest.close()


def process(args):
    # 此处应该生成关于目标的copyNumber, phi等信息
    treeNum, summaries, allMutAss, params, partialDict, isStripe, bestIdx =\
        ResultGenerator().generate(args.treeFile, args.SCNAPoolFile)

    outFilePath = args.outputFolder + "/partialData"
    output_partial_data(outFilePath, partialDict, treeNum)

    outFilePath = args.outputFolder + "/resultData"
    output_result_distribution_data(outFilePath, summaries, bestIdx)

    # 此处生成图形
    fg = FigureGenerator()
    fg.set_answer_file_path(args.answerFilePath)

    # for treeIdx in summaries.keys():
        # figPrefix = args.outputFolder + "/tree{}".format(treeIdx)
        # fg.set_output_prefix(figPrefix)
        # fg.draw_all_result(
            # summaries[treeIdx]['SCNAL'],
            # summaries[treeIdx]['tree'],
            # isStripe,
            # summaries[treeIdx]['llh'])

    figPrefix = args.outputFolder + "/bestTree{}".format(bestIdx)
    fg.set_output_prefix(figPrefix)
    fg.draw_all_result(
        summaries[bestIdx]['SCNAL'],
        summaries[bestIdx]['tree'],
        isStripe,
        summaries[bestIdx]['llh'])
    # 堆排序重写

    # 此处将生成关于目标的copyNumber, phi等信息输出到文件
    # writer = JsonWriter(args.datasetName)
    # writer.write_summaries(summaries, params, args.treeSummaryOutput)
    # writer.write_mutlist(mutlist, args.mutlistOutput)
    # writer.write_mutass(allMutAss, args.mutassOutput)

def main():
    process(args)

if __name__ == '__main__':
    main()

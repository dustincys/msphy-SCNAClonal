#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pickle as pkl

from phySCNAClonal.postprocess.pwgsresults.figure_generator import \
    FigureGenerator
from phySCNAClonal.postprocess.pwgsresults.result_generator import \
    ResultGenerator


def output_partial_data(outFilePath, pd, SCNAPool, isSegment=True):
    with open(outFilePath, 'w') as outFile:
        outFile.write("fromNode\ttoNode\tweight\tfromNodeName\ttoNodeName\n")
        for fromNode in pd.keys():
            totalLink = sum([pd[fromNode][toNode] for toNode in pd[fromNode]])
            for toNode in pd[fromNode]:
                if isSegment:
                    outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                        fromNode, toNode, pd[fromNode][toNode]*1.0/totalLink,
                        SCNAPool.segments[fromNode].name,
                        SCNAPool.segments[toNode].name))
                else:
                    outFile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                        fromNode, toNode, pd[fromNode][toNode]*1.0/totalLink,
                        SCNAPool.stripes[fromNode].name,
                        SCNAPool.stripes[toNode].name))

def process(args):
    # 此处应该生成关于目标的copyNumber, phi等信息
    summaries, allMutAss, params, partialDict, isStripe =\
        ResultGenerator().generate(args.treeFile, args.SCNAPoolFile)

    outFilePath = args.outputFolder + "/partialData.txt"
    pklFile = open(args.SCNAPoolFile, 'rb')
    SCNAPool = pkl.load(pklFile)
    pklFile.close()
    if hasattr(SCNAPool, 'segments'):
        output_partial_data(outFilePath, partialDict, SCNAPool, True)
    else:
        output_partial_data(outFilePath, partialDict, SCNAPool, False)

    # 此处生成图形
    fg = FigureGenerator()
    fg.set_answer_file_path(args.answerFilePath)

    for treeIdx in summaries.keys():
        figPrefix = args.outputFolder + "/tree{}".format(treeIdx)
        fg.set_output_prefix(figPrefix)
        fg.draw_all_result(
            summaries[treeIdx]['SCNAPool'],
            summaries[treeIdx]['tree'],
            isStripe,
            summaries[treeIdx]['llh'])

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

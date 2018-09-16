#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from phySCNAClonal.postprocess.pwgsresults.figure_generator import \
    FigureGenerator
from phySCNAClonal.postprocess.pwgsresults.result_generator import \
    ResultGenerator


def process(args):
    # 此处应该生成关于目标的copyNumber, phi等信息
    summaries, allMutAss, params, isStripe =\
        ResultGenerator().generate(args.treeFile, args.SCNAPoolFile)

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

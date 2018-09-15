#!/usr/bin/env python2
import argparse
from phySCNAClonal.postprocess.pwgsresults.result_generator import ResultGenerator
from phySCNAClonal.postprocess.pwgsresults.result_munger import ResultMunger
from phySCNAClonal.postprocess.pwgsresults.json_writer import JsonWriter


def process(args):
    # 此处应该生成关于目标的copyNumber, phi等信息
    summaries, allMutAss, params, isStripe =\
        ResultGenerator().generate(args.treeFile, args.segPoolFile)

    # 此处生成图形

    # 此处将生成关于目标的copyNumber, phi等信息输出到文件
    writer = JsonWriter(args.datasetName)
    writer.write_summaries(summaries, params, args.treeSummaryOutput)
    writer.write_mutlist(mutlist, args.mutlistOutput)
    writer.write_mutass(mutass, args.mutassOutput)


def main():
    process(args)


if __name__ == '__main__':
    main()

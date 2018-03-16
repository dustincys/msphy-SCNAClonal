#!/usr/bin/env python2
import argparse
from phySCNAClonal.postprocess.pwgsresults.result_generator import ResultGenerator
from phySCNAClonal.postprocess.pwgsresults.result_munger import ResultMunger
from phySCNAClonal.postprocess.pwgsresults.json_writer import JsonWriter


def process(args):
    # 此处应该生成关于目标的copyNumber, phi等信息
    summaries, mutlist, mutass, params = ResultGenerator().generate(
        args.treeFile, args.segPoolFile, args.includeStripeNames)

    munger = ResultMunger(summaries, mutlist, mutass)

    # 此处应该不必去除
    # summaries, mutass = munger.remove_small_nodes(args.min_ssms)
    # 去除多克隆和超级克隆，此处可以使用
    munger.remove_superclones()
    munger.remove_polyclonal_trees()

    # 此处将生成关于目标的copyNumber, phi等信息输出到文件
    writer = JsonWriter(args.datasetName)
    writer.write_summaries(summaries, params, args.treeSummaryOutput)
    writer.write_mutlist(mutlist, args.mutlistOutput)
    writer.write_mutass(mutass, args.mutassOutput)


if __name__ == '__main__':
    main()

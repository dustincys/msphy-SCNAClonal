#!/usr/bin/env python2
import argparse
from pwgsresults.result_generator import ResultGenerator
from pwgsresults.result_munger import ResultMunger
from pwgsresults.json_writer import JsonWriter


def main():
    parser = argparse.ArgumentParser(
        description='Write JSON files describing trees',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '--include-stripe-names',
        dest='includeStripeNames',
        action='store_true',
        help='Include stripe names in output (which may be sensitive data)')
    parser.add_argument(
        '--includeSegmentList',
        dest='includeSegmentList',
        action='store_true',
        help='Include segment list in output (which may be sensitive data)')
    # 此处可以设置为最小Data
    # parser.add_argument(
        # '--min-ssms',
        # dest='min_ssms',
        # type=float,
        # default=0.01,
        # help='Minimum number or percent of SSMs to retain a subclone')
    parser.add_argument('datasetName', help='Name identifying dataset')
    parser.add_argument('treeFile', help='File containing sampled trees')
    parser.add_argument(
        'treeSummaryOutput',
        help='Output file for JSON-formatted tree summaries')
    parser.add_argument(
        'mutlistOutput',
        help='Output file for JSON-formatted list of mutations')
    parser.add_argument(
        'mutassOutput',
        help=
        'Output file for JSON-formatted list of SSMs and CNVs assigned to each subclone'
    )
    args = parser.parse_args()

    summaries, mutlist, mutass, params = ResultGenerator().generate(
        args.treeFile, args.includeStripeNames, args.includeSegmentList)

    munger = ResultMunger(summaries, mutlist, mutass)

# 此处应该不必去除
    # summaries, mutass = munger.remove_small_nodes(args.min_ssms)
    # 去除多克隆和超级克隆，此处可以使用
    munger.remove_superclones()
    munger.remove_polyclonal_trees()

    writer = JsonWriter(args.datasetName)
    writer.write_summaries(summaries, params, args.treeSummaryOutput)
    writer.write_mutlist(mutlist, args.mutlistOutput)
    writer.write_mutass(mutass, args.mutassOutput)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: evolve.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-12-05 16:00:23
#       History:
# =============================================================================
'''
import argparse
import cPickle as pickle
import json
import os
import signal
import sys
import tempfile
import threading
import time
import traceback
from copy import deepcopy
from datetime import datetime

import numpy as np
from gwpy.segments import Segment, SegmentList

from phySCNAClonal.model.params import get_c_fnames, metropolis
from phySCNAClonal.model.printo import print_top_trees, print_tree_latex
from phySCNAClonal.model.stripenode import StripeNode
from phySCNAClonal.model.tssb import TSSB
from phySCNAClonal.model.usegsampler.segsupportive import MultiRangeSampler
from phySCNAClonal.model.util import boundbeta
from phySCNAClonal.model.util2 import (BackupManager, StateManager, TreeWriter,
                                       load_data, map_datum_to_node,
                                       set_node_height,
                                       set_path_from_root_to_node)


# sampleNum: number of MCMC samples
# mhItr: number of metropolis-hasting iterations
# randSeed: random seed (initialization). Set to None to choose random
# seed automatically.

def start_new_run(stateManager,
                  backupManager,
                  safeToExit,
                  runSucceeded,
                  config,
                  stripesFile,
                  stripesTextFile,
                  paramsFile,
                  topKTreesFile,
                  clonalFreqsFile,
                  burninSampleNum,
                  sampleNum,
                  mhItr,
                  mhStd,
                  writeStateEvery,
                  writeBackupsEvery,
                  randSeed,
                  tmpDir,
                  maxCopyNumber):
    state = {}

    try:
        state['rand_seed'] = int(randSeed)
    except TypeError:
        # If randSeed is not provided as command-line arg, it will be None,
        # meaning it will hit this code path.
        #
        # Use random seed in this order:
        #   1. If a seed is given on the command line, use that.
        #   2. Otherwise, if `randomSeed.txt` exists, use the seed stored there.
        #   3. Otherwise, choose a new random seed and write to randomSeed.txt.
        try:
            with open('random_seed.txt') as seedf:
                state['rand_seed'] = int(seedf.read().strip())
        except (TypeError, IOError) as E:
            # Can seed with [0, 2**32).
            state['rand_seed'] = np.random.randint(2**32)

    np.random.seed(state['rand_seed'])
    with open('random_seed.txt', 'w') as seedf:
        seedf.write('%s\n' % state['rand_seed'])

    state['stripes_file'] = stripesFile
    state['stripes_text_file'] = stripesTextFile
    state['tmp_dir'] = tmpDir
    state['top_k_trees_file'] = topKTreesFile
    state['clonal_freqs_file'] = clonalFreqsFile
    state['write_state_every'] = writeStateEvery
    state['write_backups_every'] = writeBackupsEvery

    # 此处载入数据，此处含有baseline
    stripes, baseline = load_data(state['stripes_file'])
    modify_stripes_time_tag(stripes, 3)


    state['baseline'] = baseline
    state['time_tags'] = sorted(list(set([int(item.tag) for item in
                                             stripes])))

    stripeNum = len(stripes)

    if len(stripes) == 0:
        logmsg('No stripes provided. Exiting.', sys.stderr)
        return

    # sample 个数
    # NTPS = len(codes[0].a)  # number of samples / time point
    # gene list
    # state['glist'] = [datum.name for datum in codes if len(datum.name) > 0]

    # stripe list
    state['stripe_list'] = [stripe.name for stripe in stripes]
    # max copy number
    state['max_copy_number'] = maxCopyNumber

    # MCMC settings
    state['burnin_sample_number'] = burninSampleNum
    state['sample_number'] = sampleNum
    state['dp_alpha'] = 25.0
    state['dp_gamma'] = 1.0
    state['alpha_decay'] = 0.25
    state['top_k'] = 5

    # Metropolis-Hastings settings
    state['mh_burnin'] = 0
    state['mh_itr'] = mhItr  # No. of iterations in metropolis-hastings
    state['mh_std'] = mhStd


    state['cd_llh_traces'] = np.zeros((np.power(int(state['sample_number']),
                                                int(len(state['time_tags']))),
                                       1))
    state['burnin_cd_llh_traces'] = np.zeros((
        np.power((state['burnin_sample_number']+ state['sample_number']),
                 len(state['time_tags']))- np.power(state['sample_number'],
                                                    len(state['time_tags'])),
        1))

    state['working_directory'] = os.getcwd()

    root = StripeNode(conc=0.1)

    state['tssb'] = TSSB(
        dpAlpha=state['dp_alpha'],
        dpGamma=state['dp_gamma'],
        alphaDecay=state['alpha_decay'],
        rootNode=root,
        data=stripes,
        baseline=baseline,
        maxCopyNumber=state['max_copy_number'])
    # hack...
    # 初始化把所有数据放到第一个孩子节点
    if 1:
        depth = 0
        state['tssb'].root['sticks'] = np.vstack([
            state['tssb'].root['sticks'],
            boundbeta(1, state['tssb'].dpGamma) if depth != 0 else .999999
        ])
        state['tssb'].root['children'].append({
            'node':
            state['tssb'].root['node'].spawn(),
            'main':
            boundbeta(1.0, (state['tssb'].alphaDecay**
                            (depth + 1)) * state['tssb'].dpAlpha)
            if state['tssb'].minDepth <= (depth + 1) else 0.0,
            'sticks':
            np.empty((0, 1)),
            'children': [],
            'tag': False
        })
        newNode = state['tssb'].root['children'][0]['node']
        for n in range(state['tssb'].dataNum):
            state['tssb'].assignments[n].remove_datum(n)
            newNode.add_datum(n)
            state['tssb'].assignments[n] = newNode

    for stripe in stripes:
        stripe.tssb = state['tssb']

    treeWriter = TreeWriter()
    # 此处放入压缩文件中的额外的配置文件
    # treeWriter.add_extra_file('cnv_logical_physical_mapping.json',
                               # json.dumps(cnv_logical_physical_mapping))

    if paramsFile is not None:
        with open(paramsFile) as F:
            params = json.load(F)
    else:
        params = {}
    treeWriter.add_extra_file('params.json', json.dumps(params))

    # dump to pkl
    stateManager.write_initial_state(state)

    logmsg("Starting MCMC run...")
    state['last_iteration'] = -state['burnin_sample_number'] - 1
    state['total_iteration'] = 0

    # This will overwrite file if it already exists, which is the desired
    # behaviour for a fresh run.
    with open('mcmc_samples.txt', 'w') as mcmcf:
        mcmcf.write('Iteration\tLLH\tTime\n')

    do_mcmc(stateManager,
            backupManager,
            safeToExit,
            runSucceeded,
            config,
            state,
            treeWriter,
            stripes,
            stripeNum,
            tmpDir)


def resume_existing_run(stateManager, backupManager, safeToExit,
                        runSucceeded, config):
    # If error occurs, restore the backups and try again. Never try more than two
    # times, however -- if the primary file and the backup file both fail, the
    # error is unrecoverable.
    try:
        state = stateManager.load_state()
        treeWriter = TreeWriter(resumeRun=True)
    except BaseException:
        logmsg('Restoring state failed:', sys.stderr)
        traceback.print_exc()
        logmsg('Restoring from backup and trying again.', sys.stderr)
        backupManager.restore_backup()

        state = stateManager.load_state()
        treeWriter = TreeWriter(resumeRun=True)

    np.random.set_state(state['rand_state'])  # Restore NumPy's RNG state.

    stripes, baseline = load_data(state['stripes_file'])

    print [sp.tag for sp in stripes]

    stripeNum = len(stripes)

    do_mcmc(stateManager,
            backupManager,
            safeToExit,
            runSucceeded,
            config,
            state,
            treeWriter,
            stripes,
            stripeNum,
            state['tmp_dir'])


def do_mcmc(stateManager,
            backupManager,
            safeToExit,
            runSucceeded,
            config,
            state,
            treeWriter,
            stripes,
            stripeNum,
            tmpDir):
    # here, the start iteration is saved in the state of 'last_iteration'
    # It should resample last incomplete iteration.
    startIter = state['last_iteration']

    unwrittenTreeL = []
    mcmcSampleTimesL = []
    lastMcmcSampleTime = time.time()

    # If --tmp-dir is not specified on the command line, it will by default be
    # None, which will cause mkdtemp() to place this directory under the system's
    # temporary directory. This is the desired behaviour.
    config['tmp_dir'] = tempfile.mkdtemp(prefix='pwgsdataexchange.', dir=tmpDir)

    totalIter = state['total_iteration']

    def proceedTime(timeTag, uSupportiveRanges, isBurnIn, args):
        # Here, at the outer loop, time tag is 0, start from the last iteration
        # Since we have separate resample_assignment process into multiple
        # iteration.

        totalIter, state, unwrittenTreeL, mcmcSampleTimesL, lastMcmcSampleTime,\
            config, stateManager, backupManager = (item[0] for item in args)


        if 0 == timeTag:
            iterRanges = range(startIter, state['sample_number'])
            if 0 < startIter:
                isBurnIn = False
            else:
                isBurnIn = True
        else:
            iterRanges = range(-state['burnin_sample_number'], state['sample_number'])

        for iteration in iterRanges:
            if 0 == timeTag:
                state['last_iteration'] = iteration
                state['total_iteration'] = totalIter

            tssb = state['tssb']
            safeToExit.set()

            # here resample
            if timeTag > state['time_tags'][0]:
                show_tree_structure(tssb, "iter{0}_time{1}".format(iteration,
                                                                   timeTag))

            tssb.resample_assignments(timeTag, uSupportiveRanges)

            show_tree_structure(tssb, "iter{0}_time{1}".format(iteration,
                                                               timeTag))

            if timeTag < state['time_tags'][-1]:
                tssb.mark_time_tag(timeTag)
                show_tree_structure(tssb, "iter{0}_time{1}".format(iteration,
                                                                   timeTag))
                uNext = deepcopy(uSupportiveRanges)
                uNext.remove(tssb.get_u_segL())

                argsNext = ([totalIter], [state], [unwrittenTreeL],
                        [mcmcSampleTimesL], [lastMcmcSampleTime], [config],
                        [stateManager], [backupManager])

                if not isBurnIn and 0 < Iteration:
                    proceedTime(timeTag+1, uNext, False, argsNext)
                else:
                    proceedTime(timeTag+1, uNext, True, argsNext)
            else:
                if isBurnIn:
                    logmsg(totalIter)

                tssb.cull_tree()

                # assign node ids
                wts, nodes = tssb.get_mixture()
                for i, node in enumerate(nodes):
                    node.id = i

                ##################################################
                # some useful info about the tree,
                # used by CNV related computations,
                # to be called only after resampling assignments
                set_node_height(tssb)
                set_path_from_root_to_node(tssb)
                map_datum_to_node(tssb)
                ##################################################

                state['mh_acc'] = metropolis(
                    tssb,
                    state['mh_itr'],
                    state['mh_std'],
                    state['mh_burnin'],
                    stripeNum,
                    state['stripes_text_file'],
                    state['rand_seed'],
                    state['max_copy_number'],
                    state['baseline'],
                    config['tmp_dir'])

                if float(state['mh_acc']) < 0.08 and state['mh_std'] < 10000:
                    state['mh_std'] = state['mh_std'] * 2.0
                    logmsg("Shrinking MH proposals. Now %f" % state['mh_std'])
                if float(state['mh_acc']) > 0.5 and float(state['mh_acc']) < 0.99:
                    state['mh_std'] = state['mh_std'] / 2.0
                    logmsg("Growing MH proposals. Now %f" % state['mh_std'])

                tssb.resample_sticks()
                tssb.resample_stick_orders()
                tssb.resample_hypers(dpAlpha=True, alphaDecay=True, dpGamma=True)

                lastLlh = tssb.complete_data_log_likelihood()

                if not isBurnIn:
                    cdltIdx = totalIter - state['burnin_cd_llh_traces'].shape[0]

                    state['cd_llh_traces'][cdltIdx] = lastLlh
                    if True or mod(cdltIdx, 10) == 0:
                        weights, nodes = tssb.get_mixture()
                        logmsg(' '.join([
                            str(v)
                            for v in (cdltIdx, len(nodes),
                                    state['cd_llh_traces'][cdltIdx],
                                    state['mh_acc'], tssb.dpAlpha, tssb.dpGamma,
                                    tssb.alphaDecay)
                        ]))
                    if argmax(state['cd_llh_traces'][:cdltIdx + 1]) == cdltIdx:
                        logmsg("%f is best per-data complete data likelihood so far." %
                            (state['cd_llh_traces'][cdltIdx]))
                else:
                    # here the save the llh
                    state['burnin_cd_llh_traces'][totalIter] = lastLlh

                # Can't just put tssb in unwrittenTreeL, as this object will be modified
                # on subsequent iterations, meaning any stored references in
                # unwrittenTreeL will all point to the same sample.
                serialized = pickle.dumps(tssb, protocol=pickle.HIGHEST_PROTOCOL)
                unwrittenTreeL.append((serialized, totalIter, lastLlh))
                state['tssb'] = tssb
                state['rand_state'] = np.random.get_state()

                if len([
                        C for C in state['tssb'].root['children']
                        if C['node'].has_data()
                ]) > 1:
                    logmsg('Polyclonal tree detected with %s clones.' % len(
                        state['tssb'].root['children']))

                newMcmcSampleTime = time.time()
                mcmcSampleTimesL.append(newMcmcSampleTime - lastMcmcSampleTime)
                lastMcmcSampleTime = newMcmcSampleTime

                # It's not safe to exit while performing file IO, as we don't want
                # trees.zip or the computation state file to become corrupted from an
                # interrupted write.
                safeToExit.clear()
                # here remove
                # shouldWriteBackup = totalIter % state['write_backups_every'] == 0 and totalIter != startIter
                shouldWriteBackup = totalIter % state['write_backups_every'] == 0
                shouldWriteState = totalIter % state['write_state_every'] == 0
                isLastIteration = (totalIter == state['cd_llh_traces'].shape[0]
                                   - 1 + state['burnin_cd_llh_traces'].shape[0])

                # If backup is scheduled to be written, write both it and full program
                # state regardless of whether we're scheduled to write state this
                # totalIter.
                if shouldWriteBackup or shouldWriteState or isLastIteration:
                    with open('mcmc_samples.txt', 'a') as mcmcf:
                        llhsAndTimes = [
                            (itr, llh, itr_time)
                            for (tssb, itr, llh
                                ), itr_time in zip(unwrittenTreeL, mcmcSampleTimesL)
                        ]
                        llhsAndTimes = '\n'.join([
                            '%s\t%s\t%s' % (itr, llh, itr_time)
                            for itr, llh, itr_time in llhsAndTimes
                        ])
                        mcmcf.write(llhsAndTimes + '\n')
                    treeWriter.write_trees(unwrittenTreeL)
                    stateManager.write_state(state)
                    unwrittenTreeL = []
                    mcmcSampleTimesL = []
                    if shouldWriteBackup:
                        backupManager.save_backup()

                totalIter = totalIter + 1

    args = ([totalIter], [state], [unwrittenTreeL], [mcmcSampleTimesL],
            [lastMcmcSampleTime], [config], [stateManager], [backupManager])
    proceedTime(state['time_tags'][0], MultiRangeSampler(0,1), True, args)

    backupManager.remove_backup()
    safeToExit.clear()
    # save the best tree
    print_top_trees(TreeWriter.defaultArchiveFn, state['top_k_trees_file'],
                    state['top_k'])

    # save clonal frequencies
    freq = dict([(g, []) for g in state['stripe_list']])
    stripeL = array(freq.keys(), str)
    stripeL.shape = (1, len(stripeL))
    savetxt(
        state['clonal_freqs_file'],
        vstack((stripeL, array([freq[g] for g in freq.keys()]).T)),
        fmt='%s',
        delimiter=', ')

    safeToExit.set()
    runSucceeded.set()


def test():
    tssb = pickle.load(open('ptree'))
    wts, nodes = tssb.get_mixture()
    for dat in tssb.data:
        print [dat.id, dat.__log_likelihood__(0.5)]


def run(args, safeToExit, runSucceeded, config):
    stateManager = StateManager()
    backupManager = BackupManager(
        [StateManager.defaultLastStateFn, TreeWriter.defaultArchiveFn])

    if stateManager.state_exists():
        logmsg('Resuming existing run. Ignoring command-line parameters.')
        resume_existing_run(stateManager, backupManager, safeToExit,
                            runSucceeded, config)
    else:
        # Ensure input files exist and can be read.
        try:
            stripesFile = open(args.stripesFile)
            stripesFile.close()
        except IOError as e:
            sys.stderr.write(str(e) + '\n')
            sys.exit(1)

        start_new_run(
            stateManager,
            backupManager,
            safeToExit,
            runSucceeded,
            config,
            args.stripesFile,
            args.stripesTextFile,
            args.paramsFile,
            topKTreesFile=args.topKTrees,
            clonalFreqsFile=args.clonalFreqs,
            burninSampleNum=args.burninSampleNum,
            sampleNum=args.mcmcSampleNum,
            mhItr=args.mhIterations,
            mhStd=100,
            writeStateEvery=args.writeStateEvery,
            writeBackupsEvery=args.writeBackupsEvery,
            randSeed=args.randomSeed,
            tmpDir=args.tmpDir,
            maxCopyNumber=args.maxCopyNumber)


def remove_tmp_files(tmpDir):
    if tmpDir is None:
        return
    tmpFileNames = get_c_fnames(tmpDir)
    for tmpFN in tmpFileNames:
        try:
            os.remove(tmpFN)
        except OSError:
            pass
    try:
        os.rmdir(tmpDir)
    except OSError:
        pass


def process(args):
    # Introducing threading is necessary to allow write operations to complete
    # when interrupts are received. As the interrupt halts execution of the main
    # thread and immediately jumps to the interrupt handler, we must run the
    # PhyloWGS code in a different thread, which clears the safeToExit flag
    # when in the midst of a write operation. This way, the main thread is left
    # only to handle the signal, allowing the derived thread to finish its
    # current write operation.
    safeToExit = threading.Event()

    # This will allow us to detect whether the run thread exited cleanly or not.
    # This means we can properly report a non-zero exit code if something failed
    # (e.g., something threw an exception). This is necessary because exceptions
    # in the run thread will terminate it, but can't be detected from the main
    # thread. A more robust strategy is here: http://stackoverflow.com/a/2830127.
    # Our strategy should be sufficient for the moment, though.
    runSucceeded = threading.Event()

    # We must know where temporary files are stored from within main() so that we
    # can remove them when we exit. However, we don't know this location until
    # the run thread starts, as when PWGS resumes an existing run, the parent
    # directory for the temporary files is stored in the state pickle file. Thus,
    # the run thread will set this value once it is established.
    #
    # So long as this dictionary is used only as a key-value store for primitive
    # objects, it's thread safe and doesn't require the use of a mutex. See
    # http://effbot.org/pyfaq/what-kinds-of-global-value-mutation-are-thread-safe.htm.
    # If more complex values are stored here, we must introduce a mutex.
    config = {'tmp_dir': None}

    def sigterm_handler(_signo, _stack_frame):
        logmsg('Signal %s received.' % _signo, sys.stderr)
        safeToExit.wait()
        remove_tmp_files(config['tmp_dir'])
        logmsg('Exiting now.')
        # Exit with non-zero to indicate run didn't finish.
        sys.exit(3)

    # SciNet will supposedly send SIGTERM 30 s before hard-killing the process.
    # This gives us time to clean up.
    signal.signal(signal.SIGTERM, sigterm_handler)
    # SIGINT is sent on CTRL-C. We don't want the user to interrupt a write
    # operation by hitting CTRL-C, thereby potentially resulting in corrupted
    # data being written. Permit these operations to finish before exiting.
    signal.signal(signal.SIGINT, sigterm_handler)

    runThread = threading.Thread(
        target=run, args=(args, safeToExit, runSucceeded, config))
    # Thread must be a daemon thread, or sys.exit() will wait until the thread
    # finishes execution completely.
    runThread.daemon = True
    runThread.start()

    while True:
        if not runThread.is_alive():
            break
        # I don't fully understand this. At least on the imacvm machine, calling
        # join with no timeout argument doesn't work, as the signal handler does
        # not seem to run until runThread exits. If, however, I specify *any*
        # timeout, no matter how short or long, the signal handler will run
        # *immediately* when the signal is sent -- i.e., even before the timeout
        # has expired.
        runThread.join(10)

    remove_tmp_files(config['tmp_dir'])
    if runSucceeded.is_set():
        logmsg('Run succeeded.')
        sys.exit(0)
    else:
        logmsg('Run failed.')
        sys.exit(1)


def logmsg(msg, fd=sys.stdout):
    print >> fd, '[%s] %s' % (datetime.now(), msg)


################################################################################
#                              function for debug                              #
################################################################################


def show_tree_structure(tssb, texFileName):
    print_tree_latex(tssb, texFileName+".tex", 0)
    command = "/usr/local/texlive/2017/bin/x86_64-linux/pdflatex {0}.tex && /usr/bin/okular {0}.pdf 2>&1 >/dev/null &".format(texFileName)
    os.system(command)

def modify_stripes_time_tag(stripes, times):
    sliceLen = round(len(stripes) / times)
    for index, stripe in zip(range(len(stripes)), stripes):
        stripe.tag = str(int(index / sliceLen))


def main():
    process()

if __name__ == "__main__":
    main()

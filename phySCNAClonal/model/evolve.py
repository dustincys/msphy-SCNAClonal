#!/usr/bin/env python2
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
import os
import sys
import cPickle as pickle

from numpy import *
from numpy.random import *
from tssb import *
from stripenode import *
from util import *

import numpy.random

from util2 import *
from params import *
from printo import *

import argparse
import signal
import tempfile
import threading
import traceback
import time
import json
from datetime import datetime

# numSamples: number of MCMC samples
# mhItr: number of metropolis-hasting iterations
# randSeed: random seed (initialization). Set to None to choose random
# seed automatically.


def start_new_run(stateManager,
                  backupManager,
                  safeToExit,
                  runSucceeded,
                  config,
                  stripesFile,
                  paramsFile,
                  topKTreesFile,
                  clonalFreqsFile,
                  burninSamples,
                  numSamples,
                  mhItr,
                  mhStd,
                  writeStateEvery,
                  writeBackupsEvery,
                  randSeed,
                  tmpDir):
    state = {}

    try:
        state['randSeed'] = int(randSeed)
    except TypeError:
        # If randSeed is not provided as command-line arg, it will be None,
        # meaning it will hit this code path.
        #
        # Use random seed in this order:
        #   1. If a seed is given on the command line, use that.
        #   2. Otherwise, if `random_seed.txt` exists, use the seed stored there.
        #   3. Otherwise, choose a new random seed and write to random_seed.txt.
        try:
            with open('random_seed.txt') as seedf:
                state['randSeed'] = int(seedf.read().strip())
        except (TypeError, IOError) as E:
            # Can seed with [0, 2**32).
            state['randSeed'] = randint(2**32)

    seed(state['randSeed'])
    with open('random_seed.txt', 'w') as seedf:
        seedf.write('%s\n' % state['randSeed'])

    state['stripesFile'] = stripesFile
    state['tmpDir'] = tmpDir
    state['topKTreesFile'] = topKTreesFile
    state['clonalFreqsFile'] = clonalFreqsFile
    state['writeStateEvery'] = writeStateEvery
    state['writeBackupsEvery'] = writeBackupsEvery

    # 此处载入数据
    stripes, baseline = load_data(state['stripesFile'])
    stripeNum = len(stripes)

    if len(stripes) == 0:
        logmsg('No stripes provided. Exiting.', sys.stderr)
        return

    # sample 个数
    # NTPS = len(codes[0].a)  # number of samples / time point
    # gene list
    # state['glist'] = [datum.name for datum in codes if len(datum.name) > 0]

    # stripe list
    state['stripeL'] = [stripe.stripe_name for stripe in stripes]

    # MCMC settings
    state['burnin'] = burninSamples
    state['numSamples'] = numSamples
    state['dpAlpha'] = 25.0
    state['dpGamma'] = 1.0
    state['alphaDecay'] = 0.25
    state['topK'] = 5

    # Metropolis-Hastings settings
    state['mhBurnin'] = 0
    state['mhItr'] = mhItr  # No. of iterations in metropolis-hastings
    state['mhStd'] = mhStd

    state['cdLlhTraces'] = zeros((state['numSamples'], 1))
    state['burninCdLlhTraces'] = zeros((state['burnin'], 1))
    state['workingDirectory'] = os.getcwd()

    root = StripeNode(conc=0.1)

    state['tssb'] = TSSB(
        dpAlpha=state['dpAlpha'],
        dpGamma=state['dpGamma'],
        alphaDecay=state['alphaDecay'],
        rootNode=root,
        data=stripes)
    # hack...
    if 1:
        depth = 0
        state['tssb'].root['sticks'] = vstack([
            state['tssb'].root['sticks'],
            boundbeta(1, state['tssb'].dpGamma) if depth != 0 else .999999
        ])
        state['tssb'].root['children'].append({
            'node':
            state['tssb'].root['node'].spawn(),
            'main':
            boundbeta(1.0, (state['tssb'].alphaDecay**
                            (depth + 1)) * state['tssb'].dpAlpha)
            if state['tssb'].min_depth <= (depth + 1) else 0.0,
            'sticks':
            empty((0, 1)),
            'children': []
        })
        newNode = state['tssb'].root['children'][0]['node']
        for n in range(state['tssb'].num_data):
            state['tssb'].assignments[n].remove_datum(n)
            newNode.add_datum(n)
            state['tssb'].assignments[n] = newNode

    for stripe in stripes:
        stripe.tssb = state['tssb']

    treeWriter = TreeWriter()
    # treeWriter.add_extra_file('cnv_logical_physical_mapping.json',
                               # json.dumps(cnv_logical_physical_mapping))

    if paramsFile is not None:
        with open(paramsFile) as F:
            params = json.load(F)
    else:
        params = {}
    treeWriter.add_extra_file('params.json', json.dumps(params))

    stateManager.write_initial_state(state)

    logmsg("Starting MCMC run...")
    state['last_iteration'] = -state['burnin'] - 1

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
        treeWriter = TreeWriter(resume_run=True)
    except BaseException:
        logmsg('Restoring state failed:', sys.stderr)
        traceback.print_exc()
        logmsg('Restoring from backup and trying again.', sys.stderr)
        backupManager.restore_backup()

        state = stateManager.load_state()
        treeWriter = TreeWriter(resume_run=True)

    set_state(state['rand_state'])  # Restore NumPy's RNG state.

    stripes, baseline = load_data(state['stripes_file'])
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
            state['tmpDir'])


def do_mcmc(stateManager,
            backupManager,
            safeToExit,
            runSucceeded,
            config,
            state,
            treeWriter,
            stripes,
            stripeNum,
            tmp_dir_parent):
    startIter = state['last_iteration'] + 1
    unwrittenTreeL = []
    mcmcSampleTimesL = []
    lastMcmcSampleTime = time.time()

    # If --tmp-dir is not specified on the command line, it will by default be
    # None, which will cause mkdtemp() to place this directory under the system's
    # temporary directory. This is the desired behaviour.
    config['tmpDir'] = tempfile.mkdtemp(
        prefix='pwgsdataexchange.', dir=tmp_dir_parent)

    for iteration in range(startIter, state['numSamples']):
        safeToExit.set()
        if iteration < 0:
            logmsg(iteration)

        # Referring to tssb as local variable instead of dictionary element is much
        # faster.
        tssb = state['tssb']
        tssb.resample_assignments()
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
            state['mhItr'],
            state['mhStd'],
            state['mhBurnin'],
            stripeNum,
            state['stripesFile'],
            state['randSeed'],
            config['tmpDir'])

        if float(state['mh_acc']) < 0.08 and state['mhStd'] < 10000:
            state['mhStd'] = state['mhStd'] * 2.0
            logmsg("Shrinking MH proposals. Now %f" % state['mhStd'])
        if float(state['mh_acc']) > 0.5 and float(state['mh_acc']) < 0.99:
            state['mhStd'] = state['mhStd'] / 2.0
            logmsg("Growing MH proposals. Now %f" % state['mhStd'])

        tssb.resample_sticks()
        tssb.resample_stick_orders()
        tssb.resample_hypers(dpAlpha=True, alphaDecay=True, dpGamma=True)

        lastLlh = tssb.complete_data_log_likelihood()
        if iteration >= 0:
            state['cdLlhTraces'][iteration] = lastLlh
            if True or mod(iteration, 10) == 0:
                weights, nodes = tssb.get_mixture()
                logmsg(' '.join([
                    str(v)
                    for v in (iteration, len(nodes),
                              state['cdLlhTraces'][iteration],
                              state['mh_acc'], tssb.dpAlpha, tssb.dpGamma,
                              tssb.alphaDecay)
                ]))
            if argmax(state['cdLlhTraces'][:iteration + 1]) == iteration:
                logmsg("%f is best per-data complete data likelihood so far." %
                       (state['cdLlhTraces'][iteration]))
        else:
            state['burninCdLlhTraces'][iteration
                                          + state['burnin']] = lastLlh

        # Can't just put tssb in unwrittenTreeL, as this object will be modified
        # on subsequent iterations, meaning any stored references in
        # unwrittenTreeL will all point to the same sample.
        serialized = pickle.dumps(tssb, protocol=pickle.HIGHEST_PROTOCOL)
        unwrittenTreeL.append((serialized, iteration, lastLlh))
        state['tssb'] = tssb
        state['rand_state'] = get_state()
        state['last_iteration'] = iteration

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
        shouldWriteBackup = iteration % state['writeBackupsEvery'] == 0 and iteration != startIter
        shouldWriteState = iteration % state['writeStateEvery'] == 0
        isLastIteration = (iteration == state['numSamples'] - 1)

        # If backup is scheduled to be written, write both it and full program
        # state regardless of whether we're scheduled to write state this
        # iteration.
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

    backupManager.remove_backup()
    safeToExit.clear()
    # save the best tree
    print_top_trees(TreeWriter.default_archive_fn, state['topKTreesFile'],
                    state['topK'])

    # save clonal frequencies
    freq = dict([(g, []) for g in state['stripeL']])
    stripeL = array(freq.keys(), str)
    stripeL.shape = (1, len(stripeL))
    savetxt(
        state['clonalFreqsFile'],
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


def parse_args():
    parser = argparse.ArgumentParser(
        description=
        'Run phySCNAClonal to infer subclonal composition from SCNA stripes',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-b',
        '--write-backups-every',
        dest='writeBackupsEvery',
        default=100,
        type=int,
        help=
        'Number of iterations to go between writing backups of program state')
    parser.add_argument(
        '-S',
        '--write-state-every',
        dest='writeStateEvery',
        default=10,
        type=int,
        help=
        'Number of iterations between writing program state to disk. Higher values reduce IO burden at the cost of losing progress made if program is interrupted.'
    )
    parser.add_argument(
        '-k',
        '--top-k-trees',
        dest='top_k_trees',
        default='top_k_trees',
        help='Output file to save top-k trees in text format')
    parser.add_argument(
        '-f',
        '--clonal-freqs',
        dest='clonal_freqs',
        default='clonalFrequencies',
        help='Output file to save clonal frequencies')
    parser.add_argument(
        '-B',
        '--burnin-samples',
        dest='burninSamples',
        default=1000,
        type=int,
        help='Number of burnin samples')
    parser.add_argument(
        '-s',
        '--mcmc-samples',
        dest='mcmc_samples',
        default=2500,
        type=int,
        help='Number of MCMC samples')
    parser.add_argument(
        '-i',
        '--mh-iterations',
        dest='mh_iterations',
        default=5000,
        type=int,
        help='Number of Metropolis-Hastings iterations')
    parser.add_argument(
        '-r',
        '--random-seed',
        dest='random_seed',
        type=int,
        help='Random seed for initializing MCMC sampler')
    parser.add_argument(
        '-t',
        '--tmp-dir',
        dest='tmpDir',
        help='Path to directory for temporary files')
    parser.add_argument(
        '-p',
        '--params',
        dest='paramsFile',
        help='JSON file listing run parameters, generated by the parser')
    parser.add_argument(
        'stripesFile',
        help=
        'File listing stripes(SCNA stripes). For proper format, see README.md.')
    args = parser.parse_args()
    return args


def run(safeToExit, runSucceeded, config):
    stateManager = StateManager()
    backupManager = BackupManager(
        [StateManager.default_last_state_fn, TreeWriter.default_archive_fn])

    if stateManager.state_exists():
        logmsg('Resuming existing run. Ignoring command-line parameters.')
        resume_existing_run(stateManager, backupManager, safeToExit,
                            runSucceeded, config)
    else:
        args = parse_args()
        # Ensure input files exist and can be read.
        try:
            stripesFile = open(args.stripesFile)
            stripe.close()
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
            args.paramsFile,
            topKTreesFile=args.top_k_trees,
            clonalFreqsFile=args.clonal_freqs,
            burninSamples=args.burninSamples,
            numSamples=args.mcmc_samples,
            mhItr=args.mh_iterations,
            mhStd=100,
            writeStateEvery=args.writeStateEvery,
            writeBackupsEvery=args.writeBackupsEvery,
            randSeed=args.random_seed,
            tmpDir=args.tmpDir)


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


def main():
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
    config = {'tmpDir': None}

    def sigterm_handler(_signo, _stack_frame):
        logmsg('Signal %s received.' % _signo, sys.stderr)
        safeToExit.wait()
        remove_tmp_files(config['tmpDir'])
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
        target=run, args=(safeToExit, runSucceeded, config))
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

    remove_tmp_files(config['tmpDir'])
    if runSucceeded.is_set():
        logmsg('Run succeeded.')
        sys.exit(0)
    else:
        logmsg('Run failed.')
        sys.exit(1)


def logmsg(msg, fd=sys.stdout):
    print >> fd, '[%s] %s' % (datetime.now(), msg)


if __name__ == "__main__":
    main()

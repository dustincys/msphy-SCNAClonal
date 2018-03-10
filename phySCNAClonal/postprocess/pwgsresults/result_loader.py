import gzip
import json
import zipfile


class ResultLoader(object):
    def __init__(self, treeSummFn, mutListFn, mutAssignFn):
        self._treeSummFn = treeSummFn
        self._mutListFn = mutListFn
        self._mutAssignFn = mutAssignFn

        self.mutlist = None
        self.treeSumm = None
        self.datasetName = None

        self._load_tree_data()

    def _convert_keys_to_ints(self, dic):
        keys = dic.keys()
        for key in dic.keys():
            dic[int(key)] = dic[key]
            del dic[key]

    def _load_tree_data(self):
        with gzip.GzipFile(self._treeSummFn) as treesummf:
            treeJson = json.load(treesummf)
            self.datasetName = treeJson['datasetName']
            self.treeSumm = treeJson['trees']
            self.params = treeJson['params']

            if 'tree_densities' in treeJson:
                self.tree_densities = treeJson['tree_densities']
                self._convert_keys_to_ints(self.tree_densities)
            else:
                self.tree_densities = {}

        self._convert_keys_to_ints(self.treeSumm)
        for treeIdx, treeFeatures in self.treeSumm.items():
            self._convert_keys_to_ints(treeFeatures['populations'])
            self._convert_keys_to_ints(treeFeatures['structure'])

        with gzip.GzipFile(self._mutListFn) as mutlistf:
            self.mutlist = json.load(mutlistf)
        self.num_ssms = len(self.mutlist['ssms'])

    def _load_assignments(self, mutf, treeIdx):
        mutass = json.loads(mutf.read('%s.json' % treeIdx))
        mutass = mutass['mut_assignments']
        self._convert_keys_to_ints(mutass)
        return mutass

    def load_mut_assignments(self, treeIdx):
        with zipfile.ZipFile(self._mutAssignFn) as mutf:
            return self._load_assignments(mutf, treeIdx)

    def load_all_mut_assignments(self):
        with zipfile.ZipFile(self._mutAssignFn) as mutf:
            treeIndices = [
                int(i.filename.split('.')[0]) for i in mutf.infolist()
            ]
            treeIndices.sort()
            for treeIdx in treeIndices:
                yield (treeIdx, self._load_assignments(mutf, treeIdx))

    def load_all_mut_assignments_into_memory(self):
        mutass = {I: M for (I, M) in self.load_all_mut_assignments()}
        return mutass

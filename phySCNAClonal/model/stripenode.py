from scipy.stats import beta, binom
import scipy.stats as stat
from scipy.misc import comb
from util import *
from numpy import *
from node import *

from util2 import *


class StripeNode(Node):

    init_mean = 0.5
    min_conc = 0.01
    max_conc = 0.1

    def __init__(self, parent=None, tssb=None, conc=0.1):
        super(alleles, self).__init__(parent=parent, tssb=tssb)

        if tssb is not None:
            ntps = len(tssb.data[0].a)

        # pi is a first-class citizen
        self.pi = 0.0
        self.param = 0.0
        self.param1 = 0.0
        self.pi1 = 0.0  # used in MH	to store old state

        self.path = None  # set of nodes from root to this node
        self.ht = 0.0

        if parent is None:
            self._conc = conc
            self.pi = 1.0
            self.param = 1.0

        else:
            self.pi = rand(1)*parent.pi
            parent.pi = parent.pi - self.pi
            self.param = self.pi

    def conc(self):
        if self.parent() is None:
            return self._conc
        else:
            return self.parent().conc()

    def kill(self):
        if self._parent is not None:
            self._parent._children.remove(self)
        self._parent.pi = self._parent.pi + self.pi
        self._parent = None
        self._children = None

    def logprob(self, x):
        return x[0]._log_likelihood(self.param)

    def logprob_restricted(self, x):
        lower_node, upper_node = self.find_neighbor_datum_n(x)

        l_flag = True
        u_flag = True
        if lower_node is not None:
            l_flag = self.__is_good_gap(lower_node, x, "lower")
        else:
            l_flag = True

        if upper_node is not None:
            u_flag = self.__is_good_gap(lower_node, x, "upper")
        else:
            u_flag = True
        if l_flag and u_flag:
            return self.logprob(x)
        else:
            return -float('Inf')

    def complete_logprob(self):
        return sum([self.logprob([data]) for data in self.get_data()])

    def find_neighbor_datum_n(self, x):
        datums = self.get_data()
        if x not in datums:
            datums.append(x)

        datums_sorted = sorted(datums,
            key=lambda item: 1.0*item.tumor_reads_num/item.normal_reads_num)

        idx = datums_sorted.index(x)
        if 0 == idx:
            return (None, datums_sorted[1])
        elif len(datums_sorted) - 1 == idx:
            return (datums_sorted[idx-1], None)
        else:
            return (datums_sorted[idx-1], datums_sorted[idx+1])

    def __is_good_gap(self, lower_node, upper_node, position):
        rdr_lower = 1.0*lower_node.tumor_reads_num/lower_node.normal_reads_num
        rdr_upper = 1.0*upper_node.tumor_reads_num/upper_node.normal_reads_num
        L = np.exp(rdr_upper - rdr_lower)

        if "lower" == position:
            cn = lower_node.copy_number
        elif "upper" == position:
            cn = upper_node.copy_number - 1

        if cn < 0:
            return False
        else:
            return L >= (1.0 + (self.param /
                                (cn * self.param + 2 * (1 - self.param))))

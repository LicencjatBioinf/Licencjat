import IsoSpecPy
import numpy as np
from numba import float32, double, njit  # import the types
from numba.experimental import jitclass




spec = [
    ('dist', float32),
    ('covered_probability', float32),
    ('node_matrix', double[:, :]),
]


class FormulaError(Exception):
    pass


@njit()
def bisect_left(a, x, lo=0, hi=None):
    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


@jitclass(spec)
class WGraph:
    def __init__(self, node_matrix):
        self.dist = 0
        self.covered_probability = 0.0
        self.node_matrix = node_matrix  # every row has a form [mass, probability, flow on the right edge], rows are
        # sorted by mass

    def insert_new_node(self, mass, probability):
        self.covered_probability += probability
        position = bisect_left(self.node_matrix[:, 0], mass)
        self.node_matrix = np.concatenate((self.node_matrix[:position, :], np.array([[
            mass, probability, self.node_matrix[position - 1, 2]]]), self.node_matrix[position:, :]))
        self.find_best_flow(position)

    def find_best_flow(self, position):
        left_index, right_index = position - 1, position + 1
        if position == 0:
            left_path_length = np.inf
            right_path_length = self.node_matrix[right_index, 0] - self.node_matrix[position, 0]
        elif position == np.shape(self.node_matrix)[0] - 1:
            left_path_length = self.node_matrix[position][0] - self.node_matrix[left_index][0]
            right_path_length = np.inf
        else:
            left_path_length = self.node_matrix[position][0] - self.node_matrix[left_index][0]
            right_path_length = self.node_matrix[right_index][0] - self.node_matrix[position][0]

        while self.node_matrix[position][1] > 0 and \
                (left_path_length != np.inf or right_path_length != np.inf):
            if left_path_length <= right_path_length:
                if self.node_matrix[left_index][1] > 0:
                    self.calculate_distance(start=left_index, stop=position, direction=-1)
                else:
                    old = left_index
                    left_index -= 1
                    left_path_length += self.node_matrix[old][0] - self.node_matrix[left_index][0] if left_index >= 0 \
                        else np.inf
            else:
                if self.node_matrix[right_index][1] > 0:
                    self.calculate_distance(start=position, stop=right_index)
                else:
                    old = right_index
                    right_index += 1
                    right_path_length += self.node_matrix[right_index][0] - self.node_matrix[old][0] if \
                        right_index <= np.shape(self.node_matrix)[0] - 1 else np.inf

    def calculate_distance(self, start, stop, direction=1):
        flow = np.minimum(self.node_matrix[stop, 1], self.node_matrix[start, 1])
        self.node_matrix[stop, 1] -= flow
        self.node_matrix[start, 1] -= flow
        edges_length = self.node_matrix[1:, 0] - self.node_matrix[:-1, 0]
        self.node_matrix[start:stop, 2] += flow * direction
        self.dist = np.sum(np.abs(self.node_matrix[:-1, 2] * edges_length))


class FlowGraph(object):
    __slots__ = ['WGraph', 'peakGenerator']

    def __init__(self, *args, fasta=None, node_matrix=None):
        try:
            self.peakGenerator = iter(IsoSpecPy.IsoLayeredGenerator(fasta=fasta))
        except Exception:
            raise FormulaError('invalid sequence') from None
        self.WGraph = WGraph(node_matrix)

    def add_new_peak(self):
        try:
            new_node = next(self.peakGenerator)
            self.WGraph.insert_new_node(*new_node)
        except Exception as e:
            print(e)

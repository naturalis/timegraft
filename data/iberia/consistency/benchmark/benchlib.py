"""Shared helpers of the benchmark of #56: tables, backbone trees, Faith's PD, rank correlation."""

import csv
import gzip
import math
from pathlib import Path

import dendropy

csv.field_size_limit(10**9)


def read(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows, columns=None):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    columns = columns or (list(rows[0]) if rows else [])
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def quantile(values, q):
    values = sorted(values)
    if not values:
        return math.nan
    pos = q * (len(values) - 1)
    lo, hi = math.floor(pos), math.ceil(pos)
    return values[lo] + (values[hi] - values[lo]) * (pos - lo)


def factor(a, b):
    """How far apart two positive numbers are, as a factor of at least 1."""
    return max(a / b, b / a)


class Backbone:
    """A dated tree: the age and the tips of every internal node, and the stem node of every tip."""

    def __init__(self, path):
        self.tree = dendropy.Tree.get(path=str(path), schema="newick", rooting="force-rooted",
                                      preserve_underscores=True)
        self.age, self.tips, self.parent = {}, {}, {}
        self.negative = 0
        for node in self.tree.postorder_node_iter():
            if node.is_leaf():
                self.age[node], self.tips[node] = 0.0, frozenset([node.taxon.label])
            else:
                children = node.child_nodes()
                self.age[node] = max(self.age[c] + (c.edge.length or 0.0) for c in children)
                self.tips[node] = frozenset().union(*(self.tips[c] for c in children))
                self.negative += sum(1 for c in children if (c.edge.length or 0.0) < -1e-9)
            for child in node.child_nodes():
                self.parent[child] = node
        self.leaf = {n.taxon.label: n for n in self.tree.leaf_node_iter()}
        self.by_tips = {self.tips[n]: n for n in self.tree.preorder_internal_node_iter()}

    def stem_age(self, label):
        return self.age[self.parent[self.leaf[label]]]

    def node_ages(self):
        """Age per clade (as its set of tips), internal nodes only."""
        return {tips: self.age[node] for tips, node in self.by_tips.items()}

    def mrca_age(self, labels):
        """Age of the smallest clade that holds all the labels (they must all be in the tree)."""
        labels = set(labels)
        node = self.leaf[next(iter(labels))]
        while not labels <= self.tips[node]:
            node = self.parent[node]
        return self.age[node]

    def smallest_clade_with_company(self, label):
        """The other tips of the smallest clade that holds the tip and at least one other."""
        node = self.parent[self.leaf[label]]
        return self.tips[node] - {label}

    def path_to_root(self, label):
        node = self.parent[self.leaf[label]]
        while node is not None:
            yield node
            node = self.parent.get(node)


class PDTree:
    """Faith's PD (with the path to the root) of sets of tips, on a tree read once."""

    def __init__(self, path):
        tree = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True,
                                 rooting="force-rooted")
        self.parent, self.length, self.tip = {}, {}, {}
        for i, node in enumerate(tree.preorder_node_iter()):
            node._i = i
        for node in tree.preorder_node_iter():
            if node.parent_node is not None:
                self.parent[node._i] = node.parent_node._i
                self.length[node._i] = node.edge.length or 0.0
            if node.is_leaf():
                self.tip[node.taxon.label] = node._i
        self.total = sum(self.length.values())

    def pd(self, labels):
        seen, total = set(), 0.0
        for label in labels:
            i = self.tip[label]
            while i in self.parent and i not in seen:
                seen.add(i)
                total += self.length[i]
                i = self.parent[i]
        return total


def spearman(x, y):
    def rank(values):
        order = sorted(range(len(values)), key=values.__getitem__)
        ranks = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            for k in range(i, j + 1):
                ranks[order[k]] = (i + j) / 2
            i = j + 1
        return ranks

    a, b = rank(x), rank(y)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    cov = sum((p - ma) * (q - mb) for p, q in zip(a, b))
    va, vb = sum((p - ma) ** 2 for p in a), sum((q - mb) ** 2 for q in b)
    return cov / math.sqrt(va * vb) if va and vb else math.nan

"""analysis.py -- small exact/plug-in information measures used by the tests."""
from collections import Counter
from math import log2


def entropy(counter):
    n = sum(counter.values())
    return -sum(v / n * log2(v / n) for v in counter.values())


def plugin_mi(xs, ys):
    """Plug-in mutual information (bits) between two equally long discrete sequences."""
    assert len(xs) == len(ys)
    return entropy(Counter(xs)) + entropy(Counter(ys)) - entropy(Counter(zip(xs, ys)))

from dataclasses import dataclass

import numpy as np

from GRMA.Utilities.cutils import chash


@dataclass
class ClassMinusOne:
    subclass: int
    class_num: int
    allele_num: int

    def __hash__(self):
        return self.subclass

    def __eq__(self, other):
        if isinstance(other, int):
            return self.subclass == other
        elif isinstance(other, ClassMinusOne):
            return self.subclass == other.subclass
        return False


class HashableArray:
    __slots__ = "arr", "it"

    def __init__(self, arr):
        self.it = None
        self.arr = np.array(arr, dtype=np.uint16)

    def __hash__(self):
        return chash(self.arr)

    def __add__(self, other):
        return HashableArray(np.concatenate([self.arr, other.np()], axis=0))

    def __radd__(self, other):
        return HashableArray(np.concatenate([other.np(), self.arr], axis=0))

    def __next__(self):
        return next(self.it)

    def __iter__(self):
        self.it = iter(self.arr)
        return self.it

    def __getitem__(self, item):
        return HashableArray(self.arr[item])

    def __len__(self):
        return len(self.arr)

    def __repr__(self):
        return str(self.arr)

    def np(self):
        return self.arr


import numpy as np


class Node(object):
    def __init__(self, idx: int, x: float, y: float, z: float):
        self.idx = idx
        self.coord = np.array([x, y, z])

        self._is_boundary = False

    @property
    def x(self):
        return self.coord[0]

    @property
    def y(self):
        return self.coord[1]

    @property
    def z(self):
        return self.coord[2]

    @property
    def is_boundary(self):
        return self._is_boundary

    @is_boundary.setter
    def is_boundary(self, value: bool):
        self._is_boundary = value

    def __repr__(self):
        return f"idx: {self.idx}, {self.coord}"

from Node import Node

"""
Node Order:
node0----node1----node2
"""


class Element_D1O2(object):
    def __init__(self, idx: int, node0: Node, node1: Node, node2: Node):
        self.idx = idx
        self.nodes = [node0, node1, node2]

        self._xi_a = 0.5 * (self.x0 + self.x2) - self.x1
        self._xi_b = 0.5 * (self.x2 - self.x0)
        self._xi_c = self.x1

    @property
    def node0(self):
        return self.nodes[0]

    @property
    def node1(self):
        return self.nodes[1]

    @property
    def node2(self):
        return self.nodes[2]

    @property
    def node0_idx(self):
        return self.node0.idx

    @property
    def node1_idx(self):
        return self.node1.idx

    @property
    def node2_idx(self):
        return self.node2.idx

    @property
    def has_boundary_node(self):
        return self.node0.is_boundary or self.node2.is_boundary

    @property
    def boundary_nodes(self):
        if self.has_boundary_node:
            return [self.node0, self.node2]
        else:
            return []

    @property
    def boundary_node_idxes(self):
        return [n.idx for n in self.boundary_nodes]

    @property
    def x0(self):
        return self.node0.x

    @property
    def x1(self):
        return self.node1.x

    @property
    def x2(self):
        return self.node2.x

    def __repr__(self):
        return f"idx: {self.idx}, nodes: [{self.node0_idx}, {self.node1_idx}, {self.node2_idx}]"

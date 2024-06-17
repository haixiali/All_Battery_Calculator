import numpy as np
import matplotlib.pyplot as plt

from objdict import ObjDict
from GmshMesh import GmshMesh
from palettable.colorbrewer.qualitative import Pastel1_9, Set1_9
from palettable.scientific.sequential import Acton_20


def readmsh(filename):
    return FEMesh(GmshMesh(filename))


class FEMesh(object):

    __slots__ = ['nodes', '_elements', 'groups', 'p', 'Ne', 'elements', 'edgeNodeIDs']

    def __init__(self, m: GmshMesh):

        # Nodes
        self.nodes = np.array(m.nodes).T
        if not self.nodes[2, :].any():
            self.nodes = self.nodes[0:2, :]

        # Elements
        self._elements = m.elements

        # Groups
        self.groups = ObjDict()

        def addGroup(name, blocks):
            self.groups[name] = Group(name, self, blocks)

        # Logical groups
        addGroup('edges', m.elementBlocksBy(entityDim=1))
        addGroup('faces', m.elementBlocksBy(entityDim=2))

        # Physical groups
        for n in m.physicalNamesList:
            bs = m.elementBlocksBy(physicalName=n)
            addGroup(n, bs)

        # Simple mode: provide easy access
        self.p = self.groups.faces.p
        self.Ne = self.groups.faces.Ne
        self.elements = self.groups.faces.elements
        self.edgeNodeIDs = self.groups.edges.nodeIDs

    @property
    def mode(self):
        if len(self.groups) > 2:
            return 'advanced'
        else:
            return 'simple'

    @property
    def Nn(self):
        return self.nodes.shape[1]

    @property
    def groupNames(self):
        return list(self.groups.keys())

    def findNodeAt(self, x1, x2):
        x = np.array([x1, x2])
        for i in range(0, self.Nn):
            if np.linalg.norm(self.nodes[:, i] - x) < 1e-12:
                return i
        raise Exception(f'No node at ({x1}, {x2}) found')

    def plot(self, scalars=None):

        # Collect groups to plot
        groupsToPlot = []
        for g in self.groups.values():
            if not g.isLogical or self.mode == 'simple':
                groupsToPlot.append(g)

        # Prepare plot
        fig1, ax1 = plt.subplots()
        ax1.set_aspect('equal')

        if scalars is None:

            # Faces
            cnt = 0
            for g in groupsToPlot:
                if g.dimension == 2:
                    c = g.Ne * [cnt]
                    ax1.tripcolor(
                        self.nodes[0, :], self.nodes[1, :], g.elements.T,
                        facecolors=c, cmap=Pastel1_9.mpl_colormap.reversed(), vmin=0, vmax=8,
                        edgecolor='black', linewidth=0.5
                    )
                    cnt += 1

            # Edges
            cnt = 0
            ec = np.array(Set1_9.colors) / 255
            for g in groupsToPlot:
                if g.dimension == 1:
                    x = self.nodes[0, g.elements]
                    y = self.nodes[1, g.elements]
                    plt.plot(x, y, c=ec[cnt, :])
                    cnt += 1

        else:
            cm = Acton_20.mpl_colormap.reversed()
            if len(scalars) == self.Nn:
                p = ax1.tripcolor(
                    self.nodes[0, :], self.nodes[1, :], self.elements.T, scalars, 
                    cmap=cm, shading='gouraud'
                )
                ax1.triplot(
                    self.nodes[0, :], self.nodes[1, :], self.elements.T, color='black', linewidth=0.5
                )
            else:
                p = ax1.tripcolor(
                    self.nodes[0, :], self.nodes[1, :], self.elements.T, scalars,
                    cmap=cm, edgecolor='black', linewidth=0.5
                )
            fig1.colorbar(p)

    def print(self):
        np.set_printoptions(linewidth=np.nan)
        print('Nodes')
        print(self.nodes)
        print('Elements')
        print(self._elements)
        print('Groups')
        for g in self.groups.values():
            g.print()


class Group(object):

    __slots__ = ['p', 'name', 'dimension', 'feMesh', '_elementIDs', '_nodeIDs', '_nodes', '_elements']

    def __init__(self, name, feMesh, blocks):

        # Basic properties
        self.p = ObjDict()
        self.name = name
        self.dimension = blocks[0].entityDimension

        # Pointer to Mesh
        self.feMesh = feMesh

        # Element IDs
        self._elementIDs = []
        for b in blocks:
            assert self.dimension == b.entityDimension
            self._elementIDs += b.tags

        # Cache for other things
        self._nodeIDs = None
        self._nodes = None
        self._elements = None

    @property
    def Nn(self):
        return len(self.nodeIDs)

    @property
    def nodeIDs(self):
        if self._nodeIDs is None:
            self._nodeIDs = list(np.unique(self.elements))
        return self._nodeIDs

    @property
    def nodes(self):
        if self._nodes is None:
            self._nodes = self.feMesh.nodes[:, self.nodeIDs]
        return self._nodes

    @property
    def Ne(self):
        return len(self._elementIDs)

    @property
    def elements(self):
        if self._elements is None:
            e = [self.feMesh._elements[i] for i in self._elementIDs]
            self._elements = np.array(e).T

        return self._elements

    @property
    def isLogical(self):
        return self.name in ['edges', 'faces']

    def print(self):
        print('\nGroup \'' + self.name + '\'')
        print('Dimension: ', self.dimension)
        print('NodeIDs:\n', self.nodeIDs)
        print('Nodes:\n', self.nodes)
        print('ElementIDs:\n', self._elementIDs)
        print('Elements:\n', self.elements)


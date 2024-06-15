class MeshFormat(object):
    def __init__(self, version, file_type, data_size):
        self.version = version
        self.file_type = file_type
        self.data_size = data_size

    def __repr__(self):
        return f"gmsh, version: {self.version}, file_type: {self.file_type}, data_size: {self.data_size}"


class GmshMesh(object):

    def __init__(self, filename) -> None:

        f = open(filename, 'r')
        line = f.readline().strip()

        # Read file
        while line:
            if line == '$MeshFormat':
                self.mesh_format = read_mesh_format(f)
            elif line == '$PhysicalNames':
                self.physicalNames = read_physical_names(f)
            elif line == '$Entities':
                self.entities = read_entities(f)
            elif line == '$Nodes':
                self.nodes = read_nodes(f)
            elif line == '$Elements':
                self.elements = read_elements(f)
            else:
                print('Warning: Unknown entry', line)
                readsection(f, line)

            # read section ending keyword
            cc = f.read(1)
            while cc and cc != '$':
                cc = f.read(1)
            f.readline()
            line = f.readline().strip()


"""
    def elementBlocksBy(self, entityDim=None, physicalName=None):
        blocks = []

        # By entity dimension
        if entityDim is not None:
            for eb in self.elementBlocks:
                if eb.entityDimension == entityDim:
                    blocks.append(eb)

        # By physical names
        if physicalName is not None:
            for eb in self.elementBlocks:
                e = self._entityBy(dim=eb.entityDimension, tag=eb.entityTag)
                for t in e.physicalTags:
                    pn = self._physicalNameBy(t)
                    if pn.name == physicalName:
                        blocks.append(eb)

        return blocks

    def _entityBy(self, dim, tag):
        for e in self.allEntities:
            if e.dimension == dim and e.tag == tag:
                return e

    def _physicalNameBy(self, tag):
        for pn in self.physicalNames:
            if pn.tag == tag:
                return pn

    @property
    def physicalNamesList(self):
        names = []
        for pn in self.physicalNames:
            names.append(pn.name)
        return names

    @property
    def allEntities(self):
        return self.entities.points + self.entities.curves + \
            self.entities.surfaces + self.entities.volumes
"""


def read_elements(f):
    # description of all element blocks
    num_blocks, num_elements, min_element_tag, max_element_tag = scan_line(f, '%d %d %d %d')

    # check numbering
    if max_element_tag - min_element_tag + 1 != num_elements:
        raise Exception('elements not continuously numbered')

    # return values
    elements = []
    # loop over element blocks
    for i in range(num_blocks):
        # description of element block
        entity_dim, entity_tag, element_type = scan_line(f, '%d %d %d')

        # number of elements in block
        num_elements_in_block = scan_line(f, '%d')[0]

        # read elements
        for j in range(num_elements_in_block):
            tag = scan_line(f, '%d')[0]
            node_tags = list(map(int, f.readline().split()))
            elements.append([tag, node_tags])

    # Return
    return elements


def read_nodes(f):
    # description of all node blocks
    num_blocks, num_nodes, min_node_tag, max_node_tag = scan_line(f, '%d %d %d %d')

    # check numbering
    if max_node_tag - min_node_tag + 1 != num_nodes:
        raise Exception('nodes not continuously numbered')

    # return values
    nodes = []

    # loop over node blocks
    for i in range(num_blocks):
        # description of node block
        entity_dim, entity_tag, parametric = scan_line(f, '%d %d %d')

        # number, tags and coordinates
        num_node_in_block = scan_line(f, '%d')[0]
        tags = scan_line(f, '%d', num_node_in_block)
        coordinates = scan_line(f, '%f %f %f', num_node_in_block)

        # checks
        if parametric != 0:
            raise Exception('parametric geometry not supported yet')

        # save values
        if num_node_in_block == 1:
            nodes.append([tags[0] - 1, coordinates])
        else:
            for j in range(num_node_in_block):
                nodes.append([tags[j][0] - 1, coordinates[j]])  # in gmsh, tag starts from 1

    # Return
    return nodes


def read_cur_sur_vol(f, num, dim):
    ent = []
    for i in range(num):
        tag, min_x, min_y, min_z, max_x, max_y, max_z = scan_line(f, '%d %f %f %f %f %f %f')
        num_phy_tags = scan_line(f, '%d')[0]
        phy_tags = scan_line(f, '%d ' * num_phy_tags)
        num_bounding = scan_line(f, '%d')[0]
        bounding_tags = scan_line(f, '%d ' * num_bounding)

        ent.append(
            [dim, tag, [min_x, min_y, min_z], [max_x, max_y, max_z], num_phy_tags, phy_tags, num_bounding,
             bounding_tags])

    return ent


def read_points(f, n):
    points = []
    for i in range(n):
        tag, x, y, z, num_phy_tags = scan_line(f, '%d %f %f %f %d')
        phy_tags = scan_line(f, '%d ' * num_phy_tags)

        points.append([tag, [x, y, z], num_phy_tags, phy_tags])
    return points


def read_entities(f):
    num_points, num_curves, num_surfaces, num_volumes = scan_line(f, '%d %d %d %d', 1)
    points = read_points(f, num_points)
    curves = read_cur_sur_vol(f, num_curves, 1)
    surfaces = read_cur_sur_vol(f, num_surfaces, 2)
    volumes = read_cur_sur_vol(f, num_volumes, 3)
    return points, curves, surfaces, volumes


def read_physical_names(f):
    num_physical_names = scan_line(f, '%d')[0]
    # dimension, physical_tag, name
    phy_names = scan_line(f, '%d %d %s', num_physical_names)
    return phy_names


def read_mesh_format(f):
    version, file_type, data_size = scan_line(f, '%f %d %d')
    mesh_format = MeshFormat(version, file_type, data_size)
    return mesh_format


def readsection(f, s):
    line = f.readline().strip()
    while line and line != '$End' + s[1:]:
        line = f.readline().strip()


def scan_line(f, format_str, n_lines=1):
    # token format
    token_format = [s for s in format_str.split()]

    # result
    result = []
    for i in range(n_lines):
        row = []
        for j in token_format:
            token = next_token(f, j)
            row.append(token)
        result.append(row)

    # Return
    if n_lines == 1:
        return result[0]
    else:
        return result


def next_token(f, token_type):
    token = ''

    # find beginning
    cc = f.read(1)
    while cc and cc.isspace():
        cc = f.read(1)

    if token_type == "%d" or token_type == "%f":
        # read until next whitespace
        while cc and not cc.isspace():
            token += cc
            cc = f.read(1)
    elif token_type == "%s":
        # read string in wrapped in ""
        if cc != "\"":
            raise Exception("string should be wrapped in \"\" but get " + cc)
        cc = f.read(1)  # read left "\""

        while cc and cc != "\"":  # read string
            token += cc
            cc = f.read(1)

    if token_type == '%d':
        ret = int(token)
    elif token_type == '%f':
        ret = float(token)
    elif token_type == '%s':
        ret = str(token)
    else:
        raise Exception('illegal format: ' + token_type)

    # return
    return ret

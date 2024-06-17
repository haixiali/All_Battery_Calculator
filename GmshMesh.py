class MeshFormat(object):
    def __init__(self, version, file_type, data_size):
        self.version = version
        self.file_type = file_type
        self.data_size = data_size

    def __repr__(self):
        return f"gmsh, version: {self.version}, file_type: {self.file_type}, data_size: {self.data_size}"


class PhysicalName(object):
    def __init__(self, dim, tag, name):
        self.dim = dim
        self.tag = tag
        self.name = name

    def __repr__(self):
        return f"(tag: {self.tag}, name: {self.name}, dim: {self.dim})"


class Entity123(object):
    def __init__(self, tag, min_x, min_y, min_z, max_x, max_y, max_z, physical_tags, bounding_tags):
        self.tag = tag
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z
        self.physical_tags = physical_tags
        self.bounding_tags = bounding_tags

    def __repr__(self):
        return (f"tag: {self.tag}, min_x: {self.min_x}, min_y: {self.min_y}, min_z: {self.min_z}, "
                f"max_x: {self.max_x}, max_y: {self.max_y}, max_z: {self.max_z}, "
                f"physical tags: {self.physical_tags}, bounding tags: {self.bounding_tags}")


class Point(object):
    def __init__(self, tag, x, y, z, physical_tags):
        self.tag = tag
        self.x = x
        self.y = y
        self.z = z
        self.physical_tags = physical_tags

    def __repr__(self):
        return f"(tag: {self.tag}, x: {self.x}, y: {self.y}, z: {self.z})"


class Curve(Entity123):
    def __init__(self, tag, min_x, min_y, min_z, max_x, max_y, max_z, physical_tags, point_tags):
        super().__init__(tag, min_x, min_y, min_z, max_x, max_y, max_z, physical_tags, point_tags)


class Surface(Entity123):
    def __init__(self, tag, min_x, min_y, min_z, max_x, max_y, max_z, physical_tags, curve_tags):
        super().__init__(tag, min_x, min_y, min_z, max_x, max_y, max_z, physical_tags, curve_tags)


class Volume(Entity123):
    def __init__(self, tag, min_x, min_y, min_z, max_x, max_y, max_z, physical_tags, surface_tags):
        super().__init__(tag, min_x, min_y, min_z, max_x, max_y, max_z, physical_tags, surface_tags)


class Entity(object):
    def __init__(self, points: list[Point], curves: list[Curve], surfaces: list[Surface], volumes: list[Volume]):
        self.points = points
        self.curves = curves
        self.surfaces = surfaces
        self.volumes = volumes

    def __repr__(self):
        return (f"Point: {len(self.points)}, Curve: {len(self.curves)}, "
                f"Surface: {len(self.surfaces)}, Volume: {len(self.volumes)}")


class Node(Point):
    def __init__(self, tag, x, y, z):
        super().__init__(tag, x, y, z, physical_tags=None)


class NodeBlock(object):
    def __init__(self, dim, tag, parametric, node_tags):
        self.dim = dim
        self.tag = tag
        self.parametric = parametric
        self.node_tags = node_tags

    def __repr__(self):
        return f"dim: {self.dim}, tag: {self.tag}, nodes number: {len(self.node_tags)}"


class Element(object):
    def __init__(self, tag, node_tags):
        self.tag = tag
        self.node_tags = node_tags

    def __repr__(self):
        return f"(tag: {self.tag}, node tags: {self.node_tags})"


class ElementBlock(object):
    def __init__(self, dim, tag, element_type, element_tags):
        self.dim = dim
        self.tag = tag
        self.element_type = element_type
        self.element_tags = element_tags

    def __repr__(self):
        return (f"dim: {self.dim}, tag: {self.tag}, "
                f"element type: {self.element_type}, elements number: {len(self.element_tags)}")


class GmshMesh(object):

    def __init__(self, filename) -> None:
        self.mesh_format = None
        self.physical_names = None
        self.entities = None
        self.nodes = None
        self.node_blocks = None
        self.elements = None
        self.element_blocks = None

        f = open(filename, 'r')
        line = f.readline().strip()

        # read file
        while line:
            if line == '$MeshFormat':
                self.mesh_format = read_mesh_format(f)
            elif line == '$PhysicalNames':
                self.physical_names = read_physical_names(f)
            elif line == '$Entities':
                self.entities = read_entities(f)
            elif line == '$Nodes':
                self.nodes, self.node_blocks = read_nodes(f)
            elif line == '$Elements':
                self.elements, self.element_blocks = read_elements(f)
            else:
                print('Warning: Unknown entry', line)
                read_section_end(f)

            # read section ending keyword
            read_section_end(f)

            line = f.readline().strip()

    def __repr__(self):
        return (f"format: {self.mesh_format}\n"
                f"physical names: {self.physical_names}\n"
                f"entities: {self.entities}\n"
                f"nodes: {len(self.nodes)}, {self.nodes}\n"
                f"elements: {len(self.elements)}, {self.elements}\n")


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
    elements = list[Element]()
    element_blocks = list[ElementBlock]()

    # loop over element blocks
    for i in range(num_blocks):
        # description of element block
        entity_dim, entity_tag, element_type = scan_line(f, '%d %d %d')

        # number of elements in block
        num_elements_in_block = scan_line(f, '%d')[0]

        # read elements
        element_tags = []
        for j in range(num_elements_in_block):
            tag = scan_line(f, '%d')[0] - 1  # tag starts from 1 in gmsh, convert to zero-based
            node_tags = list(map(int, f.readline().split()))
            one_element = Element(tag, node_tags)
            elements.append(one_element)

            element_tags.append(tag)  # tag starts from 1 in gmsh, convert it to zero-based

        one_block = ElementBlock(entity_dim, entity_tag, element_type, element_tags)
        element_blocks.append(one_block)

    return elements, element_blocks


def read_nodes(f):
    # description of all node blocks
    num_blocks, num_nodes, min_node_tag, max_node_tag = scan_line(f, '%d %d %d %d')

    # check numbering
    if max_node_tag - min_node_tag + 1 != num_nodes:
        raise Exception('nodes not continuously numbered')

    nodes = list[Node]()
    node_blocks = list[NodeBlock]()

    # loop over node blocks
    for i in range(num_blocks):
        # description of node block
        node_tags = []
        entity_dim, entity_tag, parametric = scan_line(f, '%d %d %d')
        num_node_in_block = scan_line(f, '%d')[0]

        if parametric != 0:
            raise Exception('parametric geometry not supported yet')

        for j in range(num_node_in_block):
            tag = scan_line(f, '%d')[0] - 1  # tag starts from 1 in gmsh, convert to zero-based
            node_tags.append(tag)
        for j in range(num_node_in_block):
            x, y, z = scan_line(f, '%f %f %f')
            one_node = Node(node_tags[j], x, y, z)
            nodes.append(one_node)

        one_block = NodeBlock(entity_dim, entity_tag, parametric, node_tags)
        node_blocks.append(one_block)

    # return
    return nodes, node_blocks


def read_entities123(f, num, dim):
    entities = []
    for i in range(num):
        tag, min_x, min_y, min_z, max_x, max_y, max_z = scan_line(f, '%d %f %f %f %f %f %f')
        num_phy_tags = scan_line(f, '%d')[0]
        phy_tags = scan_line(f, '%d ' * num_phy_tags)
        num_bounding = scan_line(f, '%d')[0]
        bounding_tags = scan_line(f, '%d ' * num_bounding)

        one = None
        if dim == 1:
            one = Curve(tag, min_x, min_y, min_z, max_x, max_y, max_z, phy_tags, bounding_tags)
        elif dim == 2:
            one = Surface(tag, min_x, min_y, min_z, max_x, max_y, max_z, phy_tags, bounding_tags)
        elif dim == 3:
            one = Volume(tag, min_x, min_y, min_z, max_x, max_y, max_z, phy_tags, bounding_tags)
        else:
            raise Exception('Dimension must be 1 or 2 or 3')

        entities.append(one)

    return entities


def read_points(f, n):
    points = list[Point]()
    for i in range(n):
        tag, x, y, z, num_phy_tags = scan_line(f, '%d %f %f %f %d')
        phy_tags = scan_line(f, '%d ' * num_phy_tags)
        one = Point(tag, x, y, z, phy_tags)
        points.append(one)

    return points


def read_entities(f):
    num_points, num_curves, num_surfaces, num_volumes = scan_line(f, '%d %d %d %d', 1)
    points = read_points(f, num_points)
    curves = read_entities123(f, num_curves, 1)
    surfaces = read_entities123(f, num_surfaces, 2)
    volumes = read_entities123(f, num_volumes, 3)

    entity = Entity(points, curves, surfaces, volumes)
    return entity


def read_physical_names(f):
    num_physical_names = scan_line(f, '%d')[0]

    physical_names = list[PhysicalName]()
    for i in range(num_physical_names):
        dim, physical_tag, name = scan_line(f, '%d %d %s')
        one = PhysicalName(dim, physical_tag, name)
        physical_names.append(one)

    return physical_names


def read_mesh_format(f):
    version, file_type, data_size = scan_line(f, '%f %d %d')
    mesh_format = MeshFormat(version, file_type, data_size)
    return mesh_format


def read_section_end(f):
    cc = f.read(1)
    while cc and cc != '$':
        cc = f.read(1)
    f.readline()


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

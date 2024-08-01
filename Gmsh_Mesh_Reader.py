from Element_D1O2 import Element_D1O2
from Node import Node
from Domain_D1 import Domain_D1


def read_elements(f):
    # description of all element blocks
    num_blocks, num_elements, min_element_tag, max_element_tag = scan_line(f, '%d %d %d %d')

    # check numbering
    if max_element_tag - min_element_tag + 1 != num_elements:
        raise Exception('elements not continuously numbered')

    # return values
    elements = list[G_Element]()
    element_blocks = list[Element_Block]()

    # loop over element blocks
    for i in range(num_blocks):
        # description of element block
        entity_dim, entity_tag, element_type = scan_line(f, '%d %d %d')

        # number of elements in block
        num_elements_in_block = scan_line(f, '%d')[0]

        # read elements
        element_tags = []
        for j in range(num_elements_in_block):
            tag = scan_line(f, '%d')[0]
            node_tags = list(map(int, f.readline().split()))
            one_element = G_Element(tag, node_tags)
            elements.append(one_element)

            element_tags.append(tag)  # tag starts from 1 in gmsh, convert it to zero-based

        one_block = Element_Block(entity_dim, entity_tag, element_type, element_tags)
        element_blocks.append(one_block)

    return elements, element_blocks


def read_nodes(f):
    # description of all node blocks
    num_blocks, num_nodes, min_node_tag, max_node_tag = scan_line(f, '%d %d %d %d')

    # check numbering
    if max_node_tag - min_node_tag + 1 != num_nodes:
        raise Exception('nodes not continuously numbered')

    nodes = list[G_Node]()
    node_blocks = list[Node_Block]()

    # loop over node blocks
    for i in range(num_blocks):
        # description of node block
        node_tags = []
        entity_dim, entity_tag, parametric = scan_line(f, '%d %d %d')
        num_node_in_block = scan_line(f, '%d')[0]

        if parametric != 0:
            raise Exception('parametric geometry not supported yet')

        for j in range(num_node_in_block):
            tag = scan_line(f, '%d')[0]
            node_tags.append(tag)
        for j in range(num_node_in_block):
            x, y, z = scan_line(f, '%f %f %f')
            one_node = G_Node(node_tags[j], x, y, z)
            nodes.append(one_node)

        one_block = Node_Block(entity_dim, entity_tag, parametric, node_tags)
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
    points = list[G_Point]()
    for i in range(n):
        tag, x, y, z, num_phy_tags = scan_line(f, '%d %f %f %f %d')
        phy_tags = scan_line(f, '%d ' * num_phy_tags)
        one = G_Point(tag, x, y, z, phy_tags)
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

    physical_names = list[Physical_Name]()
    for i in range(num_physical_names):
        dim, physical_tag, name = scan_line(f, '%d %d %s')
        one = Physical_Name(dim, physical_tag, name)
        physical_names.append(one)

    return physical_names


def read_mesh_format(f):
    version, file_type, data_size = scan_line(f, '%f %d %d')
    mesh_format = Mesh_Format(version, file_type, data_size)
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


class Mesh_Format(object):
    def __init__(self, version, file_type, data_size):
        self.version = version
        self.file_type = file_type
        self.data_size = data_size

    def __repr__(self):
        return f"gmsh, version: {self.version}, file_type: {self.file_type}, data_size: {self.data_size}"


class Physical_Name(object):
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


class G_Point(object):
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
    def __init__(self, points: list[G_Point], curves: list[Curve], surfaces: list[Surface], volumes: list[Volume]):
        self.points = points
        self.curves = curves
        self.surfaces = surfaces
        self.volumes = volumes

    def __repr__(self):
        return (f"Point: {len(self.points)}, Curve: {len(self.curves)}, "
                f"Surface: {len(self.surfaces)}, Volume: {len(self.volumes)}")


class G_Node(G_Point):
    def __init__(self, tag, x, y, z):
        super().__init__(tag, x, y, z, physical_tags=None)


class Node_Block(object):
    def __init__(self, dim, tag, parametric, node_tags):
        self.dim = dim
        self.tag = tag
        self.parametric = parametric
        self.node_tags = node_tags

    def __repr__(self):
        return f"dim: {self.dim}, tag: {self.tag}, nodes number: {len(self.node_tags)}"


class G_Element(object):
    def __init__(self, tag, node_tags):
        self.tag = tag
        self.node_tags = node_tags

    def __repr__(self):
        return f"(tag: {self.tag}, node tags: {self.node_tags})"


class Element_Block(object):
    def __init__(self, dim, tag, element_type, element_tags):
        self.dim = dim
        self.tag = tag
        self.element_type = element_type
        self.element_tags = element_tags

    def __repr__(self):
        return (f"dim: {self.dim}, tag: {self.tag}, "
                f"element type: {self.element_type}, elements number: {len(self.element_tags)}")


class Gmsh_Mesh_Reader(object):

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

    def read_to_region_dimension1(self) -> Domain_D1:
        nodes = list[Node]()
        elements = list[Element_D1O2]()
        for g_node in self.nodes:
            nodes.append(Node(g_node.tag - 1, g_node.x, g_node.y, g_node.z))

        boundary_elements = list[G_Element]()
        ele_idx = 0
        for g_ele in self.elements:
            if len(g_ele.node_tags) == 3:

                node0 = None
                node1 = None
                node2 = None
                for n in nodes:
                    if n.idx == g_ele.node_tags[0] - 1:
                        node0 = n
                    if n.idx == g_ele.node_tags[1] - 1:
                        node1 = n
                    if n.idx == g_ele.node_tags[2] - 1:
                        node2 = n

                if None not in [node0, node1, node2]:
                    # convert GMesh line element node order to our node order
                    # TODO: maybe use the same node order as GMesh
                    one_ele = Element_D1O2(ele_idx, node0, node2, node1)
                    elements.append(one_ele)
                    ele_idx += 1
                else:
                    raise Exception('Can not find nodes in element')

            elif len(g_ele.node_tags) == 1:
                boundary_elements.append(g_ele)

            else:
                raise Exception('Not 1-D 2nd order (3 nodes) or boundary (1 node) element')

        if len(boundary_elements) != 2:
            raise Exception('Not correct number (2) of boundary elements')

        tag0 = boundary_elements[0].tag - 1
        tag1 = boundary_elements[1].tag - 1
        for ele in elements:
            if tag0 in ele.boundary_node_idxes:
                if tag0 == ele.node0_idx:
                    ele.node0.is_boundary = True
                else:
                    ele.node2.is_boundary = True

            if tag1 in ele.boundary_node_idxes:
                if tag1 == ele.node0_idx:
                    ele.node0.is_boundary = True
                else:
                    ele.node2.is_boundary = True

        region = Domain_D1(0, elements)

        return region

    def __repr__(self):
        return (f"format: {self.mesh_format}\n"
                f"physical names: {self.physical_names}\n"
                f"entities: {self.entities}\n"
                f"nodes: {len(self.nodes)}, {self.nodes}\n"
                f"elements: {len(self.elements)}, {self.elements}\n")

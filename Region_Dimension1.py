from Node import Node
from Element_Dimension1_Order2 import Element_Dimension1_Order2 as Element_D1O2


class Region_Dimension1(object):
    def __init__(self, idx: int, elements: list[Element_D1O2]):
        self.idx = idx
        self.elements = elements
    
    def __repr__(self) -> str:
        i_str = "i:"
        x_str = "x:"
        for ele in self.elements:
            i_str += f" [{ele.node0_idx: >9}, {ele.node1_idx: >9}, {ele.node2_idx: >9}]"
            x_str += f" [{ele.node0.x: .2E}, {ele.node1.x: .2E}, {ele.node2.x: .2E}]"
            
        return i_str + '\n' + x_str
        
        
    
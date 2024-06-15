class Vertex(object):
    def __init__(self, x, idx):
        self.x = x
        self.idx = idx

    def __str__(self):
        return f"x={self.x}, idx={self.idx}"

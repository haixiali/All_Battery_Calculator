# 1-d node
class NodeDim1(object):
    def __init__(self, idx, x):
        self.idx = idx
        self.x = x

    def __str__(self):
        return f"idx: {self.idx}, x: {self.x:10.3f}"

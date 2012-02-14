""" This module contains the data representation for the structure of trees.
Note that we separate the structure of the branched tree from its
realization, that would contain things such as charges, positions, etc.
that evolve even when the structure is fixed.
"""
from random import random as uniform
from numpy import *

class Tree(object):
    def __init__(self):
        # We must carry a global (tree-level) index to access the parameter
        # data arrays
        self.n = 0
        self.segments = []
        self.root = None
        

    def add_segment(self, segment):
        """ Adds a segment to this tree.  Returns the index of the segment
        inside the tree. """
        if self.n == 0:
            self.root = segment

        self.segments.append(segment)
        index = self.n
        self.n += 1

        return index


    def parents(self, root_index=0):
        """ Builds an array with the indices to each segment parent.
        The root segment gets an index root_index. """
        p = zeros((self.n,), dtype='i')
        for i, segment in enumerate(self.segments):
            try:
                p[i] = segment.parent.index
            except AttributeError:
                p[i] = root_index

        return p

        
    def make_root(self):
        """ Creates a segment node to be root of this tree. """
        root = Segment()
        root.set_tree(self)
        self.root = root
        
        return root


    def zeros(self, dim=None):
        """ Returns an array that can hold all the data needed for a variable
        in this tree's segments.  For multi-dimension data, use dim. """
        if dim is None:
            return zeros((self.n,))

        else:
            return zeros((self.n, dim))


    def __iter__(self):
        return iter(self.segments)


class Segment(object):
    """ This is class of the segments composing a tree.  """
    def __init__(self):
        self.children = []
        self.parent = None
        self.tree = None

        
    def set_tree(self, tree):
        self.tree = tree
        self.index = tree.add_segment(self)


    def set_parent(self, parent):
        self.parent = parent
        self.set_tree(parent.tree)


    def get(self, a):
        """ Gets the value in array a corresponding to this segment. """

        return a[self.index]


    def set(self, a, value):
        """ Sets the value in array a corresponding to this segment. """

        a[self.index] = value


    def iter_adjacent(self):
        """ Iterates over all adjacent segments, including parent,
        children (if any) and brothers (if any). """

        if self.parent is not None:
            yield self.parent
            for brother in self.parent.children:
                if self != brother:
                    yield brother

        for child in self.children:
            yield child


        
    def add_child(self, other):
        other.set_parent(self)
        self.children.append(other)

    
def random_branching_tree(n, p):
    """ Builds a branched tree of n segments where every segment has a
    probability p of having two descendants.  This produces nice pictures
    and can be useful for testing. """
    tree = Tree()
    root = tree.make_root()

    leafs = [root]
    for i in xrange(n):
        l = leafs.pop(0)

        # Every leaf has at least one descendant
        s = Segment()
        l.add_child(s)
        leafs.append(s)

        # With probability p it has two children
        if random.uniform() < p:
            s = Segment()
            l.add_child(s)
            leafs.append(s)


    return tree


def sample_endpoints(tree):
    """ Gives endpoints to a tree structure.  Useful for plotting sample
    trees [DEBUG]. """

    r = tree.zeros(dim=3)
    deltav = {1: array([[0, 0, -1]]),
              2: array([[-1, 0, -1], [1, 0, -1]])}
    
    def recurse(leaf, v):
        if leaf.parent is None:
            leaf.set(r, (0, 0, 0))
        else:
            leaf.set(r, leaf.parent.get(r) + v)

        n = len(leaf.children)
        lr = leaf.get(r)
        
        for i, child in enumerate(leaf.children):
            vnew = (v * array([0.9, 1.0, 0.95]) +
                    (deltav[n][i] + random.uniform(-0.1, 0.1, size=3))
                    * exp(lr[1] / 100))
            
            recurse(child, vnew)

    recurse(tree.root, array([0, 0, 0]))

    return r


    
def test():
    import pylab
    tree = random_branching_tree(1000, 0.05)
    r = sample_endpoints(tree)

    for segment in tree:
        ep = segment.get(r)
        try:
            ip = segment.parent.get(r)
        except AttributeError:
            ip = array([0, 0, 0])

        pylab.plot([ip[0], ep[0]], [ip[2], ep[2]], lw=0.8, c='k')

    pylab.show()



if __name__ == '__main__':
    test()

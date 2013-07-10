""" This module contains the code to save/retrieve the state of a simulation.
"""

import sys
import time
from contextlib import closing
from numpy import *
import h5py

import tree

class DataFile(object):
    """ Class to store and retrieve simulation data. """
    def __init__(self, fname, parameters):
        self.fp = h5py.File(fname, 'w')
        self.main = self.fp.create_group('main')

        # We always write at least these two metadata
        self.main.attrs['command'] = ' '.join(sys.argv)
        self.main.attrs['timestamp'] = time.time()
        self.main.attrs['ctime'] = time.ctime()

        for key, item in parameters.iteritems():
            print key, item
            self.main.attrs[key] = item
        
        self.step = 0


    def add_step(self, t, tree, dist, phi, **kwargs):
        """ Adds a step to the file. """

        g = self.main.create_group('%.5d' % self.step)
        self.step += 1

        g.attrs['timestamp'] = time.time()
        g.attrs['t'] = t
        
        for field in dist._fields:
            g.create_dataset(field, data=getattr(dist, field), 
                             compression='gzip')

        g.create_dataset('phi', data=phi, compression='gzip')

        parents = tree.parents()
        g.create_dataset('parents', data=parents, compression='gzip')

        for key, item in kwargs.iteritems():
            g.create_dataset(key, data=item, compression='gzip')
            
        self.fp.flush()


    def close(self):
        self.fp.close()
    

def load_tree(file, step):
    """ Loads a tree from file and returns a tree object and the endpoints.
        file can be either an open h5file handler or a string with the filename.
    """
    try:
        # Let's assume that we received an open h5file
        parents = array(file['main/%s/parents' % step])
        r = array(file['main/%s/r' % step])
        q = array(file['main/%s/q' % step])

        tr = tree.Tree.from_parents(parents)
        return tr, tree.Distribution(r, q)

    except TypeError:
        with closing(h5py.File(file)) as fp:
            tr, dist = load_tree(fp, step)

        return tr, dist



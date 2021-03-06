""" This module implement the structure of oct-trees for the implementation
of the Fast Multipolar Method (FMM). 

. note::

  Most of this code is not used, since we disabled the use of the FMM for
  the simulations reported in the paper.  However, the simulated tree is
  embedded in a bounding box defined by a :class:`Box` instance.
"""

import timeit

from itertools import product

from numpy import *

try:
    import pylab
    from matplotlib.patches import Rectangle
except ImportError:
    pass

import mpolar

INWARD, INOUT, OUTWARD = -1, 0, 1

class Box(object):
    """ Class of 3d boxes.  Each box can be linked to a set of charges."""

    def __init__(self, r0, r1, parent=None, rel_coords=None,
                 electrode=None):
        """ Initializes a box.  r0 contains the smaller (x, y, z) coordinates
        and r1 the largest (x, y, z).
        """
        
        self.r0 = r0
        self.r1 = r1
        self.center = 0.5 * (r0 + r1)
        self.lengths = (r1 - r0)
        
        self.children = []
        self.parent = parent

        self.outward, self.inward = None, None

        if parent is not None:
            self.level = parent.level + 1
        else:
            self.level = 0

        if rel_coords is None or parent is None:
            self.coords = array([0, 0, 0])
        else:
            self.coords = parent.coords * 2 + rel_coords

        self.rf = None
        self.electrode = electrode

    def clear(self, recurse=True):
        self.outward, self.inward = None, None

        if recurse:
            for child in self.children:
                child.clear()

    
        
    def refine(self):
        """ Creates the 8 children of the box. """
        
        for t in ndindex(2, 2, 2):
            r0 = self.r0 + array(t) * self.lengths / 2
            r1 = r0 + self.lengths / 2

            self.children.append(Box(r0, r1, parent=self,
                                     rel_coords=array(t)))


    def set_charges(self, r, q, max_charges=None, min_length=0,
                    evaluation=True):
        """ Sets the charges of this box.
        If max_charges is not None, refines the box into smaller children
        until each leaf box contains no more than max_charges.
        If evaluation is true, assumes that the charge points will also
        be the evaluation points.
        """

        self.r = r
        self.q = q
        self.n = len(q)

        if evaluation:
            self.rv = self.r
            self.phi = zeros((self.r.shape[0],))
            
        # We do the reflection that implements the image charge method after
        # setting the evaluation points because we are not usually interested
        # in evaluating "image fields".
        if self.electrode:
            self.r, self.q = self.electrode.extend(self.r, self.q)
            self.n = len(self.q)

        if (max_charges is not None and self.n > max_charges
            and 2 * min_length < self.lengths[0]):
            indices = self._indices(self.r)
            
            if not self.children:
                self.refine()

            self.flt = empty([8, len(self.q)], dtype=bool)
            if evaluation:
                # Note that when electrode != None self.q and q are not the same
                # (the former includes the reflections).  To evaluate only
                # at self.rv we create a view of the first half of flt.
                self.fltv = self.flt[:, :len(q)].view()
                
            for i, child in enumerate(self.children):
                self.flt[i, :] = (indices == i)
                child.set_charges(self.r[self.flt[i, :], :],
                                  self.q[self.flt[i, :]],
                                  max_charges=max_charges,
                                  evaluation=evaluation)
                

    def update_charges(self, q):
        """ Recursively re-set the charges contained in the box.  This
        keeps the refinement oct-tree. """

        if self.electrode is not None:
            q = r_[q, self.electrode.images_q(self.r, q)]
        
        self.q[:] = q
        for i, child in enumerate(self.children):
            child.update_charges(q[self.flt[i, :]])


        
    def set_evaluation(self, rv):
        """ Recursively sets the points where the potential will be evaluated.
        """
        
        self.rv = rv
        self.phi = zeros((self.rv.shape[0],))
        indices = self._indices(rv)

        self.fltv = empty([8, len(self.phi)], dtype=bool)

        for i, child in enumerate(self.children):
            self.fltv[i, :] = (indices == i)
            child.set_evaluation(rv[self.fltv[i, :], :])


    def set_field_evaluation(self, rf):
        """ Recursively sets the points where the fields will be evaluated.
        """
        
        self.rf = rf
        self.field = zeros((self.rf.shape[0], 3))
        indices = self._indices(rf)

        self.fltf = empty([8, self.rf.shape[0]], dtype=bool)

        for i, child in enumerate(self.children):
            self.fltf[i, :] = (indices == i)
            #print self.fltf[i, :]
            if any(self.fltf[i, :]):
                child.set_field_evaluation(rf[self.fltf[i, :], :])
        
        
    def _indices(self, r):
        """ Finds the child indices (from 0 to 7) of the points at r. """
        bits = (2 * (r - self.r0[newaxis, :]) / self.lengths[newaxis, :])
        bits = bits.astype('i')

        # We do not want to exclude the boundaries with higher values
        bits =  where(bits < 2, bits, 1)

        p2 = array([4, 2, 1])

        # This is the index of the chid where each charge is sitting
        return dot(bits, p2)


    def collect_solutions(self, field=False):
        """ Collects the solution for each of the box's children. """

        for i, child in enumerate(self.children):
            child.collect_solutions(field=field)
            
            # In newer versions of numpy this gives a ValueError exception
            # here I am just ignoring it in order to maintain perfect 
            # compatibility with the results in the paper.  There phi
            # is anyhow ignored since it we only use these routines for the
            # electric field.  This bug is corrected in newer commits.
            try:
                self.phi[self.fltv[i, :]] = child.phi
            except ValueError:
                pass

            if field and child.rf is not None:
                # Note that if child contains field evaluation points,
                # then self do too.
                self.field[self.fltf[i, :], :] = child.field


    def is_near_neighbour(self, other):
        """ Checks whether other is a near neighbour of this box.
        Note that they must belong to the same oct-tree or the algorithm
        fails.
        """
        
        if self.level != other.level:
            return False

        return mpolar.are_near_neighbours(self.coords, other.coords)
        #absdif = abs(self.coords - other.coords)

        # Note that every box is considered a near-neighbour of herself.
        #return all(absdif <= 1)


    def is_well_separated(self, other):
        return not self.is_near_neighbour(other)


    def build_lists(self, recurse=False):
        """ Builds the lists of near-neighbours and the interaction list
        of this box.  Assumes that the near-neighbours are already calculated
        up in the tree. """

        self.interaction_list = []
        self.neighbours = [self]

        if self.parent is not None:
            for other in self.parent.neighbours:
                if not other.children:
                    # Here we must count the direct interaction between boxes
                    # at different levels.
                    self.neighbours.append(other)
                    
                for child in other.children:
                    if self.is_well_separated(child):
                        self.interaction_list.append(child)
                    elif self != child:
                        self.neighbours.append(child)
        
        if recurse:
            for child in self.children:
                child.build_lists(recurse=True)


    def expand(self, p):
        """ Directly calculates the multipolar expansion of this box
        around its center. """
        self.outward = mpolar.expand(p, self.r - self.center[newaxis, :],
                                    self.q, OUTWARD)
        


    def collect(self):
        """ Collects the multipolar expansions of this box's children,
        translate them to the center and sums them.  """

        for i, t in enumerate(ndindex(2, 2, 2)):
            rshift = (t - array([0.5, 0.5, 0.5])) * self.lengths / 2
            M_child = mpolar.shift(rshift, OUTWARD, self.children[i].outward)
            
            self.outward += M_child


    def upward(self, p):
        """ Goes through the oct-tree.  For leaves of the tree, directly
        calculates the multipole expansion; for nodes with descendants
        calculates the expansion by adding the children's expansion.
        This is called "Upward Pass" in the Greengard papers.
        """
        
        # Here we build an outward and an inward expansion for each box
        # that are initially set to zero
        self.inward = zeros((p, p), dtype='complex128')
        self.outward = zeros((p, p), dtype='complex128')
        
        if not self.children:
            self.expand(p)
            return

        for child in self.children:
            child.upward(p)

        self.collect()

            
    def collect_inward(self):
        """ Calculates the inward (local) expansions of all boxes in the
        interaction list and adds them. """

        for other in self.interaction_list:
            rshift = other.center - self.center

            M = mpolar.shift(rshift, INOUT, other.outward)
            
            # mpolar.accum(self.inward,
            #             mpolar.shift(rshift, INOUT, other.outward))
            
            self.inward[:, :] += M
                

    def eval_subtree(self, other):
        """ DEBUG purposes only. """
        #self.phi += mpolar.eval_array(other.outward,
        #                              self.r - other.center[:, newaxis],
        #                              OUTWARD)

        self.phi += mpolar.direct(other.r, other.q, self.r)
        
        for child in self.children:
            child.eval_subtree(other)
        
        
    def downward(self):
        """ Performs the "Downward Pass" of the Greengard papers. """

        self.collect_inward()

        for child in self.children:
            rshift = child.center - self.center
            child.inward[:, :] += mpolar.shift(-rshift, INWARD, self.inward)

            child.downward()
        

    def solve(self, a, field=False):
        """ Once we have the local expansion for the box and the list
        of near-neighbours, we can finally evaluate the potential.
        Note that generally this function is called only for leaf nodes. """

        self.phi[:] = mpolar.eval_array(self.inward,
                                        self.rv - self.center[newaxis, :],
                                        INWARD)

        for other in self.neighbours:
            self.phi += mpolar.direct(other.r, other.q, self.rv, a)


        if field and self.rf is not None:
            self.field[:, :] = mpolar.eval_field_array(
                self.inward, self.rf - self.center[newaxis, :],
                INWARD)
            for other in self.neighbours:
                self.field += mpolar.field_direct(other.r, other.q, self.rf, a)


    def solve_all(self, a=0.0, **kwargs):
        """ Calls solve for the leaf nodes of the sub-tree rooted at self. """
        if not self.children:
            self.solve(a, **kwargs)
        else:
            for child in self.children:
                child.solve_all(a, **kwargs)


    def plot(self, dims=[0, 1], recurse=False, **kwargs):
        r0 = self.r0[dims]
        r1 = self.r1[dims]
        lx, ly = r1 - r0
        
        rect = Rectangle(r0, lx, ly, **kwargs)
        pylab.gca().add_patch(rect)

        if recurse:
            for child in self.children:
                child.plot(dims=dims, recurse=True, **kwargs)


    def scatter(self, dims=[0, 1], **kwargs):
        if len(self.phi) == 0:
            return
        
        x = self.rv[:, dims[0]]
        y = self.rv[:, dims[1]]
        pylab.scatter(x, y, c=self.phi, **kwargs)
        

    def scatter_leafs(self, *args, **kwargs):
        if not self.children:
            self.scatter(*args, **kwargs)
        else:
            for child in self.children:
                child.scatter_leafs(*args, **kwargs)


    def __str__(self):
        #return str(self.coords)
        return "(%s)@%d" % (str(self.coords), self.level)


def containing_box(r, electrode=None):
    """ Builds a Box object that contains all k points in r[3, k].
    The Box has to be a perfect cube for the FMM to work.
    If reflect is True, inludes also the reflection of all points
    over the z=0 plane.
    """

    if electrode is not None:
        r = concatenate((r, electrode.images_r(r)), axis=0)
        
    rmin = amin(r, axis=0)
    rmax = amax(r, axis=0)
    
    lengths = rmax - rmin
    center = 0.5 * (rmax + rmin)
    sides = amax(lengths) * ones((3,))

    r0 = center - sides / 2
    r1 = center + sides / 2
    
    return Box(r0, r1, electrode=electrode)



def main():
    k = 1200
    
    r = random.uniform(-1, 1, size=(k, 3))
    q = random.uniform(-1.0, 1.0, size=k)

    # Let's make things simpler
    #r[2, :] = pi / 5
    q[:] = 1.0
    
    r0 = array([-1.0, -1.0, -1.0])
    r1 = array([1.0, 1.0, 1.0])


    #pylab.plot(r[0, :], r[1, :], 'o', mfc='k', mec='k')
    
    box = Box(r0, r1)
    box.plot(recurse=True, fill=False)

    pylab.xlim([-1, 1])
    pylab.ylim([-1, 1])

    box.set_charges(r, q, max_charges=5)
    box.build_lists(recurse=True)
    box.upward(15)
    box.downward()
    box.solve_all()
    box.collect_solutions()

    phi = mpolar.direct(r, q, r, 0.0)

    box.scatter_leafs(vmin=0, vmax=1200)
    pylab.colorbar()

    # Let's compare with the exact solution
    pylab.figure(2)
    box.plot(recurse=True, fill=False)

    pylab.xlim([-1, 1])
    pylab.ylim([-1, 1])


    pylab.scatter(r[:, 0], r[:, 1],
                  c=phi, vmin=0, vmax=1200)
    pylab.colorbar()

    err = sqrt(sum((phi - box.phi)**2)) / k
    savetxt("cmp.txt", c_[phi, box.phi])
    
    print "Error = %g" % err
    
    pylab.show()



if __name__ == '__main__':
    main()
    
    

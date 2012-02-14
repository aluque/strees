""" Testing the 1d version of the FMM. """

from numpy import *
import pylab

c = 1

def phi(x):
    return sqrt(x*x + c*c)

def s(x, y, d):
    """ The exact function.
    * x is x[n]
    * y is y[m]
    * d is d[m]
    """
    #return sum(d * phi(y - x))
    phi_ = phi(y[:, newaxis] - x[newaxis, :])
    
    return dot(d, phi_)

def expansion_term(t):
    z = zeros(t.shape)
    q = c_[1 + z,
           -t,
           c**2 / 2. + z,]
    #       (c**2 * t)/2.,
    #       -c**4 / 8. + (c**2 * t**2) / 2.,
    #       (-3 * c**4 * t + 4 * c**2 * t**3)/8.,
    #       (c**6 - 12 * c**4 * t**2 + 8 * c**2 * t**4)/16.]

    return q

def expansion(y0, y, d):
    """ Returns the coefficients of an expansion around y0.
     * y0 is a scalar
     * y is a vector y[m]
     * d is a vector d[m]

     * Returns a vector a[p], where p is the degree of the expansion
     """

    # This matrix q should have shape q[m, p]
    q = expansion_term(y - y0)
    
    # Remember that d is d[m] and q is q[m, p].  dot(a, b) runs over
    # the last axis of a and the one before last of b.
    a = dot(d, q)

    return a


def r(x, y0, a):
    """ Applies the expansion a, cetered at y0 to points in x.
    * x is a vector x[n]
    * y0 is a scalar
    * a is a vector of coefficients a[p]

    * Returns a vector r[n].
    """
    p = a.shape[0]
    Q = column_stack([(x - y0)**(1-i) for i in range(p)])
    
    # a is a[p] but Q is Q[n, p] so we have to use Q.T
    return dot(a, Q.T)

# Implementation of the binary tree
class Segment(object):
    def __init__(self, x0, x1, level=0, max_level=3):
        self.x0 = x0
        self.x1 = x1
        self.diameter = x1 - x0
        self.level = level
        self.mid = (x0 + x1) / 2.0
        self.max_level = max_level
        self.children = []
        
        if self.level < self.max_level:
            self.children.append(Segment(self.x0, self.mid,
                                        level=self.level + 1,
                                        max_level=max_level))

            self.children.append(Segment(self.mid, self.x1,
                                        level=self.level + 1,
                                        max_level=max_level))

    def filt(self, x):
        return logical_and(x > self.x0, x <= self.x1)

    def contains(self, x):
        return (x > self.x0) and (x <= self.x1)
    
    def calculate_expansions(self, y, d):
        flt = self.filt(y)
        self.y, self.d = y[flt], d[flt]
        
        self.a = expansion(self.mid, self.y, self.d)
        
        for child in self.children:
            child.calculate_expansions(self.y, self.d)


    def interacting(self, other):
        """ Recursively finds all descendant segments of self
        that are interacting with other. """

        if self.is_well_separated(other):
            yield self
            return
        
        if not self.children:
            yield self
            return
        
        for child in self.children:
            for segment in child.interacting(other):
                yield segment
                
    def contribute(self, other, x, f):
        print "Interacting %s <-> %s" % (str(self), str(other))
        flt = other.filt(x)
        
        if self.is_well_separated(other):
            # Use the expansion
            f[flt] += sign(other.mid - self.mid) * r(x[flt], self.mid, self.a)
            return

        # Calculate the full interaction.
        f[flt] += s(x[flt], self.y, self.d)
        print x[flt]


    def is_well_separated(self, other):
        return abs(other.mid - self.mid) > (2 * self.diameter)
    
    
    def smallest_containing(self, x):
        if not self.contains(x):
            return None

        if not self.children:
            return self

        for child in self.children:
            s = child.smallest_containing(x)
            if s:
                return s

    def leafs(self):
        """ Returns an iterator over the tree's leafs. """
        if not self.children:
            yield self

        for child in self.children:
            for l in child.leafs():
                yield l

    def calculate_full(self, x):
        f = zeros(shape(x))
        for leaf in self.leafs():
            for source in self.interacting(leaf):
                source.contribute(leaf, x, f)

        return f
    
            
    def __str__(self):
        return "(%f, %f]@%d" % (self.x0, self.x1, self.level)

    

def main():
    L = 100.0
    y = random.uniform(0, L, size=50)
    d = 0.5 - random.uniform(size=50)
    
    a = expansion(0.5, y, d)

    x = linspace(0, 120, 1000)

    root = Segment(0.0, L, max_level=6)
    root.calculate_expansions(y, d)
    leaf = root.smallest_containing(0.7128)

    print leaf
    print [str(T) for T in root.leafs()]
    print [str(T) for T in root.interacting(leaf)]
    
    #pylab.plot(x, r(x, 0.5, a), lw=1.5, c='b')
    pylab.plot(x, s(x, y, d), lw=1.5, c='k')
    pylab.plot(x, root.calculate_full(x), lw=1.5, c='r')
    
    pylab.plot(y, d, 'o')
    
    pylab.show()


if __name__ == '__main__':
    main()

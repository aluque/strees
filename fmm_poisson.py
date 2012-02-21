import sys
import timeit

from numpy import *
import pylab

from refinement import Box, containing_box
import mpolar

MAX_CHARGES_PER_BOX = 200
MULTIPOLAR_TERMS = 4
CONDUCTOR_THICKNESS = 0.01
EXTERNAL_FIELD = array([0.0, 0.0, -5.0])
TIP_MOBILITY = 1.0
END_TIME = 1.0
TIME_STEP = 1e-4

def main():
    d = loadtxt(sys.argv[1])
    r = d[:, 0:3].T
    q = d[:, 3]
    phi0 = d[:,4]
    phi1 = d[:,5]
    
    # r = random.uniform(-1, 1, size=(3, k))
    # q = random.uniform(-1.0, 1.0, size=k)

    box = containing_box(r)
    box.set_charges(r, q, max_charges=MAX_CHARGES_PER_BOX)
    box.build_lists(recurse=True)
    box.update_charges(q)
    box.upward(MULTIPOLAR_TERMS)
    box.downward()
    box.solve_all(a=CONDUCTOR_THICKNESS, field=True)
    
    box.collect_solutions(field=True)

    phi = mpolar.direct(r, q, r, CONDUCTOR_THICKNESS)
    err = sum((phi - box.phi)**2) / len(phi)
    print "FMM error = %g" % err
    print "max(q) = ", max(q)
    
    pylab.plot(r[2, :], phi0, 'o', mew=0, ms=2.5, label="FMM.orig")
    pylab.plot(r[2, :], phi1, 'o', mew=0, ms=2.5, label="direct.orig")
    pylab.plot(r[2, :], box.phi, 'o', mew=0, ms=2.5, label="FMM.new")
    pylab.plot(r[2, :], phi, 'o', mew=0, ms=2.5, label="direct.new")
    pylab.legend()
    pylab.show()
    
    # box.scatter(dims=[0, 2], s=5.0, faceted=False)
    # pylab.show()


if __name__ == '__main__':
    main()


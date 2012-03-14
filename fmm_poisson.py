import sys
import timeit

from numpy import *
import pylab
import h5py

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
    # d = loadtxt(sys.argv[1])
    # r = d[:, 0:3].T
    # q = d[:, 3]
    # phi0 = d[:,4]
    # phi1 = d[:,5]
    
    fname, step = sys.argv[1:3]
    
    fp = h5py.File(fname, "r")
    g = fp['main']
    r = array(g[step]['r']).T
    q = array(g[step]['q'])
    phi0 = array(g[step]['phi'])

    CONDUCTOR_THICKNESS = fp['main'].attrs['conductor_thickness']
    MAX_CHARGES_PER_BOX = fp['main'].attrs['max_charges_per_box']
    MULTIPOLAR_TERMS = fp['main'].attrs['multipolar_terms']

    print "%d charges" % len(q)
    
    # r = random.uniform(-1, 1, size=(3, k))
    # q = random.uniform(-1.0, 1.0, size=k)

    box = containing_box(r)
    box.set_charges(r, q,
                    max_charges=MAX_CHARGES_PER_BOX,
                    min_length=2 * CONDUCTOR_THICKNESS)
    
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
    print r.shape, phi0.shape, q.shape
    
    # pylab.plot(r[2, :], phi0, 'o', mew=0, ms=2.5, label="FMM.orig")
    # pylab.plot(r[2, :], phi1, 'o', mew=0, ms=2.5, label="direct.orig")
    # pylab.plot(r[2, :], box.phi, 'o', mew=0, ms=2.5, label="FMM.new")
    # pylab.plot(r[2, :], phi, 'o', mew=0, ms=2.5, label="direct.new")

    pylab.figure(figsize=(7, 14))
    pylab.subplot(2, 1, 1)
    pylab.plot(r[2, :], phi - box.phi, 'o', mew=0, ms=2.5, label="error")
    pylab.legend()

    pylab.subplot(2, 1, 2)
    pylab.plot(r[2, :], phi, 'o', mew=0, ms=2.5, label="direct")
    pylab.plot(r[2, :], box.phi, 'o', mew=0, ms=2.5, label="fmm")
    pylab.legend()
    
    pylab.show()
    
    # box.scatter(dims=[0, 2], s=5.0, faceted=False)
    # pylab.show()


if __name__ == '__main__':
    main()


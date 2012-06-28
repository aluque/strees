""" Calculation of the errors in solveing the potential and dq/dt. """

from numpy import *
from numpy.linalg import norm
import h5py
import pylab

from refinement import Box, containing_box
import mpolar
import datafile
from plotter import plot_projections, bounding_box

def main():
    from optparse import OptionParser
    from contextlib import closing

    parser = OptionParser()

    (opts, args) = parser.parse_args()

    steps = args[1:]
    fname = args[0]

    for step in steps:
        with closing(h5py.File(fname)) as fp:
            tr, r = datafile.load_tree(fp, step)
            q = array(fp['main'][step]['q'])
            M = tr.ohm_matrix(r)

            r = r.T
            
            CONDUCTOR_THICKNESS = fp['main'].attrs['conductor_thickness']
            MAX_CHARGES_PER_BOX = fp['main'].attrs['max_charges_per_box']
            MULTIPOLAR_TERMS = fp['main'].attrs['multipolar_terms']
            EXTERNAL_FIELD = fp['main'].attrs['external_field']

            EXTERNAL_FIELD_VECTOR = array([0.0, 0.0, EXTERNAL_FIELD])

            box = containing_box(r, reflect=True)
            box.set_charges(r, q,
                            max_charges=MAX_CHARGES_PER_BOX,
                            min_length=16 * CONDUCTOR_THICKNESS)
            
            box.build_lists(recurse=True)
            box.update_charges(q)
            box.upward(MULTIPOLAR_TERMS)
            box.downward()
            box.solve_all(a=CONDUCTOR_THICKNESS, field=True)
            
            box.collect_solutions(field=True)
                
            phi = mpolar.direct(r, q, r, CONDUCTOR_THICKNESS)


            dq0 = M.dot(phi - dot(r.T, EXTERNAL_FIELD_VECTOR))
            dq1 = M.dot(box.phi - dot(r.T, EXTERNAL_FIELD_VECTOR))

            pylab.figure(figsize=(7, 14))
            pylab.subplot(2, 1, 1)
            pylab.plot(abs(dq0), abs(dq0 - dq1), 'o', mew=0, ms=2.5,
                       label="abs error")
            pylab.loglog()
            pylab.legend()

            pylab.subplot(2, 1, 2)
            pylab.plot(abs(dq0), abs((dq0 - dq1) / dq0), 'o', mew=0,
                       ms=2.5, label="rel error")
            pylab.loglog()
            pylab.legend()

            pylab.figure(figsize=(12, 10))
            r0, r1 = bounding_box(r)

            plot_projections(r.T, abs(dq0 - dq1), r0, r1, log=False)
            pylab.show()

if __name__ == '__main__':
    main()


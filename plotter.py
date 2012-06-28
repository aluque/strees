import sys
import os, os.path
from optparse import OptionParser
from matplotlib.colors import LogNorm

from numpy import *
import h5py
import datafile

try:
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.colors
except ImportError:
    pass


class Plotter(object):
    def __init__(self, fname):
        self.fp = h5py.File(fname, "r")
        self.main = self.fp['main']
        self.run_name = self.main.attrs['run_name']
        self.external_field = self.main.attrs['external_field']
        self.external_field_vector = array([0.0, 0.0, self.external_field])

        self.ref_step = None

        
    @property
    def steps(self):
        return self.main.keys()


    def __getitem__(self, key):
        return self.main[key]


    def set_ref(self, step):
        self.ref_step = step
        
        r = array(self.main[step]['r']).T
        self.r0, self.r1 = bounding_box(r)

        q = array(self.main[step]['q'])
        # Focus only on the positive charges
        aqmax = amax(q)
        self.qmin, self.qmax = -aqmax, aqmax
        #self.qmin, self.qmax = amin(q), amax(q)
        
        
    def plot(self, step):
        r = array(self.main[step]['r']).T
        q = array(self.main[step]['q'])
        if self.ref_step is None:
            self.set_ref(step)
        
        plot_projections(r.T, q, self.r0, self.r1,
                         vmin=self.qmin, vmax=self.qmax, plot3d=True)


    def plot_field(self, step):
        tr, r = datafile.load_tree(self.fp, step)
        phi = array(self.main[step]['phi'])
        
        p = tr.parents()
        l = tr.lengths(r)
        midpoints = tr.midpoints(r)
        
        # As we save the now (Tue Jun 12 17:09:04 2012) the latest
        # points do not have a calculated potential.
        n = len(phi)
        p, l, r = p[:n], l[:n], r[:n]
        midpoints = midpoints[:n]
        
        phi = phi - dot(r, self.external_field_vector)
        
        efields = (phi - phi[p]) / l
        try:
            plot_projections(midpoints, abs(efields), self.r0, self.r1,
                             vmin=None, vmax=None, plot3d=True)
        except ValueError:
            pass
        
        
    
def main():
    parser = OptionParser()

    parser.add_option("--ref", dest="ref", type="str",
                      help="The reference step", default=None)

    parser.add_option("--show", dest="show", action="store_true",
                      help="Open the matplotlib window?", default=False)

    parser.add_option("--field", dest="field", action="store_true",
                      help="Plot the electric field instead of the charge?",
                      default=False)

    parser.add_option("--print-parameters", dest="print_parameters",
                      action="store_true",
                      help="The reference step", default=None)

    parser.add_option("--print-times", dest="print_times",
                      action="store_true",
                      help="Print real and simulated times for each step",
                      default=False)

    (opts, args) = parser.parse_args()

    fname = args[0]
    plotter = Plotter(fname)
    
    try:
        os.mkdir(plotter.run_name)
    except OSError:
        pass
    
    steps = args[1:]
    if not steps:
        steps = plotter.steps
        
    if opts.ref is not None:
        plotter.set_ref(opts.ref)
    else:
        plotter.set_ref(steps[-1])
    

    print "%s [%s] (%d steps)" % (plotter.run_name,
                                  plotter.main.attrs['ctime'], len(steps))
    

    if opts.print_parameters:
        for key, item in plotter.main.attrs.iteritems():
            print "%-30s =\t%s" % (key, repr(item))

    if opts.print_times:
        for i, step in enumerate(steps):
            t = plotter.main[step].attrs['t']
            timestamp = plotter.main[step].attrs['timestamp']
            print "%s\t%f\t%f" % (step, t, timestamp)


    pylab.figure(figsize=(14, 10))
    for i, step in enumerate(steps):
        print ("[%s]" % step),
        sys.stdout.flush()

        if not ((i +1) % 10):
            print ''
        
        if not opts.field:
            plotter.plot(step)
        else:
            plotter.plot_field(step)
        
        if opts.show:
            pylab.show()
        
        pylab.savefig(os.path.join(plotter.run_name,
                                   '%s_%s.png' % (plotter.run_name, step)))



def bounding_box(r):
    rmin = amin(r, axis=1)
    rmax = amax(r, axis=1)
    
    lengths = rmax - rmin
    center = 0.5 * (rmax + rmin)
    sides = amax(lengths) * ones((3,))

    r0 = center - sides / 2
    r1 = center + sides / 2

    return r0, r1


def plot_projections(r, q, r0, r1, vmin=None, vmax=None, log=False,
                     plot3d=False):
    X, Y, Z = 0, 1, 2
    names = ["X", "Y", "Z"]
    
    axes = [(X, Z), (Y, Z), (X, Y)]
    pylab.subplots_adjust(left=0.1, wspace=0.35, hspace=0.2,
                          right=0.85, top=0.95)
    cmap = pylab.get_cmap("jet")
    #cmap = charge_cmap()
    
    for i, (d1, d2) in enumerate(axes):
        # For gray use axisbg='#eeefef'
        ax = pylab.subplot(2, 2, i, axisbg='#404060')
        ax.clear()
        ax.grid(ls='-', lw=1.0, c='#707090', zorder=-20)
        norm = None if not log else LogNorm()
        
        pylab.scatter(r[:, d1], r[:, d2], c=q,
                      s=9.5, faceted=False, vmin=vmin, vmax=vmax,
                      cmap=cmap, zorder=20, norm=norm),
        #pylab.colorbar()
        
        ax.set_xlabel(names[d1])
        ax.set_ylabel(names[d2])

        ax.set_xlim([r0[d1], r1[d1]])
        ax.set_ylim([r0[d2], r1[d2]])

    ax = pylab.axes([0.90, 0.1, 0.025, 0.85])
    pylab.colorbar(cax=ax)

    if plot3d:
        ax = pylab.subplot(2, 2, 3, projection='3d')
        ax.scatter(r[:, 0], r[:, 1], r[:, 2], zdir='z', c=q,
                   s=9.5, faceted=False, vmin=vmin, vmax=vmax,
                   cmap=cmap, zorder=20, norm=norm)

        ax.set_xlim3d([r0[X], r1[X]])
        ax.set_ylim3d([r0[Y], r1[Y]])
        ax.set_zlim3d([r0[Z], r1[Z]])



def charge_cmap():
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
             'green': ((0.0, 0.0, 0.0),
                       (0.5, 0.9, 0.9),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.0, 1.0),
                      (0.5, 0.5, 0.5),
                      (1.0, 0.0, 0.0))}

    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 256)
    return my_cmap


if __name__ == '__main__':
    main()
    

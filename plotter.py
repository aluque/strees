import sys
import os, os.path
from optparse import OptionParser
from matplotlib.colors import LogNorm
import scipy.constants as co

from numpy import *
import h5py
import datafile

import cmaps

try:
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.colors
except ImportError:
    pass


kV_cm = co.kilo / co.centi

class Plotter(object):
    def __init__(self, fname, **kwargs):
        self.fp = h5py.File(fname, "r")
        self.main = self.fp['main']
        self.run_name = self.main.attrs['run_name']
        self.external_field = self.main.attrs['external_field']
        self.external_field_vector = array([0.0, 0.0, self.external_field])
        self.ref_step = None

        self.branches = 0
        
        options = dict(single=False, axisbg='#c7c0c0')
        options.update(kwargs)
        
        for key, value in options.iteritems():
            setattr(self, key, value)

        print self.single
        
    @property
    def steps(self):
        return self.main.keys()


    def __getitem__(self, key):
        return self.main[key]


    def set_ref(self, step):
        self.ref_step = step
        r, q = self.charge_dens(step)
        
        self.r0, self.r1 = bounding_box(r)

        # Focus only on the positive charges
        aqmax = amax(q)
        #self.qmin, self.qmax = -aqmax, aqmax
        self.qmin, self.qmax = amin(q), amax(q)
        

    def charge_dens(self, step):
        tr, r = datafile.load_tree(self.fp, step)
        q = array(self.main[step]['q'])
        l = tr.lengths(r)

        # We drop the charges at the endpoints of each branch 
        t = tr.terminals()
        n = len(t)
        q, r, l = q[:-n], r[:-n], l[:-n]
        
        q = q / l
        q = where(isfinite(q), q, 0)
        return r, q / (co.nano / co.centi)
    
    def plot(self, step):
        if self.ref_step is None:
            self.set_ref(step)

        r, q = self.charge_dens(step)
        
        plot_projections(r / co.centi, q,
                         self.r0 / co.centi, self.r1 / co.centi,
                         vmin=self.qmin, vmax=self.qmax, plot3d=True,
                         reduce_range=5, dynamic=True,
                         label="Linear charge density [nC/cm]",
                         single=self.single,
                         axisbg=self.axisbg)
        # with reduce_range=35, we see the sign everywhere.

        savetxt("charge.dat", q)

        
    def plot_field(self, step):
        tr, r = datafile.load_tree(self.fp, step)
        phi = array(self.main[step]['phi'])
        
        p = tr.parents()
        l = tr.lengths(r)
        midpoints = tr.midpoints(r)
        t = tr.terminals()
        self.branches = len(t)
        # print "%d branches." % len(t)
        # As we save now (Tue Jun 12 17:09:04 2012) the latest
        # points do not have a calculated potential.
        n = len(phi)

        p, l, r = p[:n], l[:n], r[:n]
        midpoints = midpoints[:n]
        
        print self.external_field_vector
        phi = phi - dot(r, self.external_field_vector)
        
        efields = (phi - phi[p]) / l
        
        try:
            plot_projections(midpoints / co.centi,
                             -efields / kV_cm,
                             self.r0 / co.centi, self.r1 / co.centi,
                             vmin=None, vmax=None, plot3d=True, log=False,
                             dynamic=True,
                             label="Electric field, $E$ [kV/cm]",
                             single=self.single,
                             axisbg=self.axisbg)
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

    parser.add_option("--format", dest="format", type="str",
                      help="Format of the output figures", default='png')

    parser.add_option("--single", dest="single", action="store_true",
                      help="Plot only one projection", default=False)

    parser.add_option("--axisbg", dest="axisbg", type="str",
                      help="Background color", default='#eaeaea')

    (opts, args) = parser.parse_args()

    fname = args[0]
    plotter = Plotter(fname, single=opts.single, axisbg=opts.axisbg)
    
    try:
        os.mkdir(plotter.run_name)
    except OSError:
        pass
    
    steps = args[1:]
    if not steps:
        steps = plotter.steps
        
    steps = [s if s != 'last' else plotter.steps[-1] for s in steps]
    
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


    pylab.figure(figsize=(13, 10.5))
    for i, step in enumerate(steps):
        
        if not opts.field:
            plotter.plot(step)
        else:
            plotter.plot_field(step)
        
        print ("[%s (%d)]" % (step, plotter.branches)),
        sys.stdout.flush()

        if not ((i + 1) % 10):
            print ''

        if opts.show:
            pylab.show()
        
        pylab.savefig(os.path.join(plotter.run_name,
                                   '%s_%s.%s' % (plotter.run_name, step,
                                                 opts.format)))



def bounding_box(r):
    rmin = amin(r, axis=0)
    rmax = amax(r, axis=0)
    
    lengths = rmax - rmin
    center = 0.5 * (rmax + rmin)
    sides = amax(lengths) * ones((3,))

    r0 = center - sides / 2
    r1 = center + sides / 2

    return r0, r1


def plot_projections(r, q, r0, r1, vmin=None, vmax=None, log=False,
                     plot3d=False, dynamic=False, label=None,
                     reduce_range=None, single=False,
                     axisbg='#404060'):
    X, Y, Z = 0, 1, 2
    names = ["X [cm]", "Y [cm]", "Z [cm]"]
    
    axes = [(X, Z), (Y, Z), (X, Y)]
    if single:
        axes = [(Y, Z)]
        
    pylab.subplots_adjust(left=0.1, wspace=0.35, hspace=0.2,
                          right=0.85, top=0.95)
    cmap = pylab.get_cmap("jet")
    #cmap = charge_cmap()
    extend = 'neither'
    
    if vmin is None or vmax is None:
            vmin = nanmin(q)
            vmax = nanmax(q)

    if reduce_range is not None:
        vmin, vmax = vmin, vmax / reduce_range
        extend = 'max'
        
    if dynamic:
        cmap = cmaps.get_colormap('bluered', dynamic=True)
        cmap.center = -vmin / (vmax - vmin)
    

    iplot = [1, 2, 3]
    for i, (d1, d2) in enumerate(axes):
        # For gray use axisbg='#eeefef'
        if not single:
            ax = pylab.subplot(2, 2, iplot[i], axisbg=axisbg) # was #404060
        else:
            ax = pylab.subplot(1, 1, 1, axisbg=axisbg) # was #404060
            
        ax.clear()
        ax.grid(ls='-', lw=1.0, c='#c0c0c0', zorder=-20)
        # Thu Aug 23 15:58:10 2012
        # I used this and --axisbg='aaaaaa' for the reconnection plot:
        #ax.grid(ls='-', lw=1.0, c='#888888', zorder=-20)
        norm = None if not log else LogNorm()
        
        pylab.scatter(r[:, d1], r[:, d2], c=q,
                      s=19.5, faceted=False, vmin=vmin, vmax=vmax,
                      cmap=cmap, zorder=20, norm=norm),
        #pylab.colorbar()
        
        ax.set_xlabel(names[d1])
        ax.set_ylabel(names[d2])

        ax.set_xlim([r0[d1], r1[d1]])
        ax.set_ylim([r0[d2], r1[d2]])

    ax = pylab.axes([0.90, 0.1, 0.025, 0.85])
    cbar = pylab.colorbar(cax=ax, extend=extend)
    if label is not None:
        cbar.set_label(label)
        
    if plot3d and not single:
        ax = pylab.subplot(2, 2, 4, projection='3d')
        ax.scatter(r[::-1, 0],
                   r[::-1, 1],
                   r[::-1, 2], zdir='z', c=q[::-1],
                   s=9.5, faceted=False, vmin=vmin, vmax=vmax,
                   cmap=cmap, zorder=20, norm=norm)

        #ax.set_xlim3d([nanmin(r[:, X]), nanmax(r[:, X])])
        #ax.set_ylim3d([nanmin(r[:, Y]), nanmax(r[:, Y])])
        #ax.set_zlim3d([nanmin(r[:, Z]), nanmax(r[:, Z])])
        
        ax.set_xlim3d([r0[X], r1[X]])
        ax.set_xlabel("X [cm]")
        ax.set_ylim3d([r0[Y], r1[Y]])
        ax.set_ylabel("Y [cm]")
        ax.set_zlim3d([r0[Z], r1[Z]])
        ax.set_zlabel("Z [cm]")



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
    

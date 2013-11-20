import sys
import os, os.path
from optparse import OptionParser
from matplotlib.colors import LogNorm
import matplotlib as mpl
import scipy.constants as co
from scipy.stats import scoreatpercentile
from numpy import *
import h5py
import datafile

import cmaps
import tree

try:
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.colors
except ImportError:
    pass

AIR_N_STP = co.value('Loschmidt constant (273.15 K, 101.325 kPa)')
Td = 1e-21

class Variable(object):
    def __init__(self, func, doc=None, name="", units="?"):
        self.func = func
        self.name = name
        self.units = units
        self.doc = doc

    def __call__(self, sim):
        r = self.func(sim)
        # Apply the filter only if not a scalar
        return r
        

VARIABLES = {}
def variable(**kwargs):
    """ A decorator to define variables from functions. """
    def deco(f):
        kwargs.update({'doc': f.__doc__})
        VARIABLES[f.func_name] = Variable(f, **kwargs)
        return f
    return deco


class DataContainer(object):
    def __init__(self, fname):
        self.fp = h5py.File(fname, "r")
        self.main = self.fp['main']
        
    def __getitem__(self, key):
        return self.main.attrs[key]

    def load_step(self, step):
        self.tr, self.dist = datafile.load_tree(self.fp, step)        
        
        self.dist.r[:, 2] += self['ionosphere_height']

        self.r = self.dist.r
        self.q = self.dist.q

        self.phi = array(self.main[step]['phi'])
        self.t = self.main[step].attrs['t']

    @property
    def steps(self):
        return self.main.keys()



class Plot(object):
    """ This is an abstract class to handle all kinds of plots."""
    AXES_NAMES = ['$x$', '$y$', '$z$']
    LABEL_FORMAT = "{axes} [{units}]"

    def __init__(self, axes=None, makefig=True, dynamic=False, 
                 truncate=False, log=False, lscale=1, lunits='m'):
        if axes is None:
            axes = pylab.gca()

        self.axes = axes
        self.makefig = makefig
        self.dynamic = dynamic
        self.vmin, self.vmax = None, None
        self.cmap = pylab.get_cmap("jet")
        self.norm = None if not log else LogNorm()
        self.lscale = lscale
        self.lunits = lunits
        self.truncate = truncate

    def set_limits(self, vmin=None, vmax=None, v=None, r=None):
        if not vmin is None and not vmax is None and v is None:
            self.vmin, self.vmax = vmin, vmax
        elif vmin is None and vmax is None and v is not None:
            if not self.truncate:
                self.vmin, self.vmax = nanmin(v), nanmax(v)
            else:
                self.vmin = -15 * scoreatpercentile(-v[v < 0], 50.)
                self.vmax =  5  * scoreatpercentile( v[v > 0], 50.)

        else:
            raise ValueError("set_limit called with inconsistent parameters")

        if self.dynamic:
            self.cmap = cmaps.get_colormap('bluered', dynamic=True)
            self.cmap.center = -self.vmin / (self.vmax - self.vmin)
            
        if r is not None:
            self.r0, self.r1 = bounding_box(r)
            self.r0 /= self.lscale
            self.r1 /= self.lscale

    def plot(self):
        raise NotImplementedError("You must subclass Plot")


    def finish(self):
        raise NotImplementedError("You must subclass Plot")



class Plot2D(Plot):
    """ A single 2D projection. """
    def __init__(self, projection, **kwargs):
        super(Plot2D, self).__init__(**kwargs)
        indices = {'x':0, 'y':1, 'z':2}
        iprojection = [indices[a.lower()] for a in projection]
        
        self.xaxis, self.yaxis, self.zaxis = iprojection
    

    def plot(self, r, v):
        self.axes.clear()
        self.axes.grid(ls='-', lw=1.0, c='#c0c0c0', zorder=-20)
        isort = argsort(r[:, self.zaxis])

        if self.vmin is None:
            self.set_limits(v=v, r=r)

        self.mappable = self.axes.scatter(r[isort, self.xaxis] / self.lscale, 
                                          r[isort, self.yaxis] / self.lscale, 
                                          c=v[isort],
                                          s=14.5, faceted=False,
                                          vmin=self.vmin, vmax=self.vmax,
                                          cmap=self.cmap, zorder=20, 
                                          norm=self.norm)

        xlabel = self.LABEL_FORMAT.format(axes=self.AXES_NAMES[self.xaxis],
                                          units=self.lunits)
        ylabel = self.LABEL_FORMAT.format(axes=self.AXES_NAMES[self.yaxis],
                                          units=self.lunits)

        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)

        self.axes.set_xlim([self.r0[self.xaxis], self.r1[self.xaxis]])
        self.axes.set_ylim([self.r0[self.yaxis], self.r1[self.yaxis]])

    def finish(self):
        self.axes.colorbar(self.mappable)


class Plot3D(Plot):
    """ A single 3D projection. """
    def __init__(self, **kwargs):
        super(Plot3D, self).__init__(**kwargs)

    def plot(self, r, v):
        self.axes.clear()
        #self.axes.grid(ls='-', lw=1.0, c='#c00000', zorder=0)
        isort = argsort(dot(array([-1, 1, 0]), r.T))

        if self.vmin is None:
            self.set_limits(v=v, r=r)

        self.axes.scatter(r[isort, 0] / self.lscale, 
                          r[isort, 1] / self.lscale, 
                          r[isort, 2] / self.lscale, 
                          c=v[isort],
                          s=9.5, zdir='z', faceted=False,
                          vmin=self.vmin, vmax=self.vmax,
                          cmap=self.cmap, norm=self.norm)
        
        xlabel = self.LABEL_FORMAT.format(axes=self.AXES_NAMES[0],
                                          units=self.lunits)
        ylabel = self.LABEL_FORMAT.format(axes=self.AXES_NAMES[1],
                                          units=self.lunits)
        zlabel = self.LABEL_FORMAT.format(axes=self.AXES_NAMES[2],
                                          units=self.lunits)
        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)
        self.axes.set_zlabel(zlabel)

        self.axes.set_xlim3d([self.r0[0], self.r1[0]])
        self.axes.set_ylim3d([self.r0[1], self.r1[1]])
        self.axes.set_zlim3d([self.r0[2], self.r1[2]])

    def finish(self):
        self.axes.colorbar(self.mappable, 
                           extend='both' if self.truncate else 'neither')


class CombinedPlot(Plot):
    """ Combined 2D projections """
    def __init__(self, **kwargs):
        super(CombinedPlot, self).__init__(**kwargs)
        axisbg='#eaeaea'

        pylab.subplots_adjust(left=0.1, wspace=0.35, hspace=0.2,
                              right=0.85, top=0.95)

        self.axes = [pylab.subplot(2, 2, i + 1, axisbg=axisbg) for i
                     in xrange(3)]
        self.axes.append(pylab.subplot(2, 2, 4, axisbg=axisbg, projection='3d'))

        self.subplots = [Plot2D('xzy', axes=self.axes[0], **kwargs),
                         Plot2D('yzx', axes=self.axes[1], **kwargs),
                         Plot2D('xyz', axes=self.axes[2], **kwargs),
                         Plot3D(axes=self.axes[3], **kwargs)]
        self.cax = pylab.axes([0.88, 0.1, 0.025, 0.85])

    def set_limits(self, *args, **kwargs):
        for p in self.subplots:
            p.set_limits(*args, **kwargs)
        
    def plot(self, r, v):
        for p in self.subplots:
            p.plot(r, v)
        
        self.mappable = self.subplots[0].mappable

    def finish(self, time=None):
        self.cbar = pylab.colorbar(self.mappable, 
                                   cax=self.cax,
                                   extend='both' if self.truncate 
                                   else 'neither')

        if time is not None:
            p, f = prefix_and_factor(time)
            pylab.figtext(0.975, 0.025, "t = %.2g %ss" % (time / f, p),
                          color="#883333",
                          ha='right', va='bottom', size='x-large')

def prefix_and_factor(x):
    prefixes = {
        -12: 'p',
        -9: 'n',
        -6: 'u',
        -3: 'm',
         0: '',
         3: 'k'}

    n = 3 * round(log10(x) / 3)
    if n in prefixes:
        return prefixes[n], 10**n
    
    return '', 1


def iter_steps(s):
    """ Iterate over timesteps according to a string that can be
    e.g. '00034..47' or '00087'. """
    
    try:
        fro, to = s.split('..')
    except ValueError:
        yield s
        return

    l = len(fro)
    if ',' in to:
        to, step = to.split(',')
    else:
        step = 1

    ifro, ito, istep = int(fro), int(to), int(step)
    for i in xrange(ifro, ito, istep):
        yield ("%.*d" % (l, i))


def bounding_box(r, expand_r0=None, expand_r1=None):
    rmin = amin(r, axis=0)
    rmax = amax(r, axis=0)
    
    lengths = rmax - rmin
    center = 0.5 * (rmax + rmin)
    sides = amax(lengths) * ones((3,))

    if expand_r0 is None:
        expand_r0 = array([1, 1, 1])

    if expand_r1 is None:
        expand_r1 = array([1, 1, 1])

    r0 = center - expand_r0 * sides / 2
    r1 = center + expand_r1 * sides / 2

    return r0, r1


UNITS = {
    'mm': co.milli,
    'cm': co.centi,
    'm': 1,
    'km': co.kilo}

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input file")
    parser.add_argument("steps", 
                        help="Step/steps to plot (empty for all)",
                        default=None)
    parser.add_argument("var", 
                        default="field",
                        help="The variable to plot")
    parser.add_argument("--show", action='store_true',
                        help="Show the resulting plot",
                        default=False)
    parser.add_argument("--units", action='store',
                        help="Show the resulting plot",
                        choices=UNITS.keys(),
                        default='cm')
    parser.add_argument("--ref", action='store',
                        help="Reference step by default, the latest one",
                        default=None)
    parser.add_argument("--log", action='store_true',
                        help="Use a logarithmic color scale",
                        default=False)
    parser.add_argument("--truncate", action='store_true',
                        help="Truncate the color scale",
                        default=False)
    parser.add_argument("-d", "--outdir",
                        help="Output directory (may contain {rid} and {var})", 
                        action='store', default='{rid}')
    parser.add_argument("-o", "--output",
                        help="Output file (may contain {rid} {step} and {var})", 
                        action='store', default='{rid}_{step}_{var}.png')

    args = parser.parse_args()

    var = VARIABLES[args.var]
    datacontainer = DataContainer(args.input)

    rid = os.path.splitext(args.input)[0]
    outdir = args.outdir.format(rid=rid, var=args.var)
    try:
        os.mkdir(outdir)
    except OSError:
        pass

    if args.steps is None or args.steps == 'all':
        steps = datacontainer.steps
    elif args.steps == 'latest':
        steps = [datacontainer.steps[-1]]
    else:
        steps = list(iter_steps(args.steps))

    if args.ref is None:
        ref = steps[-1]
    else:
        ref = args.ref

    pylab.figure(figsize=(13, 10.5))
    plot = CombinedPlot(log=args.log, dynamic=True, 
                        lunits=args.units, lscale=UNITS[args.units],
                        truncate=args.truncate)

    if len(steps) > 1 and args.show:
        print "For your own good, --show is incompatible with more than 1 plot."
        sys.exit(-1)

    # Set the reference
    datacontainer.load_step(ref)
    r, v = var(datacontainer)
    plot.set_limits(v=v, r=r)

    for step in steps:
        datacontainer.load_step(step)
        r, v = var(datacontainer)
        plot.plot(r, v)
        plot.finish(time=datacontainer.t)
        if args.show:
            pylab.show()
        else:
            ofile = args.output.format(step=step, var=args.var, rid=rid)
            ofile = os.path.join(outdir, ofile)

            pylab.savefig(ofile)
            #pylab.close()
            print "File '%s' saved" % ofile



@variable(name="$E$", units="V/m")
def field(data):
    p = data.tr.parents()
    l = data.tr.lengths(data.dist)
    midpoints = data.tr.midpoints(data.dist)
    t = data.tr.terminals()
    n = len(t)

    p, l, r, phi = p[:-n], l[:-n], data.r[:-n], data.phi[:-n]
    midpoints = midpoints[:-n]

    efields = (phi - phi[p]) / l

    return midpoints, -efields


@variable(name="$E/n$", units="Td")
def en(data):
    p = data.tr.parents()
    l = data.tr.lengths(data.dist)
    midpoints = data.tr.midpoints(data.dist)
    t = data.tr.terminals()
    n = len(t)

    p, l, r, phi = p[:-n], l[:-n], data.r[:-n], data.phi[:-n]
    midpoints = midpoints[:-n]

    n = AIR_N_STP * exp(-midpoints[:, 2] / data['scale_height'])
    
    efields = (phi - phi[p]) / l
    en = efields / n / Td
    return midpoints, -en


@variable(name="$\lambda$", units="C/m")
def charge(data):
    p = data.tr.parents()
    l = data.tr.lengths(data.dist)
    midpoints = data.tr.midpoints(data.dist)
    t = data.tr.terminals()
    n = len(t)

    q, r, l = data.q[:-n], data.r[:-n], l[:-n]
    midpoints = midpoints[:-n]

    q = q / l
    q = where(isfinite(q), q, 0)
    return r, q


@variable(name="$\phi$", units="V")
def phi(data):
    return data.r, data.phi


@variable(name="$\sigma$", units="m/$\Omega$")
def sigma(data):
    midpoints = data.tr.midpoints(data.dist)
    return midpoints, data.dist.s


@variable(name="$E_\perp$", units="V/m")
def eperp(data):
    p = data.tr.parents()
    l = data.tr.lengths(data.dist)
    qsegment = 0.5 * (data.q + data.q[p])
    lmbd = qsegment / l
    midpoints = data.tr.midpoints(data.dist)

    er = (2 * data['maxwell_factor'] * lmbd / data.dist.a)
    t = data.tr.terminals()
    n = len(t)

    return midpoints[:-n], er[:-n]


if __name__ == '__main__':
    main()

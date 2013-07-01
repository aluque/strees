from matplotlib.colors import Colormap
import matplotlib.cm as cm
import matplotlib.colors as colors

from numpy import where

def get_colormap(cmap, dynamic=False):
    """ If cmap is not None, sets the colormap by te specified string"""

    cmaps = {'hot': cm.hot,
             'autumn': cm.autumn,
             'bone': cm.bone,
             'cool': cm.cool,
             'copper': cm.copper,
             'flag': cm.flag,
             'gray': cm.gray,
             'hsv': cm.hsv,
             'jet': cm.jet,
             'pink': cm.pink,
             'prism': cm.prism,
             'spring': cm.spring,
             'summer': cm.summer,
             'winter': cm.winter,
             'YlOrRd': cm.YlOrRd,

             'invjet': invjet,
             'invhot': invhot,
             'blue': blue,
             'hottest': hottest,
             'invhottest': InvertedColormap(hottest),
             'darkhot': darkhot,
             'invdarkhot': InvertedColormap(darkhot),
             'darkred': darkred,
             'invdarkred': InvertedColormap(darkred),
             'coolhot': coolhot,
             'reddish': reddish,
             'bluered': bluered,
             #'spectra': spectra
             }

    if cmap:
        try:
            cmap = cmaps[cmap]
        except KeyError:
            warn("Colormap `%s' unknown: using default colormap" % cmap)
            return None

    if dynamic:
        return DynamicColormap(cmap)
    else:
        return cmap
    
class InvertedColormap(Colormap):
    def __init__(self, cm):
        self._cm = cm

    def __getattr__(self, name):
        return getattr(self._cm, name)
            
    def __call__(self, X, alpha=1.0, **kwargs):
        return self._cm(1.0 - X, alpha=alpha, **kwargs)


class DynamicColormap(Colormap):
    def __init__(self, basis, center=0.5):
        self.center = center
        self.basis = basis
        self.__dict__.update(self.basis.__dict__)
        
    #def __getattr__(self, name):
    #    return getattr(self.basis, name)

    #def __setattr__(self, name, value):
    #    setattr(self.basis, name, value)
        
    def __call__(self, X, **kwargs):
        #if self.center <= 0.5:
        #    Y = 0.5 * (1 + (X - self.center) / (1 - self.center))
        #else:
        #    Y = 0.5 * X / self.center
        Y = where(X < self.center,
                  0.5 * X / self.center,
                  0.5 * (1 + (X - self.center) / (1 - self.center)))
        
        #if X <= self.center:
        #    Y = 0.5 * X / self.center
        #else:
        #    Y = 0.5 * (1 + (X - self.center) / (1 - self.center))

        return self.basis(Y, **kwargs)


def combine_data(low, high):
    newd = {}
    for color in 'blue', 'red', 'green':
        l = [(0.5 * (1.0 - x), y, z) for x, y, z in low[color]]
        l.reverse()
        h = [(0.5 * (1.0 + x), y, z) for x, y, z in high[color]]
        newd[color] = tuple(l + h)

    return newd

def swap_data(data, col1, col2):
    newd = data.copy()
    newd[col1] = newd[col2]
    newd[col2] = data[col1]

    return newd



# Some colormaps definitions

_blue_data = {'blue': ((0., 0.0416, 0.0416),
                       (0.365079, 1.000000, 1.000000),
                       (1.0, 1.0, 1.0)),
             'green': ((0., 0., 0.),
                       (0.365079, 0.000000, 0.000000),
                       (0.746032, 1.000000, 1.000000),
                       (1.0, 1.0, 1.0)),
             'red':  ((0., 0., 0.),
                      (0.746032, 0.000000, 0.000000),
                      (1.0, 1.0, 1.0))}                  


_darkhot_data = {'red':   ((0., 0.0416, 0.0416),(0.565079, 1.000000, 1.000000),(1.0, 1.0, 1.0)),
             'green': ((0., 0., 0.),(0.565079, 0.000000, 0.000000),
                       (0.9, 1.000000, 1.000000),(1.0, 1.0, 1.0)),
             'blue':  ((0., 0., 0.),(0.9, 0.000000, 0.000000),(1.0, 1.0, 1.0))}                  

_hottest_data = {'red': ((0., 0.0416, 0.0416),
                         (0.365079, 1.000000, 1.000000),
                         (1.0, 1.0, 1.0)),
             'green': ((0., 0., 0.),
                       (0.646032, 0.200000, 0.200000),
                       (0.846032, 0.400000, 0.400000),
                       (1.0, 1.0, 1.0)),
             'blue':  ((0., 0., 0.),
                       (0.446032, 0.000000, 0.000000),
                       (1.0, 1.0, 1.0))}                  

_darkred_data = {'red': ((0., 0.6, 0.6), (1.0, 0.6, 0.6)),
                 'green': ((0., 0., 0.), (1.0, 1.0, 0.)),
                 'blue': ((0., 0., 0.), (1.0, 1.0, 0.))}
                 

_reddish_data = {'red': ((0., 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.6, 0.6)),
                 'green': ((0., 1.0, 1.0), (0.6, 0.6, 0.6), (1.0, 0., 0.)),
                 'blue': ((0., 1.0, 1.0), (0.4, 0., 0.), (1.0, 0, 0))}
                 
_redhot_data = {'red': ((0., 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.6, 0.6)),
                 'green': ((0., 1.0, 1.0), (0.6, 0.6, 0.6), (1.0, 0., 0.)),
                 'blue': ((0., 1.0, 1.0), (0.4, 0., 0.), (1.0, 0, 0))}


_greyred_data = {'red': ((0., 1.0, 1.0), (0.25, 1.0, 1.0), (1.0, 0.5, 0.5)),
                 'green': ((0., 0.9, 0.9), (0.25, 0.2, 0.2), (1.0, 0., 0.)),
                 'blue': ((0., 0.4, 0.4), (0.25, 0.0, 0.0), (1.0, 0, 0))}

_bluegrey_data = {'red': ((0., 1.0, 1.0), (0.15, 0.0, 0.0), (1.0, 0., 0.)),
                  'green': ((0., 0.9, 0.9), (0.33, 0.0, 0.0), (1.0, 0., 0.)),
                  'blue': ((0., 0.4, 0.4), (0.33, 1.0, 1.0), (1.0, 0.5, 0.5))}

_bluish_data = swap_data(_reddish_data, 'red', 'blue')
_coolhot_data = combine_data(_bluish_data, _reddish_data)
_bluered_data = combine_data(_bluegrey_data, _greyred_data)

invjet = InvertedColormap(cm.jet)
invhot = InvertedColormap(cm.hot)
blue   = colors.LinearSegmentedColormap('blue', _blue_data, cm.LUTSIZE)
hottest = colors.LinearSegmentedColormap('hottest', _hottest_data, cm.LUTSIZE)
darkhot = colors.LinearSegmentedColormap('darkhot', _darkhot_data, cm.LUTSIZE)
darkred = colors.LinearSegmentedColormap('darkred', _darkred_data, cm.LUTSIZE)
coolhot = colors.LinearSegmentedColormap('coolhot', _coolhot_data, cm.LUTSIZE)
reddish = colors.LinearSegmentedColormap('reddish', _reddish_data, cm.LUTSIZE)
bluered = colors.LinearSegmentedColormap('bluered', _bluered_data, cm.LUTSIZE)


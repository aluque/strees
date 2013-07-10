""" This is an auxiliary module for reading input .ini files.
It is based on Python's stdlib 
`ConfigParser <http://docs.python.org/2/library/configparser.html/>`_

This module also provides functionality for setting run sets where
one or more of the parameters run over a list of values.
"""

import os, os.path, socket
from ConfigParser import SafeConfigParser, RawConfigParser, NoOptionError
from warnings import warn
from itertools import product
import re

import numpy

def guess_type(s):
    """ Converts s into int or float if possible.  If not, leave it as a
    string. This will give us problems if the user wants to e.g. use
    filenames or run names that can be parser as a number.  So have to
    implement a better approach whereby parameter names are associated with
    types.
    """

    try:
        return int(s)
    except ValueError:
        pass

    try:
        return float(s)
    except ValueError:
        pass

    return s


# These are decorators to check allowed values.
def positive(func):
    """ A decorator to constraint the parameter to positive values. """
    def f(s):
        r = func(s)
        if not r >= 0:
            raise ValueError("%s must be positive" % func.func_name)
        return r
    
    f.__doc__ = func.__doc__
    f.func_name = func.func_name

    return f


def nonnegative(func):
    """ A decorator to constraint the parameter to nonnegative values. """
    def f(s):
        r = func(s)
        if r < 0:
            raise ValueError("%s must be positive" % func.func_name)
        return r
    
    f.__doc__ = func.__doc__
    f.func_name = func.func_name

    return f

def boolean(func):
    """ A decorator to check that the parameter can be interpreted as a
    boolean. """
    choices = {'true', 'false'}
    def f(s):
        if not s.lower() in choices:
            raise ValueError("%s must be a boolean" % func.func_name)

        r = func(s)
        return r
    
    f.__doc__ = func.__doc__
    f.func_name = func.func_name

    return f

# A decorator to set a default value, stored in the dictionary PARAM_DEFAULTS
PARAM_DEFAULTS = {}
def default(value):
    def deco(func):
        global PARAM_DEFAULTS
        PARAM_DEFAULTS[func.func_name] = value
        return func

    return deco


def load_input(fname, parameters, d=None, upper=False, raw=False):
    """ Loads an input file and stores its values in the dictionary
    d.   If upper is true, transforms the parameter names to upper case. """

    if d is None:
        d = dict()
    
    config = SafeConfigParser()
    defaults = dict(home=os.environ['HOME'],
                    user=os.environ['LOGNAME'],
                    cwd=os.getcwd(),
                    input=os.path.splitext(os.path.basename(fname))[0],
                    input_dir=os.path.split(os.path.realpath(fname))[0],
                    hostname=socket.gethostname())
                    
    config.read(fname)
    
    for section in ['global', 'parameters']:
        for name, value in config.items(section, vars=defaults, raw=raw):
            print name, repr(value)
            try:
                func = getattr(parameters, name)
                value = expand(value, func)
                print "%-30s =\t%-10s [%s]" % (name, repr(value), func.__doc__)
                
            except AttributeError:
                warn("'%s' is not defined as a parameter." % name)
                value = guess_type(value)

            if upper:
                name = name.upper()

            d[name] = value

    r = PARAM_DEFAULTS.copy()
    r.update(d)
    return r


RE_LIST = re.compile(r'@\((.+)\)')
RE_LOOP = re.compile(r'@@\(([\d.-]+):([\d.-]+)(:([\d.-]+))?\)')
def expand(s, parser):
    """ Expands special characters to produce e.g. lists of many parameters.
    """
    # All expansions start with the symbol '@':
    if not '@' in s:
        return parser(s)

    m = RE_LIST.match(s)
    if m:
        return [parser(x) for x in m.group(1).split()]

    m = RE_LOOP.match(s)
    if m:
        if m.group(4) is not None:
            a, b, c = parser(m.group(1)), parser(m.group(2)), parser(m.group(4))
        else:
            a, b, c = parser(m.group(1)), parser(m.group(2)), None            
        r = numpy.arange(a, b, c)
        print r
        print [parser(x) for x in r]
        return [parser(str(x)) for x in r]



def expand_dict(d):
    """ Takes a dictionary that may contain a few lists as values and returns
    an iterator over dictionaries where each of these elements is iterated. """

    keys, lists = zip(*[(k, v) for k, v in d.iteritems()
                        if isinstance(v, list)])
    for tpl in product(*lists):
        d0 = d.copy()
        for k, v in zip(keys, tpl):
            d0[k] = v

        yield d0
        
    
def expand_input(ifile, parameters):
    """ Expands an input file by expanding some of its arguments.
    """
    d = load_input(ifile, parameters, d={}, upper=False, raw=True)
    base = os.path.splitext(ifile)[0]
    
    config = RawConfigParser()
    config.read(ifile)
    l = []

    for i, d0 in enumerate(expand_dict(d)):
        fname = '%s_%.4d.ini' % (base, i)
        l.append(fname)
        with open(fname, 'w') as fp:
            for k, v in d0.iteritems():
                for sect in config.sections():
                    try:
                        old = config.get(sect, k)
                        config.set(sect, k, v)
                    except NoOptionError:
                        pass

            config.write(fp)
    
    return l


def main():
    import sys
    import parameters
    
    expand_input(sys.argv[1], parameters)

    # d = load_input(sys.argv[1], parameters, d={}, upper=True)
    # for d0 in expand_dict(d):
    #     print d0
    

if __name__ == '__main__':
    main()
    

#! /usr/bin/python

import distutils.sysconfig
from distutils.core import setup, Extension
import numpy as np

numod = Extension('mpolar',
                  sources = ['mpolarmod.c'],
                  extra_objects = ['misc.o', 'multipol.o', 'efield.o'],
                  include_dirs = [np.get_include(),
                                  distutils.sysconfig.get_python_inc()]
                  )

setup (ext_modules = [numod],
       name = 'mpolar',
       version = '1.0',
       description = 'Multipolar expansions',
       author = 'Alejandro Luque',
       author_email = 'aluque@iaa.es')

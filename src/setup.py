#! /usr/bin/python

from distutils.core import setup, Extension

numod = Extension('mpolar',
                  sources = ['mpolarmod.c'],
                  extra_objects = ['misc.o', 'multipol.o', 'efield.o'],
                  include_dirs = ['/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/',
                                  '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7/']
                  )

setup (ext_modules = [numod],
       name = 'mpolar',
       version = '1.0',
       description = 'Multipolar expansions',
       author = 'Alejandro Luque',
       author_email = 'aluque@iaa.es')

User guide
==========

Installation
------------

Required Python libraries
^^^^^^^^^^^^^^^^^^^^^^^^^

``strees`` uses some Python libraries for the numerical computations as well
as its output format and plotting.  You have to have them installed and
accesible to run ``strees``:

  * NumPy and SciPy.  See `here <http://scipy.org/install.html>`_ for installation instructions for your operating system.
  * `Matplotlib <http://matplotlib.org/>`_
  * `h5py <http://code.google.com/p/h5py/>`_

For Linux and Mac OS X we recommend using the package repositories of your
distribution or MacPorts/Homebrew for Macs.

Compilation
^^^^^^^^^^^

``strees`` is mostly written in Python but in order to speed up the code, 
C is used in the calculation of electrostatic interactions.  For this
purpose we built the *mpolar* library, which implements both 
Fast Multipole Method (FMM) calculations and direct, :math:`O(N^2)` 
electrostatic computations.

This means that you have to compile the C code to obtain a ``mpolar`` 
library.  You can do this by typing::

  make

from the ``/src`` sub-directory.

Before you run the code, make sure that the ``mpolar`` library is accessible
to the dynamic linker of your operating system.  In Linux you achive that with::

  export LD_LIBRARY_PATH=/path/to/your/mpolar.so

In Mac OS X, you do the same with::

  export DYLD_LIBRARY_PATH=/path/to/your/mpolar.so



Running a simulation
--------------------

Once the code is set up, you use the script grow_tree.py to start a 
simulation::

  python /path/to/grow_tree.py simulation.ini

where the file ``simulation.ini`` contains the input parameters of the 
simulation.  See :doc:`parameters` for a list of all valid parameters with a
short description.



Input files
-----------

Input files for the code are ``.ini`` files with key/value pairs.  Here is
an example::

  [global]
  run_name = %(input)s
  out_file = %(input_dir)s/%(run_name)s.h5
  desc = A canonical simulation with the standard parameter set

  [parameters]
  external_field = -1.5e6
  tip_mobility = 0.09
  end_time = 3.5e-07
  time_step = 2.5e-10
  branching_probability = 100.0
  branching_sigma = 0.0001
  conductor_thickness = 0.001
  conductance = 9.6e-07
  maxwell_factor = 9e9
  fmm_threshold = 500000
  max_charges_per_box = 200
  multipolar_terms = 10
  random_seed = 11
  electrode_geometry = planar

The input file contains two sections: 

  * ``[global]`` contains the parameters ``run_name``, ``out_file`` and 
    ``desc``.  In this section you can use substitutions.  The code
    provides substitutions variables 
      * ``%(input)``, containing the *basename* of the input file (i.e.
        if the input file is ``simulation.ini`` ``%(input)`` is replaced
        by ``simulation``.
      * ``%(input_dir)`` is replaced by the directory containing the input
        file.

  * ``[parameters]`` contains all other parameters.  In general you can
    specify these parameters in any unit system.  In the example above
    we used SI units, so we set ``maxwel_factor=9e9``.  This is 
    :math:`1/4\pi \varepsilon_0`.


Output files
------------

The code produces a single output file per simulation, specified
in the :func:`parameters.out_file` parameter.  The output is a 
`HDF5 <http://www.hdfgroup.org/HDF5/>`_ file, organized as follows:

The file contains a single group, called ``main``.  This group contains
attributes containing the input parameters used in the simulation as well
as information about the running environment:

``command``
   The command that was used to start the simulation.

``timestamp``
   The machine time when the simulation was started.

``ctime``
   Human-readable version of ``timestamp``

``user``
   Login name of the user that run the simulation.

``hostname``
   The name of the computer where the simulation run,


The group ``main`` contains
a sub-group per saved timestep, with names ``00000``, ``00001`` and so forth.
Each of these sub-groups contains a snapshot of the tree at the given time.
The attributes of the groups are:

``t``
   Simulated time of the snapshot.

``timestamp``
   Machine time when the snapshot was saved


The data is contained in these fields:

``r``
   An :math:`N \times 3` array containing the locations of all nodes 
   :math:`i, i=0, 1, \dots N-1`.

``q``
   The charge of each node.

``phi``
   The electrostatic potential of each node.

``parents``
   An array containing the index of the parent node of each node.  The parent
   node is the immediate one further up the tree.

``error``
   An estimation of the numerical error in the potentials :math:`\phi`.

``error_dq``
   An estimation of the numerical error in :math:`dq_i/dt`


Plotting
--------

The script ``plotter.py`` can be used to plot the results of a 
simulation for a quick inspection.  Use it with::

   python /path/to/plotter.py outputfile.h5 [steps] [options]

Where

  * ``step`` indicates the step that you want to plot, such as ``00300``.
    You can also specify ``latest`` to see the latest snapshot.  If you
    do not provide a step, the program will plot *all* the available steps.

  * ``options`` stands for a combination of these options:

       ``-h, --help``         
         Show a help message and exit.
       ``--ref=REF``
         The reference step
       ``--show``
         Open the matplotlib window?
       ``--field``
         Plot the electric field instead of the charge?
       ``--print-parameters``
         The reference step
       ``--print-times``
         Print real and simulated times for each step
       ``--format=FORMAT``
         Format of the output figures
       ``--single``
         Plot only one projection
       ``--axisbg=AXISBG``
         Background color

If you do not set the ``--show`` option, the code will produce an output
file named ``outputfile/outputfile_XXXXX.png``, where *XXXXX* is the step
index.


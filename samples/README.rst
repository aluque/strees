Sample parameter files
======================

This directory contains some sample parameter files, including all the files 
for the simulations that have been included in the paper.  The files are organized in the following sub-directories:

``sample_long``
  A set of 25 simulations with the standard parameters but changing the random 
  seed.  The one used for figure 3, 4, and 6 in the paper is 
  ``sample_long_0011.ini``.  Figures 9 and 10, showing a reconnection event
  are based on ``sample_long_0022.ini``.

``conductivity``
  Simulations to study the effect of the channel conductivity :math:`\sigma`.  
  They were used to generate figure 7.

``sigma``
  Simulations to investigate the effect of the initial separation 
  between branches :math:`\ell_{sib}`.  These were used to generate
  figure 8.

``reconnect``
  Simulations to investigate reconnection with different values of the
  average branching length and different random seeds.   Figure 11 was 
  produced using these simulations.

``convergence``
  These simulations were used too test the convergence of the code with
  decreasing timesteps.  They were used for figure A1.

``emin``
  These are simulations corresponding to the modification of the model 
  suggested in appendix B.  Figure B1 is based on ``sample_long_0017.ini``.



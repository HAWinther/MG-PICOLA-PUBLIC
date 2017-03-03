# MG-PICOLA-PUBLIC
MG-PICOLA v1.0, March 2017
Author: Hans A. Winther

The code extends the [COLA (COmoving Lagrangian Acceleration)](https://arxiv.org/abs/1301.0322) method for simulating cosmological structure formation from LCDM to theories that exhibit scale-dependent growth at the level of 2LPT. The code includes the fast approximate screening method of [Winther & Ferreira (2014)](https://arxiv.org/abs/1403.6492). For the scientific paper explaining the approach see [Winther, Koyama, Manera, Wright and Zhao (2017)](https://arxiv.org/abs/1703.00879).

This code is based on the [L-PICOLA](https://github.com/CullanHowlett/l-picola) code written by Cullan Howlett & Marc Manera. For a documentation of L-PICOLA see [L-PICOLA/Documentation](https://github.com/CullanHowlett/l-picola/tree/master/Documentation)

 - Requires the [FFTW3 library](http://www.fftw.org/download.html) (needs to be compiled with --enable-float to use the SINGLE_PRECISION option).

 - Requires [GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/)

 - Compile the code as make -f Makefile.model and run as mpirun -np 1 MG_PICOLA_MODEL paramfile.txt. See [paramfiles](paramfiles) for some example parameter-files.

 - To implement a new model see [src_v3/user_defined_functions.h](src_v3/user_defined_functions.h). The code has f(R), DGP and general (m(a), beta(a)) models included.

 - A simple code to extract power-spectra from output-files (GADGET / ASCII) can be found in SimplePofk. Will add this and a halo-finder to the code at some point. See [MatchMaker](https://github.com/damonge/MatchMaker) for a FoF halo-finder that can be run on the GADGET output-files.

 - The scale-dependent version needs the define -DSCALEDEPENDENT. This version requires several Fourier transforms per time-step which makes the code ~5 times slower.

 - The lightcone version of the code has not been tested. Neither has the non-gaussian IC generation.

 - Some optimizations should be done in the SCALEDEPENDENT versions with respect to Output(). Currently we recompute the LPT fields twice here. Should be changed.

MG-PICOLA is distributed under the GNU Public License v3 (see COPYING for details).

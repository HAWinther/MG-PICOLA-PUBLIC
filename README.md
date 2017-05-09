# MG-PICOLA-PUBLIC
MG-PICOLA v1.1, May 2017
Author: Hans A. Winther

The code extends the [COLA (COmoving Lagrangian Acceleration)](https://arxiv.org/abs/1301.0322) method for simulating cosmological structure formation from LCDM to theories that exhibit scale-dependent growth at the level of 2LPT. The code includes the fast approximate screening method of [Winther & Ferreira (2014)](https://arxiv.org/abs/1403.6492). For the scientific paper explaining the approach see [Winther, Koyama, Manera, Wright and Zhao (2017)](https://arxiv.org/abs/1703.00879).

This code is based on the [L-PICOLA](https://github.com/CullanHowlett/l-picola) code written by Cullan Howlett & Marc Manera. For a documentation of L-PICOLA see [L-PICOLA/Documentation](https://github.com/CullanHowlett/l-picola/tree/master/Documentation)
 
 - See [doc/note.pdf](doc/note.pdf) for some (preliminary) documentation.

 - Requires the [FFTW3 library](http://www.fftw.org/download.html) (needs to be compiled with --enable-float to use the SINGLE\_PRECISION option).

 - Requires [GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/)

 - Compile the code as make -f Makefile.model and run as mpirun -np 1 MG\_PICOLA\_MODEL paramfile.txt. See [paramfiles](paramfiles) for some example parameter-files.

 - To implement a new model see [src\_v3/user\_defined\_functions.h](src_v3/user_defined_functions.h). The code have f(R), DGP and general (m(a), beta(a)) models included.

 - If the parameter input\_pofk\_is\_for\_lcdm = 1 then the code assumes the linear P(k) we read in (or from fittin formula) is for LCDM and uses the ratio of growth-factor in MG and LCDM to compute the MG power-spectrum used to set the IC. If also input\_sigma8\_is\_for\_lcdm = 1 then the normalization of sigma8 is done for LCDM; thi is useful to make MG and LCDM runs with the same IC.

 - A simple code to extract power-spectra from output-files (GADGET / ASCII) can be found in SimplePofk.

 - The code can now also compute the power-spectrum at every time-step by compiling with the option COMPUTE\_POFK. This is activated by adding pofk\_compute\_ever\y\_step = 1 to the parameterfile. Other options one must add is : pofk\_nbins (number of bins), pofk\_bintype (0 linear, 1 logarithmic), pofk\_subtract\_shotnoise (1 do it, 0 don't), pofk\_kmin (h/Mpc), pofk\_kmax  (h/Mpc). Put all of these to negative values to use fiducial values (nbins = Nmesh, linear, subtract shotnoise, kmin = 2 * pi / Box, kmax = 2 * pi/Box * Nmesh).

 - The code can now also compute RSD power-spectrum multipoles (P0, P2 and P4). Add the extra parameter pofk\_compute\_rsd\_pofk = 1 to the the options above. See doc/ for more information.
 
 - The code has been linked to [MatchMaker](https://github.com/damonge/MatchMaker) which allows us to do halo-finding on the fly. See doc/ for how to use it. Not properly tested yet so use with care.

 - The scale-dependent version needs the define SCALEDEPENDENT. This version requires several Fourier transforms per time-step which makes the code ~5 times slower.

 - The lightcone version of the code have not been tested, but should work fine for the scale-independent version of the code (i.e. using LCDM growth-factors).

 - Some optimizations should be done in the SCALEDEPENDENT versions with respect to Output(). Currently we recompute the LPT fields twice here. Should be changed.

MG-PICOLA is distributed under the GNU Public License v3 (see COPYING for details).

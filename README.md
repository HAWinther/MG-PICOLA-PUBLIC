# MG-PICOLA-PUBLIC
MG-PICOLA v0.9, May 2017
Author: Hans A. Winther

The code extends the [COLA (COmoving Lagrangian Acceleration)](https://arxiv.org/abs/1301.0322) method for simulating cosmological structure formation from LCDM to theories that exhibit scale-dependent growth at the level of 2LPT. The code includes the fast approximate screening method of [Winther & Ferreira (2014)](https://arxiv.org/abs/1403.6492). For the scientific paper explaining the approach see [Winther, Koyama, Manera, Wright and Zhao (2017)](https://arxiv.org/abs/1703.00879).

This code is based on the [L-PICOLA](https://github.com/CullanHowlett/l-picola) code written by Cullan Howlett & Marc Manera. For a documentation of L-PICOLA see [L-PICOLA/Documentation](https://github.com/CullanHowlett/l-picola/tree/master/Documentation)
 
 - See the [doc/note.pdf](doc/note.pdf) for some (preliminary) documentation. This code is under development and features are currently being added. A stable validated version will come soon.

 - Requires the [FFTW3 library](http://www.fftw.org/download.html) (needs to be compiled with --enable-float to use the SINGLE\_PRECISION option).

 - Requires [GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/)

 - Compile the code as make Makefile MODEL=FOFR (DGP, BRANSDICKE, MBETA, etc.) and run as mpirun -np 1 MG\_PICOLA\_MODEL paramfile.txt. See [paramfiles](paramfiles) for some example parameter-files and the [documentation](doc/note.pdf) for more information.

 - To implement a new model see [src\_v3/user\_defined\_functions.h](src_v3/user_defined_functions.h). The code have LCDM, LCDM with (wa,w0) parametrization, Jordan-Brans-Dicke, f(R), nDGP and general (m(a), beta(a)) models included so far + massive neutrinos.

 - If the parameter input\_pofk\_is\_for\_lcdm = 1 then the code assumes the linear P(k) we read in (or from fitting formula) is for LCDM and uses the ratio of growth-factor in MG and LCDM to compute the MG power-spectrum used to set the IC. If also input\_sigma8\_is\_for\_lcdm = 1 then the normalization of sigma8 is done for LCDM; thi is useful to make MG and LCDM runs with the same IC.

 - A very simple code to extract power-spectra from output-files (RAMSES / GADGET / ASCII) can be found in SimplePofk (separate version for OpenMP and MPI). Includes a code to estimate redshift space multipoles using the global parallel plane approximation.

 - The code can now also compute the power-spectrum at every time-step by compiling with the option COMPUTE\_POFK. This is activated by adding pofk\_compute\_ever\y\_step = 1 to the parameterfile. Other options one must add is : pofk\_nbins (number of bins), pofk\_bintype (0 linear, 1 logarithmic), pofk\_subtract\_shotnoise (1 do it, 0 don't), pofk\_kmin (h/Mpc), pofk\_kmax  (h/Mpc). Put all of these to negative values to use fiducial values (nbins = Nmesh, linear, subtract shotnoise, kmin = 2 * pi / Box, kmax = 2 * pi/Box * Nmesh).

 - The code can now also compute RSD power-spectrum multipoles (P0, P2 and P4). Add the extra parameter pofk\_compute\_rsd\_pofk = 1 or 2 to the the options above. See the [documentation](doc/note.pdf) for more information. Not properly tested yet!
 
 - The code has been linked to [MatchMaker](https://github.com/damonge/MatchMaker) which allows us to do halo-finding on the fly. See the [documentation](doc/note.pdf) for how to use it. Not properly tested yet!

 - The scale-dependent version needs the define SCALEDEPENDENT. This version requires several Fourier transforms per time-step which makes the code ~3-4 times slower.

 - The lightcone version of the code have not been tested, but should work fine for the scale-independent version of the code (i.e. using LCDM growth-factors).

 - Added memory monitoring (using wrappers around malloc and free).

MG-PICOLA is distributed under the GNU Public License v3 (see COPYING for details).

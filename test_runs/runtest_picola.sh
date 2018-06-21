#!/bin/bash

#=============================================
# This is for running a simple
# test of the code for comparing LCDM sims of
# MGPICOLA to that of PICOLA.
#=============================================

#==================================
# Paths and options
#==================================
picolaoutputdir="output_picola"
mymgpicolaexec_dgp="MG_PICOLA_DGP"
mymgpicolaexec_fofr="MG_PICOLA_FOFR"
myepspdf="epspdf"
plotname="pofk_picola_mgpicola_testrun"
recompile="true"
makeplot="true"
runsim="true"
ncpu="4"

#==================================
# Function to generate parameterfile
#==================================
function make_parameter_file_dgp(){
  FileBase="$1"
  OutputDir="$2"

  paramfile="
modified_gravity_active         0
rcH0_DGP                        1.0
Rsmooth                         1.0
include_screening               1 
use_lcdm_growth_factors         1
input_pofk_is_for_lcdm          1
input_sigma8_is_for_lcdm        1
amplitude_fixed_initial_condition 0
inverted_initial_condition        0
OutputDir                       $OutputDir
FileBase                        $FileBase
OutputRedshiftFile              data_picola_comparison/output_redshifts_n20.dat
NumFilesWrittenInParallel       1  
UseCOLA                         1           
Buffer                          1.5 
Nmesh                           128      
Nsample                         128      
Box                             512.0   
Init_Redshift                   29.0    
Seed                            1
SphereMode                      0       
WhichSpectrum                   1
WhichTransfer                   0         
FileWithInputSpectrum           data_picola_comparison/input_pofk.dat
FileWithInputTransfer           -
Omega                           0.3175
OmegaBaryon                     0.049
HubbleParam                     0.671
Sigma8                          0.7115
PrimordialIndex                 0.966
UnitLength_in_cm                3.085678e24
UnitMass_in_g                   1.989e43   
UnitVelocity_in_cm_per_s        1e5        
InputSpectrum_UnitLength_in_cm  3.085678e24

pofk_compute_every_step         1 
pofk_compute_rsd_pofk           0
pofk_nbins                      64 
pofk_bintype                    0 
pofk_subtract_shotnoise         1 
pofk_kmin                       0.
pofk_kmax                       0.785375
"
  
  echo "$paramfile"                         
}

function make_parameter_file_fofr(){
  FileBase="$1"
  OutputDir="$2"

  paramfile="
modified_gravity_active         0
fofr0                           1.0
nfofr                           1.0
include_screening               1 
use_lcdm_growth_factors         1
input_pofk_is_for_lcdm          1
input_sigma8_is_for_lcdm        1
OutputDir                       $OutputDir
FileBase                        $FileBase
OutputRedshiftFile              data_picola_comparison/output_redshifts_n20.dat
NumFilesWrittenInParallel       1  
UseCOLA                         1           
Buffer                          1.5 
Nmesh                           128      
Nsample                         128      
Box                             512.0   
Init_Redshift                   29.0    
Seed                            1
SphereMode                      0       
WhichSpectrum                   1
WhichTransfer                   0         
FileWithInputSpectrum           data_picola_comparison/input_pofk.dat
FileWithInputTransfer           -
Omega                           0.3175
OmegaBaryon                     0.049
HubbleParam                     0.671
Sigma8                          0.7115
PrimordialIndex                 0.966
UnitLength_in_cm                3.085678e24
UnitMass_in_g                   1.989e43   
UnitVelocity_in_cm_per_s        1e5        
InputSpectrum_UnitLength_in_cm  3.085678e24

pofk_compute_every_step         1 
pofk_compute_rsd_pofk           0
pofk_nbins                      64 
pofk_bintype                    0 
pofk_subtract_shotnoise         1 
pofk_kmin                       0.
pofk_kmax                       0.785375
"

  echo "$paramfile"                         
}

#==================================
# Run code
#==================================
if [[ "$runsim" == "true" ]]; then
  # Recompile code?
  if [[ "$recompile" == "true" ]]; then
    cd ../
    make clean; make MODEL=DGP OPT="-DCOMPUTE_POFK"
    cp $mymgpicolaexec_dgp test_runs
    cd test_runs
    
    cd ../
    make clean; make MODEL=FOFR OPT="-DCOMPUTE_POFK"
    cp $mymgpicolaexec_fofr test_runs
    cd test_runs
  fi

  # Compile code if executable is not found
  if [[ ! -e $mymgpicolaexec_dgp ]]; then
    cd ../
    make clean; make MODEL=DGP OPT="-DCOMPUTE_POFK"
    cp $mymgpicolaexec_dgp test_runs
    cd test_runs
  fi
  
  if [[ ! -e $mymgpicolaexec_fofr ]]; then
    cd ../
    make clean; make MODEL=FOFR OPT="-DCOMPUTE_POFK"
    cp $mymgpicolaexec_fofr test_runs
    cd test_runs
  fi

  # Make output directory
  if [[ ! -e $picolaoutputdir ]]; then
    mkdir $picolaoutputdir
  else
    rm $picolaoutputdir/*z0p000*
  fi

  # Run LCDM simulation (DGP version)
  paramfile=$( make_parameter_file_dgp lcdm_dgp $picolaoutputdir )
  echo "$paramfile" > $picolaoutputdir/lcdm.inp
  mpirun -np $ncpu $mymgpicolaexec_dgp $picolaoutputdir/lcdm.inp
  
  # Run LCDM simulations (f(R) version)
  paramfile=$( make_parameter_file_fofr lcdm_fofr $picolaoutputdir )
  echo "$paramfile" > $picolaoutputdir/lcdm.inp
  mpirun -np $ncpu $mymgpicolaexec_fofr $picolaoutputdir/lcdm.inp
fi

#==================================
# Make plot
#==================================
if [[ "$makeplot" == "true" ]]; then
  
  gnuplotscript="
    set terminal post eps color enhanced font \"Helvetica,20\" linewidth 5 dashlength 4
    set output \"$plotname.eps\"
    
    set log x
    set xrange[2*pi/512.0:pi/512.0 * 128]
    set yrange[0.95:1.05]
    set ylabel \"P_{MG-PICOLA}(k) / P_{PICOLA}(k)\"
    set xlabel \"k   (h/Mpc)\"
    
    plot \\
         '<paste data_picola_comparison/picola_pofk_z0.000_CDM_n20.txt $picolaoutputdir/pofk_lcdm_dgp_z0.000_CDM.txt'  u 3:(\$4/\$2/512.**3) w l ti \"DGP version\" dashtype 1, \\
         '<paste data_picola_comparison/picola_pofk_z0.000_CDM_n20.txt $picolaoutputdir/pofk_lcdm_fofr_z0.000_CDM.txt' u 3:(\$4/\$2/512.**3) w l ti \"f(R) version\" dashtype 2, \\
         1 lt 0 not

    set output"

  # Run gnuplot
  echo "$gnuplotscript" | gnuplot
  
  # Convert to PDF
  $myepspdf $plotname.eps
fi


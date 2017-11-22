#!/bin/bash

#=============================================
# This is for running a simple
# test of the code for f(R). 
# Runs the code, computes P(k) and
# makes a plot of P / P_LCDM 
#=============================================

#==================================
# Paths and options
#==================================
picolaoutputdir="output_fofr"
mymgpicolaexec="MG_PICOLA_FOFR"
mympirun="mpirun-openmpi-gcc6"
myepspdf="epspdf"
plotname="pofk_fofr_testrun"
recompile="true"
makeplot="true"
runsim="true"
ncolasteps="10"
fofr0="1e-5"
ngrid="64"
npart="64"
ncpu="4"
box="200.0"

#==================================
# Function to generate parameterfile
#==================================
function make_parameter_file(){
  modified_gravity_active="$1"
  include_screening="$2"
  use_lcdm_growth_factors="$3"
  input_sigma8_is_for_lcdm="$4"
  fofr0="$5"
  box="$6"
  ngrid="$7"
  npart="$8"
  FileBase="$9"
  OutputDir="${10}"

  paramfile="
modified_gravity_active         $modified_gravity_active 
fofr0                           $fofr0
nfofr                           1.0
include_screening               $include_screening 
use_lcdm_growth_factors         $use_lcdm_growth_factors 
input_pofk_is_for_lcdm          1
input_sigma8_is_for_lcdm        $input_sigma8_is_for_lcdm 
OutputDir                       $OutputDir
FileBase                        $FileBase
OutputRedshiftFile              $OutputDir/output_redshifts.dat
NumFilesWrittenInParallel       1  
UseCOLA                         1           
Buffer                          2.25    
Nmesh                           $ngrid      
Nsample                         $npart      
Box                             $box   
Init_Redshift                   19.0    
Seed                            5001    
SphereMode                      0       
WhichSpectrum                   1
WhichTransfer                   0         
FileWithInputSpectrum           ../files/input_power_spectrum.dat
FileWithInputTransfer           -
Omega                           0.267 
OmegaBaryon                     0.049 
HubbleParam                     0.71  
Sigma8                          0.8   
PrimordialIndex                 0.966 
UnitLength_in_cm                3.085678e24
UnitMass_in_g                   1.989e43   
UnitVelocity_in_cm_per_s        1e5        
InputSpectrum_UnitLength_in_cm  3.085678e24

pofk_compute_every_step         1 
pofk_compute_rsd_pofk           0
pofk_nbins                      0 
pofk_bintype                    0 
pofk_subtract_shotnoise         1 
pofk_kmin                       0.
pofk_kmax                       0.
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
    make clean; make MODEL=FOFR OPT="-DCOMPUTE_POFK"
    cp $mymgpicolaexec test_runs/
    cd test_runs
  fi

  # Compile code if executable is not found
  if [[ ! -e $mymgpicolaexec ]]; then
    cd ../
    make clean; make MODEL=FOFR OPT="-DCOMPUTE_POFK"
    cp $mymgpicolaexec test_runs/
    cd test_runs
  fi

  # Make output directory
  if [[ ! -e $picolaoutputdir ]]; then
    mkdir $picolaoutputdir
  else
    rm $picolaoutputdir/*z0p000*
  fi
  
  # Make step/output file
  echo "0, $ncolasteps" > $picolaoutputdir/output_redshifts.dat

  # Run LCDM simulation
  paramfile=$( make_parameter_file 0 0 1 1 1e-10 $box $ngrid $npart lcdm $picolaoutputdir )
  echo "$paramfile" > $picolaoutputdir/lcdm.inp
  mpirun -np $ncpu $mymgpicolaexec $picolaoutputdir/lcdm.inp

  # Run f(R) simulation (no screening)
  paramfile=$( make_parameter_file 1 0 0 1 $fofr0 $box $ngrid $npart fofr $picolaoutputdir )
  echo "$paramfile" > $picolaoutputdir/fofr.inp
  mpirun -np $ncpu $mymgpicolaexec $picolaoutputdir/fofr.inp
  
  # Run f(R) simulation (with screening)
  paramfile=$( make_parameter_file 1 1 0 1 $fofr0 $box $ngrid $npart fofr_screen $picolaoutputdir )
  echo "$paramfile" > $picolaoutputdir/fofr_screen.inp
  mpirun -np $ncpu $mymgpicolaexec $picolaoutputdir/fofr_screen.inp
  
  # Run f(R) simulation with LCDM growthfac (no screening)
  paramfile=$( make_parameter_file 1 0 1 1 $fofr0 $box $ngrid $npart fofr_Dlcdm $picolaoutputdir )
  echo "$paramfile" > $picolaoutputdir/fofr_Dlcdm.inp
  mpirun -np $ncpu $mymgpicolaexec $picolaoutputdir/fofr_Dlcdm.inp
  
  # Run f(R) simulation with LCDM growthfac (with screening)
  paramfile=$( make_parameter_file 1 1 1 1 $fofr0 $box $ngrid $npart fofr_Dlcdm_screen $picolaoutputdir )
  echo "$paramfile" > $picolaoutputdir/fofr_Dlcdm_screen.inp
  mpirun -np $ncpu $mymgpicolaexec $picolaoutputdir/fofr_Dlcdm_screen.inp
fi

#==================================
# Make plot
#==================================
if [[ "$makeplot" == "true" ]]; then
  
  gnuplotscript="
    set terminal post eps color enhanced font \"Helvetica,20\" linewidth 5 dashlength 4
    set output \"$plotname.eps\"
    
    set log x
    set xrange[2*pi/$box:pi/$box * $ngrid]
    set yrange[0.95:1.5]
    set ylabel \"P(k) / P_{{/Symbol L}CDM}(k)\"
    set xlabel \"k   (h/Mpc)\"
    
    label1 = \"{/Symbol L}CDM growth factor; with screening\"
    label2 = \"{/Symbol L}CDM growth factor;   no screening\"
    label3 = \"f(R) growth factor; with screening\"
    label4 = \"f(R) growth factor;   no screening\"
    label5 = \"N-body result (fitting function)\"
    
    fit_f5(x) = ((1.075 + 0.063*log(x/0.1) - 1)*( 12.0*x/(1+12.0*x)  ) + 1.0) 
    
    plot \\
         '<paste $picolaoutputdir/pofk_fofr_Dlcdm_screen_z0.000_CDM.txt $picolaoutputdir/pofk_lcdm_z0.000_CDM.txt'  u (\$1):(\$2/\$6) w l ti label1 dashtype 1, \\
         '<paste $picolaoutputdir/pofk_fofr_Dlcdm_z0.000_CDM.txt        $picolaoutputdir/pofk_lcdm_z0.000_CDM.txt'  u (\$1):(\$2/\$6) w l ti label2 dashtype 2, \\
         '<paste $picolaoutputdir/pofk_fofr_screen_z0.000_CDM.txt       $picolaoutputdir/pofk_lcdm_z0.000_CDM.txt'  u (\$1):(\$2/\$6) w l ti label3 dashtype 3, \\
         '<paste $picolaoutputdir/pofk_fofr_z0.000_CDM.txt              $picolaoutputdir/pofk_lcdm_z0.000_CDM.txt'  u (\$1):(\$2/\$6) w l ti label4 dashtype 4, \\
         fit_f5(x) ti label5 dashtype 5 lc -1, \\
         1 lt 0 not
    
    set output"

  # Run gnuplot
  echo "$gnuplotscript" | gnuplot
  
  # Convert to PDF
  $myepspdf $plotname.eps
fi


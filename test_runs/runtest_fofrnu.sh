#!/bin/bash

#==================================
# Paths and options
#==================================
picolaoutputdir="output_fofrnu"
mymgpicolaexec="MG_PICOLA_FOFRNU"
mympirun="mpirun-openmpi-gcc6"
myepspdf="epspdf"
plotname="pofk_fofrnu_testrun"
recompile="false"
makeplot="true"
runsim="false"
ncolasteps="20"
inputfiledir="inputfiles_fofrnu"
ncpu="4"
box="512.0"
ngrid="128.0"

#==================================
# Run code
#==================================
if [[ "$runsim" == "true" ]]; then
  s
  # Recompile code?
  if [[ "$recompile" == "true" ]]; then
    cd ../
    make clean; make MODEL=FOFRNU OPT="-DCOMPUTE_POFK"
    cp $mymgpicolaexec test_runs/
    cd test_runs
  fi

  # Compile code if executable is not found
  if [[ ! -e $mymgpicolaexec ]]; then
    cd ../
    make clean; make MODEL=FOFRNU OPT="-DCOMPUTE_POFK"
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

  # Run LCDM simulations
  mpirun -np $ncpu $mymgpicolaexec $inputfiledir/lcdm_nu0.0.inp
  mpirun -np $ncpu $mymgpicolaexec $inputfiledir/lcdm_nu0.2.inp
  mpirun -np $ncpu $mymgpicolaexec $inputfiledir/lcdm_nu0.4.inp
  mpirun -np $ncpu $mymgpicolaexec $inputfiledir/lcdm_nu0.6.inp
  
  # Run f(R) simulations
  mpirun -np $ncpu $mymgpicolaexec $inputfiledir/fofr_nu0.0.inp
  mpirun -np $ncpu $mymgpicolaexec $inputfiledir/fofr_nu0.2.inp
  mpirun -np $ncpu $mymgpicolaexec $inputfiledir/fofr_nu0.4.inp
  mpirun -np $ncpu $mymgpicolaexec $inputfiledir/fofr_nu0.6.inp
fi

#==================================
# Make plot
#==================================
if [[ "$makeplot" == "true" ]]; then
  
  gnuplotscript="
    set terminal post eps color enhanced font \"Helvetica,20\" linewidth 5 dashlength 4
    set output \"$plotname.eps\"

    set key left
    set log x
    set xrange[2*pi/$box:pi/$box * $ngrid * 1.5]
    set yrange[0.2:1.7]
    set ylabel \"P(k) / P_{{/Symbol L}CDM}(k)_{m = 0.0 eV}\"
    set xlabel \"k   (h/Mpc)\"
    
    label1 = \"f(R) m = 0.0 eV\"
    label2 = \"f(R) m = 0.2 eV\"
    label3 = \"f(R) m = 0.4 eV\"
    label4 = \"{/Symbol L}CDM m = 0.0 eV\"
    label5 = \"{/Symbol L}CDM m = 0.2 eV\"
    label6 = \"{/Symbol L}CDM m = 0.4 eV\"
    label7 = \"{/Symbol L}CDM m = 0.6 eV\"
    
    plot \\
         '<paste $picolaoutputdir/pofk_lcdm_nu0.0_z0.000_total.txt $picolaoutputdir/pofk_fofr_nu0.0_z0.000_total.txt' u 1:(\$6/\$2) w l ti label1 dashtype 1, \\
         '<paste $picolaoutputdir/pofk_lcdm_nu0.0_z0.000_total.txt $picolaoutputdir/pofk_fofr_nu0.2_z0.000_total.txt' u 1:(\$6/\$2) w l ti label2 dashtype 2, \\
         '<paste $picolaoutputdir/pofk_lcdm_nu0.0_z0.000_total.txt $picolaoutputdir/pofk_fofr_nu0.4_z0.000_total.txt' u 1:(\$6/\$2) w l ti label3 dashtype 3, \\
         1 lt 0 ti label4, \\
         '<paste $picolaoutputdir/pofk_lcdm_nu0.0_z0.000_total.txt $picolaoutputdir/pofk_lcdm_nu0.2_z0.000_total.txt' u 1:(\$6/\$2) w l ti label5 dashtype 4, \\
         '<paste $picolaoutputdir/pofk_lcdm_nu0.0_z0.000_total.txt $picolaoutputdir/pofk_lcdm_nu0.4_z0.000_total.txt' u 1:(\$6/\$2) w l ti label6 dashtype 5, \\
         '<paste $picolaoutputdir/pofk_lcdm_nu0.0_z0.000_total.txt $picolaoutputdir/pofk_lcdm_nu0.6_z0.000_total.txt' u 1:(\$6/\$2) w l ti label7 dashtype 6
    
    set output"

  # Run gnuplot
  echo "$gnuplotscript" | gnuplot
  
  # Convert to PDF
  $myepspdf $plotname.eps
fi


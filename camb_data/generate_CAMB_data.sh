#!/bin/bash

#===========================================#
#                                           #
# Script to generate CAMB data needed for   #
# the massive neutrino module in MG-PICOLA  #
# alongside the transfer-info file that     #
# MG-PICOLA reads to get the paths to the   #
# CAMB transfer files                       #
#                                           #
#===========================================#

#==========================
# Set the parameters
#==========================
suffix="nu0.2"
output_folder="$HOME/local/MG-PICOLA-PUBLIC/camb_data/example_data_${suffix}"
camb_template="$HOME/local/CAMB/HighLExtrapTemplate_lenspotentialCls.dat"
camb_executable="$HOME/local/CAMB/camb"
camb_parameter_file="$HOME/local/MG-PICOLA-PUBLIC/camb_data/sample_${suffix}.ini"
camb_output_root="lcdm_${suffix}"
camb_input_filename="camb_${camb_output_root}.ini"
picola_output_filename="picola_transfer_info_${suffix}.txt"
number_redshifts="50"
z_min="0.0"
z_max="50.0"

#==========================
# Make output folder
#==========================
if [[ ! -e $output_folder ]]; then
  mkdir $output_folder
fi
cd $output_folder

#=================================
# Copy over camb code and template
#=================================
if [[ -e $camb_executable ]]; then
  cp $camb_executable camb
else
  echo "Error: cannot find the camb executable [$camb_template]"
  exit
fi
if [[ -e $camb_template ]]; then
  cp $camb_template HighLExtrapTemplate_lenspotentialCls.dat
else
  echo "Error: cannot find the camb highL template [$camb_template]"
  exit
fi

#=================================
# Get the other CAMB variables
#=================================
if [[ -e $camb_parameter_file ]]; then
  camb_param="$(cat $camb_parameter_file)"
else
  echo "Error: cannot find the camb parameter file [$camb_parameter_file]"
  exit
fi
tmp_file="output_root              = $camb_output_root\ntransfer_num_redshifts   = $number_redshifts"

#=================================
# Make PICOLA transfer-info file
#=================================
picola_file="$(printf "$output_folder $number_redshifts")"

#=================================
# Generate the z-values
#=================================
for i in $(seq 1 $number_redshifts)
do
  
  #=================================
  # Compute current value of z (log-spacing)
  #=================================
  cur_z=$(echo $z_min $z_max $i $number_redshifts | awk '{printf "%8.3f\n", exp( log($1+1) + (log($2+1) - log($1+1)) * ($4 - $3) / ($4-1.0) ) - 1.0 }' )
  cur_z=$(echo $cur_z | awk '{$1=$1;print}') 

  #=================================
  # Print to CAMB input-file
  #=================================
  tmp_file=$(printf "$tmp_file\ntransfer_redshift($i)    = $cur_z")
  tmp_file=$(printf "$tmp_file\ntransfer_filename($i)    = transfer_z${cur_z}.dat")
  tmp_file=$(printf "$tmp_file\ntransfer_matterpower($i) = matterpower_z${cur_z}.dat")

  #=================================
  # For PICOLA (reverse order)
  #=================================
  cur_z=$(echo $z_min $z_max $i $number_redshifts | awk '{printf "%8.3f\n", exp( log($1+1) + (log($2+1) - log($1+1)) * ($3-1.0) / ($4-1.0) ) - 1.0 }' )
  cur_z=$(echo $cur_z | awk '{$1=$1;print}') 
  picola_file=$(printf "$picola_file\n${camb_output_root}_transfer_z${cur_z}.dat  $cur_z")
done

#=================================
# Output picola transfer info file
#=================================
printf "$picola_file" > $picola_output_filename

#=================================
# Output camb input file
#=================================
printf "$camb_param\n $tmp_file" > $camb_input_filename

#=================================
# Run camb
#=================================
./camb $camb_input_filename > camb.log

#=================================
# Clean up
#=================================
rm camb
rm HighLExtrapTemplate_lenspotentialCls.dat

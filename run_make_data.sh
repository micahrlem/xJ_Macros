#!/bin/bash

# Define input files
input_files=(
  output_run2pp_ana468_2024p012_v001_data_calofit_Run47.root
  output_run2pp_ana468_2024p012_v001_data_calofit_Run48.root
  output_run2pp_ana468_2024p012_v001_data_calofit_Run49.root
  output_run2pp_ana468_2024p012_v001_data_calofit_Run50.root
  output_run2pp_ana468_2024p012_v001_data_calofit_Run51.root
  output_run2pp_ana468_2024p012_v001_data_calofit_Run52.root
  output_run2pp_ana468_2024p012_v001_data_calofit_Run53.root
)

# Loop over each file and each R value
for infile in "${input_files[@]}"; do
  for R in 0.4 0.2; do
    echo "Running make_Data_xJ_Flat with R = $R and file = $infile"
    root -l -b -q "make_Data_xJ_Flat.C(\"binning_original.config\", $R, \"$infile\")"
  done
done

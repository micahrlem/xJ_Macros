#!/bin/bash

# Step 1: Merge runs 47 to 53 individually
hadd -f -k output_run2pp_ana468_2024p012_v001_data_calofit_Run47.root /sphenix/user/mmeskowit/Dissertion_Analysis/DataFiles/outputs/Run47*/OutDir*/output_data.root
hadd -f -k output_run2pp_ana468_2024p012_v001_data_calofit_Run48.root /sphenix/user/mmeskowit/Dissertion_Analysis/DataFiles/outputs/Run48*/OutDir*/output_data.root
hadd -f -k output_run2pp_ana468_2024p012_v001_data_calofit_Run49.root /sphenix/user/mmeskowit/Dissertion_Analysis/DataFiles/outputs/Run49*/OutDir*/output_data.root
hadd -f -k output_run2pp_ana468_2024p012_v001_data_calofit_Run50.root /sphenix/user/mmeskowit/Dissertion_Analysis/DataFiles/outputs/Run50*/OutDir*/output_data.root
hadd -f -k output_run2pp_ana468_2024p012_v001_data_calofit_Run51.root /sphenix/user/mmeskowit/Dissertion_Analysis/DataFiles/outputs/Run51*/OutDir*/output_data.root
hadd -f -k output_run2pp_ana468_2024p012_v001_data_calofit_Run52.root /sphenix/user/mmeskowit/Dissertion_Analysis/DataFiles/outputs/Run52*/OutDir*/output_data.root
hadd -f -k output_run2pp_ana468_2024p012_v001_data_calofit_Run53.root /sphenix/user/mmeskowit/Dissertion_Analysis/DataFiles/outputs/Run53*/OutDir*/output_data.root

# Step 2: Combine all intermediate outputs using a wildcard
hadd -f -k output_run2pp_ana468_2024p012_v001_data_calofit_All.root output_run2pp_ana468_2024p012_v001_data_calofit_Run*.root

#!/bin/bash

CONFIG="binning_original.config"
NITER=6
DUMMY_FILE="dummy.root"  # Used if Jet20 doesn't exist (Herwig case)

# Make a tiny dummy root file once if needed
if [[ ! -f "$DUMMY_FILE" ]]; then
  root -l -b -q -e 'TFile f("dummy.root","RECREATE"); f.Close();' >/dev/null 2>&1
fi

declare -a VARIATIONS=("nominal" "negJER" "posJER" "negJES" "posJES")
declare -a GENERATORS=("pythia8" "Herwig")
declare -a RADII=("R2" "R4")

for VAR in "${VARIATIONS[@]}"; do
  for GEN in "${GENERATORS[@]}"; do
    for RVAL in "${RADII[@]}"; do

      # Skip negJER for R2
      #if [[ "$RVAL" == "R2" && "$VAR" == "negJER" ]]; then
        #echo "Skipping ${VAR} for ${RVAL} (not available)."
       # continue
      #fi

      # Build filenames (match your list)
      if [[ "$VAR" == "nominal" ]]; then
        FILE10="Xj_2D_Response_Bins19_${RVAL}_Play_Smear_Flatten_${GEN}-Jet10-Run21-multiR.root"
        FILE20="Xj_2D_Response_Bins19_${RVAL}_Play_Smear_Flatten_${GEN}-Jet20-Run21-multiR.root"
        FILE30="Xj_2D_Response_Bins19_${RVAL}_Play_Smear_Flatten_${GEN}-Jet30-Run21-multiR.root"
      else
        FILE10="Xj_2D_Response_Bins19_${RVAL}_Play_Smear_Flatten_${VAR}_${GEN}-Jet10-Run21-multiR.root"
        FILE20="Xj_2D_Response_Bins19_${RVAL}_Play_Smear_Flatten_${VAR}_${GEN}-Jet20-Run21-multiR.root"
        FILE30="Xj_2D_Response_Bins19_${RVAL}_Play_Smear_Flatten_${VAR}_${GEN}-Jet30-Run21-multiR.root"
      fi

      # Radius-dependent data file
      DATA_FILE="Xj_2D_Response_Bins19_${RVAL}_Data_Flatten_output_run2pp_ana468_2024p012_v001_data_All.Root"

      # Herwig: use dummy Jet20 and set Pythia=false
      if [[ "$GEN" == "Herwig" ]]; then
        INFILE1="$FILE10"   # Jet10
        INFILE2="$FILE30"   # Jet30
        INFILE3="$DUMMY_FILE"  # Jet20 dummy
        PYTHIA_FLAG="false"
      else
        INFILE1="$FILE10"
        INFILE2="$FILE30"
        INFILE3="$FILE20"
        PYTHIA_FLAG="true"
      fi

      echo "Running ${VAR} ${GEN} ${RVAL} with files:"
      echo "  Jet10 (infile1): $INFILE1"
      echo "  Jet30 (infile2): $INFILE2"
      echo "  Jet20 (infile3): $INFILE3"
      echo "  Data           : $DATA_FILE"
      echo "  Pythia flag    : $PYTHIA_FLAG"
      echo

      root -q -b "Combine_Skim_Unfold_Write_xJ_Play.C(\"${INFILE1}\", \"${INFILE2}\", \"${INFILE3}\", \"${DATA_FILE}\", \"${CONFIG}\", ${NITER}, ${PYTHIA_FLAG})"

    done
  done
done

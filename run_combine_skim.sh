#!/bin/bash

CONFIG="binning_original.config"
NITER=6
DUMMY_FILE="dummy.root"  # Used if Jet20 doesn't exist

declare -a VARIATIONS=("nominal" "negJER" "posJER" "negJES" "posJES")
declare -a GENERATORS=("pythia8" "Herwig")
declare -a RADII=("R2" "R4")

for VAR in "${VARIATIONS[@]}"; do
  for GEN in "${GENERATORS[@]}"; do
    for RVAL in "${RADII[@]}"; do

      # Adjust file tag for nominal or variation
      if [[ "$VAR" == "nominal" ]]; then
        FILE10="Xj_2D_Response_Bins19_${RVAL}_Smear_Flatten_${GEN}-Jet10-Run21-multiR.root"
        FILE20="Xj_2D_Response_Bins19_${RVAL}_Smear_Flatten_${GEN}-Jet20-Run21-multiR.root"
        FILE30="Xj_2D_Response_Bins19_${RVAL}_Smear_Flatten_${GEN}-Jet30-Run21-multiR.root"
      else
        FILE10="Xj_2D_Response_Bins19_${RVAL}_Smear_Flatten__${VAR}_${GEN}-Jet10-Run21-multiR.root"
        FILE20="Xj_2D_Response_Bins19_${RVAL}_Smear_Flatten__${VAR}_${GEN}-Jet20-Run21-multiR.root"
        FILE30="Xj_2D_Response_Bins19_${RVAL}_Smear_Flatten__${VAR}_${GEN}-Jet30-Run21-multiR.root"
      fi

      # Radius-dependent data file
      DATA_FILE="Xj_2D_Response_Bins19_${RVAL}_Data_Flatten_output_run2pp_ana468_2024p012_v001_data_All.Root"

      # Use dummy file for Herwig Jet20
      if [[ "$GEN" == "Herwig" ]]; then
        FILE20="$DUMMY_FILE"
        PYTHIA_FLAG="false"
      else
        PYTHIA_FLAG="true"
      fi

      echo "Running ${VAR} ${GEN} ${RVAL} with files:"
      echo "  Jet10: $FILE10"
      echo "  Jet20: $FILE20"
      echo "  Jet30: $FILE30"
      echo "  Data : $DATA_FILE"
      echo "  Pythia flag: $PYTHIA_FLAG"

      root -q -b "Combine_Skim_Unfold_Write_xJ_Pythia.C(\"${FILE10}\", \"${FILE30}\", \"${FILE20}\", \"${DATA_FILE}\", \"${CONFIG}\", ${NITER}, ${PYTHIA_FLAG})"
    done
  done
done

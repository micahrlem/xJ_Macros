#!/bin/bash

CONFIG="binning_original.config"
NITER=6
DUMMY_FILE="dummy.root"  # Jet20 not used for Herwig
GENERATOR="Herwig"
PYTHIA_FLAG="false"

declare -a VARIATIONS=("nominal" "negJER" "posJER" "negJES" "posJES")
declare -a RADII=("R2" "R4")

for VAR in "${VARIATIONS[@]}"; do
  for RVAL in "${RADII[@]}"; do

    if [[ "$VAR" == "nominal" ]]; then
      FILE10="Xj_2D_Response_Bins19_${RVAL}_Smear_Flatten_Herwig_Reweight_Herwig-Jet10-Run21-multiR.root"
      FILE30="Xj_2D_Response_Bins19_${RVAL}_Smear_Flatten_Herwig_Reweight_Herwig-Jet30-Run21-multiR.root"
    else
      FILE10="Xj_2D_Response_Bins19_${RVAL}_Smear_Flatten_Herwig_Reweight__${VAR}_Herwig-Jet10-Run21-multiR.root"
      FILE30="Xj_2D_Response_Bins19_${RVAL}_Smear_Flatten_Herwig_Reweight__${VAR}_Herwig-Jet30-Run21-multiR.root"
    fi

    DATA_FILE="Xj_2D_Response_Bins19_${RVAL}_Data_Flatten_output_run2pp_ana468_2024p012_v001_data_All.Root"

    echo "Running ${VAR} ${GENERATOR} ${RVAL} with files:"
    echo "  Jet10: $FILE10"
    echo "  Jet20: $DUMMY_FILE (not used)"
    echo "  Jet30: $FILE30"
    echo "  Data : $DATA_FILE"

    root -l -q -b "Combine_Skim_Unfold_Write_xJ_Pythia.C(\"${FILE10}\", \"${FILE30}\", \"${DUMMY_FILE}\", \"${DATA_FILE}\", \"${CONFIG}\", ${NITER}, ${PYTHIA_FLAG})"
  done
done

#!/bin/bash

CONFIGS=(
  binning_negJER.config
  binning_negJES.config
  binning_original.config
  binning_posJER.config
  binning_posJES.config
)

RVALUES=(0.4 0.2)

FILES=(
  pythia8-Jet10-Run21-multiR.root
  pythia8-Jet20-Run21-multiR.root
  pythia8-Jet30-Run21-multiR.root
  Herwig-Jet10-Run21-multiR.root
  Herwig-Jet30-Run21-multiR.root
)

for config in "${CONFIGS[@]}"; do
  for R in "${RVALUES[@]}"; do
    # Skip the invalid combo
    if [[ "$config" == "binning_negJER.config" && "$R" == "0.2" ]]; then
      echo "Skipping $config with R=$R (invalid combo)"
      continue
    fi
    for file in "${FILES[@]}"; do
      echo "Running with config=$config, R=$R, file=$file"
      root -l -b -q "make_xj_Response_smear_flatten_play.C(\"$config\", $R, \"$file\")"
    done
  done
done

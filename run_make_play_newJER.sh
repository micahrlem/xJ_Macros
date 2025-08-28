CONFIGS=(
  binning_negJER.config
  binning_posJER.config
)

RVALUES=(0.2)

FILES=(
  pythia8-Jet10-Run21-multiR.root
  pythia8-Jet20-Run21-multiR.root
  pythia8-Jet30-Run21-multiR.root
  Herwig-Jet10-Run21-multiR.root
  Herwig-Jet30-Run21-multiR.root
)

for config in "${CONFIGS[@]}"; do
  for R in "${RVALUES[@]}"; do
    for file in "${FILES[@]}"; do
      echo "Running with config=$config, R=$R, file=$file"
      root -l -b -q "make_xj_Response_smear_flatten_play.C(\"$config\", $R, \"$file\")"
    done
  done
done

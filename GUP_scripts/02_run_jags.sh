#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# run_jags.sh  —  runs the GUP JAGS R script for a single grid cell
# Usage:  bash run_jags.sh <ecoreg_cellid>
# e.g.:   bash run_jags.sh 2.2.2_668891
#
# Note: R --vanilla --args causes commandArgs() to return:
#   args[1] = /path/to/R
#   args[2] = --vanilla
#   args[3] = <ecoreg_cellid>   ← what the R script indexes into
# ─────────────────────────────────────────────────────────────────────────────

ARG=$1
RSCRIPT="/projectnb/modislc/users/seamorez/HLS_Pheno/scripts/LSP_clim_sensitivity/GUP_scripts/02_run_jags.R"

echo "Starting JAGS run for: ${ARG}"
echo "Time: $(date)"

R --vanilla --args ${ARG} < ${RSCRIPT}

echo "Finished: ${ARG}"
echo "Time: $(date)"
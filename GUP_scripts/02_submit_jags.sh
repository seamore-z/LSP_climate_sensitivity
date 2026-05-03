#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# submit_jags.sh  —  submits one JAGS batch job per .rda input file
#                    for a given ecoregion
# Usage: bash submit_jags.sh
# Edit TT to change ecoregion
# ─────────────────────────────────────────────────────────────────────────────

# ── Configuration — edit these ───────────────────────────────────────────────
TT='2.2.2'
INPUT_DIR="/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/jags_input/${TT}"
RUN_SCRIPT="/projectnb/modislc/users/seamorez/HLS_Pheno/scripts/LSP_clim_sensitivity/GUP_scripts/02_run_jags.sh"
LOG_DIR="/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/logs/${TT}"
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "${LOG_DIR}"

# Confirm input files exist before submitting anything
n_files=$(ls "${INPUT_DIR}/${TT}"_*.rda 2>/dev/null | wc -l)
if [ "$n_files" -eq 0 ]; then
    echo "ERROR: no .rda files found in ${INPUT_DIR} matching ${TT}_*.rda"
    exit 1
fi
echo "Found ${n_files} input files for ecoregion ${TT}"

n_submitted=0

for f in "${INPUT_DIR}/${TT}"_*.rda; do
    [ -f "$f" ] || continue
    
    base=$(basename "$f" .rda)
    vp="${base##*_}"
    ARG="${TT}_${vp}"
    
    echo "  Submitting: ${ARG}"
    
    qsub_opts=(
      -P modislc
      -N "jags_${vp}"
      -l h_rt=12:00:00
      -pe omp 8
      -l mem_per_core=8G
      -l 'h=!(scc-dd1|scc-dd2|scc-dd3|scc-dd4)'
      -o "${LOG_DIR}/jags_${vp}.log"
      -e "${LOG_DIR}/jags_${vp}.err"
    )
    
    qsub "${qsub_opts[@]}" "${RUN_SCRIPT}" "${ARG}"
    
    n_submitted=$((n_submitted + 1))
done

echo "Submitted ${n_submitted} / ${n_files} jobs for ecoregion ${TT}"
echo "Logs → ${LOG_DIR}"
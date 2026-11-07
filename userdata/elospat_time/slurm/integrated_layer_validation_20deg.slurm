#!/bin/bash
#SBATCH --job-name=rms_layer_val20
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/integrated_layer_validation_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/integrated_layer_validation_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
MATLAB=/opt/mathworks/R2021b/bin/matlab
OUT="$ROOT/results/integrated_layer_validation_20deg_job${SLURM_JOB_ID}"

mkdir -p "$OUT"
cd "$ROOT"
if [[ -n "${RESUME_ZERO_DIR:-}" ]]; then
    "$MATLAB" -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$ROOT/planning'); run_integrated_layer_validation_20deg('$OUT', false, '$RESUME_ZERO_DIR');"
else
    "$MATLAB" -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$ROOT/planning'); run_integrated_layer_validation_20deg('$OUT');"
fi

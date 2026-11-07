#!/bin/bash
#SBATCH --job-name=elospat_nominal
#SBATCH --output=/mnt/beegfs/ydang/project1/logs/elospat_nominal_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/logs/elospat_nominal_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail

: "${ELOSPAT_ROOT:=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time}"
: "${MATRAD_ROOT:=/mnt/beegfs/ydang/project1/data/matRad}"
: "${RUN_ROOT:=$ELOSPAT_ROOT/results/nominal_then_robust_${SLURM_JOB_ID}}"
: "${NOMINAL_OUT:=$RUN_ROOT/nominal}"
: "${NOMINAL_ARTIFACT_POINTER:=$RUN_ROOT/nominal_artifact_path.txt}"
: "${ARC_SPACING:=20}"

export ELOSPAT_ROOT MATRAD_ROOT RUN_ROOT NOMINAL_OUT NOMINAL_ARTIFACT_POINTER ARC_SPACING

LOG=/mnt/beegfs/ydang/project1/logs
mkdir -p "$LOG" "$RUN_ROOT" "$NOMINAL_OUT"

module purge
module load matlab

cd "$ELOSPAT_ROOT"

echo "[$(date)] Starting nominal-only RMS-ELO-SPAT"
echo "RUN_ROOT=$RUN_ROOT"
echo "NOMINAL_OUT=$NOMINAL_OUT"
echo "NOMINAL_ARTIFACT_POINTER=$NOMINAL_ARTIFACT_POINTER"
echo "ARC_SPACING=$ARC_SPACING"

matlab -softwareopengl -batch "addpath('$ELOSPAT_ROOT/slurm'); run_slurm_nominal_only"

echo "[$(date)] Nominal-only RMS-ELO-SPAT completed"

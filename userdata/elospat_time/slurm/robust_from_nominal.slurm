#!/bin/bash
#SBATCH --job-name=elospat_robust
#SBATCH --output=/mnt/beegfs/ydang/project1/logs/elospat_robust_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/logs/elospat_robust_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail

: "${ELOSPAT_ROOT:=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time}"
: "${MATRAD_ROOT:=/mnt/beegfs/ydang/project1/data/matRad}"
: "${RUN_ROOT:=$ELOSPAT_ROOT/results/nominal_then_robust_${SLURM_JOB_ID}}"
: "${ROBUST_OUT:=$RUN_ROOT/robust}"
: "${NOMINAL_ARTIFACT_POINTER:=$RUN_ROOT/nominal_artifact_path.txt}"
: "${SAVE_FULL_DIJ:=false}"

export ELOSPAT_ROOT MATRAD_ROOT RUN_ROOT ROBUST_OUT NOMINAL_ARTIFACT_POINTER SAVE_FULL_DIJ

LOG=/mnt/beegfs/ydang/project1/logs
mkdir -p "$LOG" "$RUN_ROOT" "$ROBUST_OUT"

module purge
module load matlab

cd "$ELOSPAT_ROOT"

echo "[$(date)] Starting robust evaluation from nominal artifact"
echo "RUN_ROOT=$RUN_ROOT"
echo "ROBUST_OUT=$ROBUST_OUT"
echo "NOMINAL_ARTIFACT_POINTER=$NOMINAL_ARTIFACT_POINTER"
echo "SAVE_FULL_DIJ=$SAVE_FULL_DIJ"

matlab -softwareopengl -batch "addpath('$ELOSPAT_ROOT/slurm'); run_slurm_robust_from_nominal"

echo "[$(date)] Robust evaluation completed"

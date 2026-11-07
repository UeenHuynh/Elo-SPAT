#!/bin/bash
#SBATCH --job-name=elospat_rob39_42
#SBATCH --output=/mnt/beegfs/ydang/project1/logs/elospat_rob39_42_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/logs/elospat_rob39_42_%j.err
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail

: "${ELOSPAT_ROOT:=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time}"
: "${MATRAD_ROOT:=/mnt/beegfs/ydang/project1/data/matRad}"
: "${RUN_ROOT:=$ELOSPAT_ROOT/results/full_nominal_then_robust_${SLURM_JOB_ID}}"
: "${NOMINAL_ROOT:=$RUN_ROOT/nominal}"
: "${ROBUST_ROOT:=$RUN_ROOT/robust}"
: "${TASK39_NOMINAL_ARTIFACT_POINTER:=$RUN_ROOT/task39_nominal_artifact_path.txt}"
: "${SAVE_FULL_DIJ:=false}"

export ELOSPAT_ROOT MATRAD_ROOT RUN_ROOT NOMINAL_ROOT ROBUST_ROOT TASK39_NOMINAL_ARTIFACT_POINTER SAVE_FULL_DIJ

LOG=/mnt/beegfs/ydang/project1/logs
mkdir -p "$LOG" "$RUN_ROOT" "$ROBUST_ROOT"

module purge
module load matlab

cd "$ELOSPAT_ROOT"

echo "[$(date)] Starting full robust Tasks 39-42"
echo "RUN_ROOT=$RUN_ROOT"
echo "NOMINAL_ROOT=$NOMINAL_ROOT"
echo "ROBUST_ROOT=$ROBUST_ROOT"
echo "TASK39_NOMINAL_ARTIFACT_POINTER=$TASK39_NOMINAL_ARTIFACT_POINTER"
echo "SAVE_FULL_DIJ=$SAVE_FULL_DIJ"

matlab -softwareopengl -batch "addpath('$ELOSPAT_ROOT/slurm'); run_slurm_full_robust_39_42"

echo "[$(date)] Full robust Tasks 39-42 completed"

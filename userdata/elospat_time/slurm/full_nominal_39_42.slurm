#!/bin/bash
#SBATCH --job-name=elospat_nom39_42
#SBATCH --output=/mnt/beegfs/ydang/project1/logs/elospat_nom39_42_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/logs/elospat_nom39_42_%j.err
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail

: "${ELOSPAT_ROOT:=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time}"
: "${MATRAD_ROOT:=/mnt/beegfs/ydang/project1/data/matRad}"
: "${RUN_ROOT:=$ELOSPAT_ROOT/results/full_nominal_then_robust_${SLURM_JOB_ID}}"
: "${NOMINAL_ROOT:=$RUN_ROOT/nominal}"
: "${TASK39_NOMINAL_ARTIFACT_POINTER:=$RUN_ROOT/task39_nominal_artifact_path.txt}"
: "${ARC_SPACING:=20}"

export ELOSPAT_ROOT MATRAD_ROOT RUN_ROOT NOMINAL_ROOT TASK39_NOMINAL_ARTIFACT_POINTER ARC_SPACING

LOG=/mnt/beegfs/ydang/project1/logs
mkdir -p "$LOG" "$RUN_ROOT" "$NOMINAL_ROOT"

module purge
module load matlab

cd "$ELOSPAT_ROOT"

echo "[$(date)] Starting full nominal Tasks 39-42"
echo "RUN_ROOT=$RUN_ROOT"
echo "NOMINAL_ROOT=$NOMINAL_ROOT"
echo "TASK39_NOMINAL_ARTIFACT_POINTER=$TASK39_NOMINAL_ARTIFACT_POINTER"
echo "ARC_SPACING=$ARC_SPACING"

matlab -softwareopengl -batch "addpath('$ELOSPAT_ROOT/slurm'); run_slurm_full_nominal_39_42"

echo "[$(date)] Full nominal Tasks 39-42 completed"

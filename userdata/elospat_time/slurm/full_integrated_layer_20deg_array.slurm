#!/bin/bash
#SBATCH --job-name=rms_full_layer20
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --array=0-2%1
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/full_integrated_layer_20deg_%A_%a.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/full_integrated_layer_20deg_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --requeue

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
MATLAB=/opt/mathworks/R2021b/bin/matlab
case "$SLURM_ARRAY_TASK_ID" in
  0) ALPHA=0.1; TOKEN=0p1 ;;
  1) ALPHA=0.3; TOKEN=0p3 ;;
  2) ALPHA=1;   TOKEN=1 ;;
  *) echo "Unexpected array task: $SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac
OUT="$ROOT/results/full_integrated_layer_20deg_job${SLURM_ARRAY_JOB_ID}/alpha_${TOKEN}"
mkdir -p "$OUT"
cd "$ROOT"
"$MATLAB" -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$ROOT/planning'); run_full_integrated_layer_20deg('$OUT',$ALPHA);"

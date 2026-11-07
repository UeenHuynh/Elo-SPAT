#!/bin/bash
#SBATCH --job-name=rms_spot_cont
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --array=0-1%1
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/stage2_spot_continuation_20deg_%A_%a.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/stage2_spot_continuation_20deg_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --requeue

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
MATLAB=/opt/mathworks/R2021b/bin/matlab
case "$SLURM_ARRAY_TASK_ID" in
  0) ALPHA=0.01; TOKEN=0p01 ;;
  1) ALPHA=0.03; TOKEN=0p03 ;;
  *) echo "Unexpected array task: $SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac
OUT="$ROOT/results/stage2_spot_continuation_20deg_job${SLURM_ARRAY_JOB_ID}/alphaS_${TOKEN}"
mkdir -p "$OUT"
cd "$ROOT"
"$MATLAB" -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$ROOT/planning'); run_stage2_spot_continuation_20deg('$OUT',$ALPHA);"

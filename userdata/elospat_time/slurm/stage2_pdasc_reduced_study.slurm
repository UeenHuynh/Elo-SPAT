#!/bin/bash
#SBATCH --job-name=rms_pdasc
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/stage2_pdasc_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/stage2_pdasc_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
MATLAB=/opt/mathworks/R2021b/bin/matlab
OUT="$ROOT/results/stage2_pdasc_reduced_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT"
"$MATLAB" -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$ROOT/planning'); run_stage2_pdasc_reduced_study('$OUT');"

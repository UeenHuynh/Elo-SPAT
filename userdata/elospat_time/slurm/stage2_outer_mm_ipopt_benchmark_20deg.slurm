#!/bin/bash
#SBATCH --job-name=rms_outer_mm
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/stage2_outer_mm_ipopt_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/stage2_outer_mm_ipopt_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
MATLAB=/opt/mathworks/R2021b/bin/matlab
OUT="$ROOT/results/stage2_outer_mm_ipopt_20deg_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT"
"$MATLAB" -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$ROOT/planning'); run_stage2_outer_mm_ipopt_benchmark_20deg('$OUT');"

#!/bin/bash
#SBATCH --job-name=rms_convsens20
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/convergence_sensitivity_20deg_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/convergence_sensitivity_20deg_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATLAB=/opt/mathworks/R2021b/bin/matlab
OUT="$ROOT/results/convergence_sensitivity_20deg_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT"
"$MATLAB" -softwareopengl -batch "addpath('$ROOT/planning'); run_convergence_sensitivity_analysis_20deg('$OUT');"

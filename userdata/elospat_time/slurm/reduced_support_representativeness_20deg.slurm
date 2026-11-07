#!/bin/bash
#SBATCH --job-name=rms_support_repr20
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/support_representativeness_20deg_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/support_representativeness_20deg_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
MATLAB=/opt/mathworks/R2021b/bin/matlab
OUT="$ROOT/results/support_representativeness_20deg_job${SLURM_JOB_ID}"
PRIOR="$ROOT/results/integrated_layer_validation_20deg_job21835"

mkdir -p "$OUT"
cd "$ROOT"
"$MATLAB" -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$ROOT/planning'); run_reduced_support_representativeness_20deg('$OUT', '$PRIOR');"

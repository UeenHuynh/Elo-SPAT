#!/bin/bash
#SBATCH --job-name=pareto_boundary
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/pareto_boundary_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/pareto_boundary_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=190G
set -euo pipefail
PROJECT=/mnt/beegfs/ydang/project1/data/Elo-SPAT
ROOT="$PROJECT/userdata/elospat_time"
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
EXP="$ROOT/experiments/pareto_safe_boundary_continuation"
OUT="$ROOT/results/pareto_safe_boundary_continuation_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT/results/survivor_compression_boundary_job22813"
sha256sum -c FROZEN_RESULT.sha256
cd "$ROOT/results/adaptive_survivor_boundary_job22926"
sha256sum -c FROZEN_RESULT.sha256
cd "$PROJECT"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "addpath('$MATRAD');matRad_rc;addpath('$ROOT/experiments/repeated_dormant_replacement');addpath('$EXP');run_pareto_safe_boundary_continuation('$OUT');"

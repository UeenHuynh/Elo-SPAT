#!/bin/bash
#SBATCH --job-name=repair_response_map
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/repair_response_map_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/repair_response_map_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=190G
set -euo pipefail
PROJECT=/mnt/beegfs/ydang/project1/data/Elo-SPAT
ROOT="$PROJECT/userdata/elospat_time"
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
EXP="$ROOT/experiments/repair_response_calibration"
OUT="$ROOT/results/repair_response_calibration_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT/results/survivor_compression_boundary_job22813"
sha256sum -c FROZEN_RESULT.sha256
cd "$ROOT/results/bdt_aware_active_gate_job23262"
sha256sum -c FROZEN_RESULT.sha256
cd "$PROJECT"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "addpath('$MATRAD');matRad_rc;addpath('$ROOT/experiments/repeated_dormant_replacement');addpath('$ROOT/experiments/bdt_aware_active_gate');addpath('$EXP');run_repair_response_calibration('$OUT');"

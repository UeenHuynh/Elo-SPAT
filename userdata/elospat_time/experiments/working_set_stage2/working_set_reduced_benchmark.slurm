#!/bin/bash
#SBATCH --job-name=rms_workset
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/working_set_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/working_set_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
EXP="$ROOT/experiments/working_set_stage2"
PG="$ROOT/results/stage2_pdasc_reduced_job22283/pgmm_selector_result.mat"
OUT="$ROOT/results/working_set_reduced_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$EXP'); run_working_set_reduced_benchmark('$OUT','$PG');"

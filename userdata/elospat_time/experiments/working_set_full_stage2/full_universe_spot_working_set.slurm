#!/bin/bash
#SBATCH --job-name=rms_ws_full
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/full_ws_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/full_ws_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
EXP="$ROOT/experiments/working_set_full_stage2"
OUT="$ROOT/results/full_universe_working_set_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$EXP'); run_full_universe_spot_working_set('$OUT');"

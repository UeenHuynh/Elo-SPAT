#!/bin/bash
#SBATCH --job-name=tasks40_42
#SBATCH --output=/mnt/beegfs/ydang/project1/logs/tasks40_42_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/logs/tasks40_42_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail
module purge
module load matlab

ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
PLAN="$ROOT/planning"
EVAL="$ROOT/evaluation"
LOG=/mnt/beegfs/ydang/project1/logs

mkdir -p "$LOG"

cd "$ROOT"

matlab -softwareopengl -batch "cd('$ROOT'); addpath('$MATRAD'); matRad_rc; addpath('$PLAN'); run_task40_parameter_sweep(); run_task41_arc_spacing_validation(); addpath('$EVAL'); run_task42_documentation_check();"

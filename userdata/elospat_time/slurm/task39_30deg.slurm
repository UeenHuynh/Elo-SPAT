#!/bin/bash
#SBATCH --job-name=task39_30deg
#SBATCH --output=/mnt/beegfs/ydang/project1/logs/task39_30deg_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/logs/task39_30deg_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail
module purge
module load matlab

ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
PLAN="$ROOT/planning"
LOG=/mnt/beegfs/ydang/project1/logs

mkdir -p "$LOG"

cd "$ROOT"

matlab -softwareopengl -batch "cd('$ROOT'); addpath('$MATRAD'); matRad_rc; addpath('$PLAN'); run_task39_robust_integration(struct('arcSpacing',30));"

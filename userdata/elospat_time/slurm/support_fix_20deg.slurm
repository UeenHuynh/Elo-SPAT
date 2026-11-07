#!/bin/bash
#SBATCH --job-name=elospat_fix20
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/support_fix_20deg_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/support_fix_20deg_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
PLAN="$ROOT/planning"
OUT="$ROOT/results/support_fix_20deg"
MATLAB=/opt/mathworks/R2021b/bin/matlab

mkdir -p "$OUT"
cd "$ROOT"

"$MATLAB" -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$PLAN'); run_task41_arc_spacing_validation(20, [], '$OUT');"

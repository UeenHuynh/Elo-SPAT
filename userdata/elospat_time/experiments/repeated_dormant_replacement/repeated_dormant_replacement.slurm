#!/bin/bash
#SBATCH --job-name=repeated_dormant
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/repeated_dormant_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/repeated_dormant_%j.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=190G
set -euo pipefail
PROJECT=/mnt/beegfs/ydang/project1/data/Elo-SPAT
ROOT="$PROJECT/userdata/elospat_time"
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
EXP="$ROOT/experiments/repeated_dormant_replacement"
OUT="$ROOT/results/repeated_dormant_replacement_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT/results/full_integrated_layer_20deg_job22180/alpha_0p3"
sha256sum -c FROZEN_LAYER_OPERATING_POINT.sha256
cd "$ROOT/results/gate_constrained_fixed_support_job22558"
sha256sum -c FROZEN_RESULT.sha256
cd "$PROJECT"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "addpath('$MATRAD');matRad_rc;addpath('$EXP');files=dir(fullfile('$EXP','*.m'));bad=false;for i=1:numel(files),m=checkcode(fullfile(files(i).folder,files(i).name),'-id');bad=bad||any(strcmp({m.id},'SYNER'));end;assert(~bad);run_repeated_dormant_replacement('$OUT');"

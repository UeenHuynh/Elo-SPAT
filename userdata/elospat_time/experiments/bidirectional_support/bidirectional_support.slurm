#!/bin/bash
#SBATCH --job-name=rms_bidir_ws
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/bidirectional_support_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/bidirectional_support_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
EXP="$ROOT/experiments/bidirectional_support"
OUT="$ROOT/results/bidirectional_support_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT/results/full_integrated_layer_20deg_job22180/alpha_0p3"
sha256sum -c FROZEN_LAYER_OPERATING_POINT.sha256
cd "$ROOT"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$EXP'); files=dir(fullfile('$EXP','*.m')); bad=false; for i=1:numel(files),m=checkcode(fullfile(files(i).folder,files(i).name),'-id');fprintf('CHECKCODE %s: %d\\n',files(i).name,numel(m));bad=bad||any(strcmp({m.id},'SYNER'));end;assert(~bad);run_bidirectional_support_experiment('$OUT');"

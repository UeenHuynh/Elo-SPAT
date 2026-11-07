#!/bin/bash
#SBATCH --job-name=gateaware_e1
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/gate_aware_event1_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/gate_aware_event1_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=190G
set -euo pipefail
PROJECT=/mnt/beegfs/ydang/project1/data/Elo-SPAT
ROOT="$PROJECT/userdata/elospat_time"
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
EXP="$ROOT/experiments/gate_aware_replacement"
OUT="$ROOT/results/gate_aware_event1_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT/results/full_integrated_layer_20deg_job22180/alpha_0p3"
sha256sum -c FROZEN_LAYER_OPERATING_POINT.sha256
cd "$PROJECT"
sha256sum -c "$ROOT/results/dormant_replacement_job22325/ARTIFACTS.sha256"
cd "$ROOT"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$EXP'); m=checkcode(fullfile('$EXP','run_gate_aware_event1.m'),'-id'); for i=1:numel(m),fprintf('%s L%d: %s\n',m(i).id,m(i).line,m(i).message);end; assert(~any(strcmp({m.id},'SYNER'))); run_gate_aware_event1('$OUT');"

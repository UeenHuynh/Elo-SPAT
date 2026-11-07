#!/bin/bash
#SBATCH --job-name=rms_secant_stress
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/collective_secant_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/collective_secant_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
EXP="$ROOT/experiments/secant_witness_circuit"
OUT="$ROOT/results/collective_secant_stress_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT/results/full_integrated_layer_20deg_job22180/alpha_0p3"
sha256sum -c FROZEN_LAYER_OPERATING_POINT.sha256
cd "$ROOT"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "addpath('$MATRAD'); matRad_rc; addpath('$EXP'); files={'run_collective_secant_stress_test.m','solve_collective_exchange_qp.m','test_collective_exchange_core.m'}; bad=false; for i=1:numel(files), m=checkcode(fullfile('$EXP',files{i}),'-id'); fprintf('CHECKCODE %s: %d messages\\n',files{i},numel(m)); bad=bad||any(strcmp({m.id},'SYNER')); end; assert(~bad,'MATLAB syntax error detected.'); test_secant_witness_core; test_collective_exchange_core; run_collective_secant_stress_test('$OUT');"

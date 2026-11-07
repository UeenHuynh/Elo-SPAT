#!/bin/bash
#SBATCH --job-name=gateaware_pre
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/gate_aware_preflight_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/gate_aware_preflight_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
set -euo pipefail
EXP=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/experiments/gate_aware_replacement
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "m=checkcode(fullfile('$EXP','run_gate_aware_event1.m'),'-id'); for i=1:numel(m),fprintf('%s L%d C%d: %s\n',m(i).id,m(i).line,m(i).column,m(i).message);end; assert(~any(strcmp({m.id},'SYNER'))); disp('GATE_AWARE_PREFLIGHT_PASS');"

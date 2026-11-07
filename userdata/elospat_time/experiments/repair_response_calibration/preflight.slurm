#!/bin/bash
#SBATCH --job-name=repair_map_preflight
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/repair_map_preflight_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/repair_map_preflight_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
set -euo pipefail
PROJECT=/mnt/beegfs/ydang/project1/data/Elo-SPAT
EXP="$PROJECT/userdata/elospat_time/experiments/repair_response_calibration"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "m=checkcode(fullfile('$EXP','run_repair_response_calibration.m'),'-id');disp(struct2table(m));assert(~any(strcmp({m.id},'SYNER')));"

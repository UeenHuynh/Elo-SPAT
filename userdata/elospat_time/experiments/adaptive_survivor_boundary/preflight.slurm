#!/bin/bash
#SBATCH --job-name=adaptive_boundary_preflight
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/adaptive_boundary_preflight_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/adaptive_boundary_preflight_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
set -euo pipefail
PROJECT=/mnt/beegfs/ydang/project1/data/Elo-SPAT
EXP="$PROJECT/userdata/elospat_time/experiments/adaptive_survivor_boundary"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "m=checkcode(fullfile('$EXP','run_adaptive_survivor_boundary.m'),'-id');disp(struct2table(m));assert(~any(strcmp({m.id},'SYNER')));"

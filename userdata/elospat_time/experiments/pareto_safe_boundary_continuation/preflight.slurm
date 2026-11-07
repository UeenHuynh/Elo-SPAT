#!/bin/bash
#SBATCH --job-name=pareto_boundary_preflight
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/pareto_boundary_preflight_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/pareto_boundary_preflight_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
set -euo pipefail
PROJECT=/mnt/beegfs/ydang/project1/data/Elo-SPAT
EXP="$PROJECT/userdata/elospat_time/experiments/pareto_safe_boundary_continuation"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "m=checkcode(fullfile('$EXP','run_pareto_safe_boundary_continuation.m'),'-id');disp(struct2table(m));assert(~any(strcmp({m.id},'SYNER')));"

#!/bin/bash
#SBATCH --job-name=dormant_report
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/dormant_report_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/dormant_report_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
set -euo pipefail
EXP=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/experiments/dormant_replacement
OUT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/dormant_replacement_job22325
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "addpath('$EXP'); finalize_stopped_dormant_experiment('$OUT');"

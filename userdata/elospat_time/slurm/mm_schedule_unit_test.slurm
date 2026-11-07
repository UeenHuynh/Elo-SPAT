#!/bin/bash
#SBATCH --job-name=rms_mm_unit
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/mm_schedule_unit_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/mm_schedule_unit_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT
cd "$ROOT"
/opt/mathworks/R2021b/bin/matlab -nodisplay -nosplash -batch "run('userdata/elospat_time/tests/test_integrated_optimizer_resume.m');"

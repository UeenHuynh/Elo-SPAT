#!/bin/bash
#SBATCH --job-name=rms_outer_unit
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/outer_mm_ipopt_unit_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/outer_mm_ipopt_unit_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT
MATRAD=/mnt/beegfs/ydang/project1/data/matRad
cd "$ROOT"
/opt/mathworks/R2021b/bin/matlab -nodisplay -nosplash -batch "addpath('$MATRAD'); matRad_rc; run('userdata/elospat_time/tests/test_weighted_matrad_ipopt.m');"

#!/bin/bash
#SBATCH --job-name=rms_mmu_audit
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/stage2_mmu_units_audit_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/stage2_mmu_units_audit_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
MATLAB=/opt/mathworks/R2021b/bin/matlab
OUT="$ROOT/results/stage2_mmu_units_audit_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$ROOT"
"$MATLAB" -softwareopengl -batch "addpath('$ROOT/planning'); audit_stage2_mmu_units('$OUT');"

#!/usr/bin/env bash
set -euo pipefail

ROOT="${ELOSPAT_ROOT:-/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time}"
MATRAD="${MATRAD_ROOT:-/mnt/beegfs/ydang/project1/data/matRad}"
ARC_SPACING="${ARC_SPACING:-20}"
SAVE_FULL_DIJ="${SAVE_FULL_DIJ:-false}"
RUN_TAG="${RUN_TAG:-nominal_then_robust_$(date +%Y-%m-%d_%H%M%S)}"
RUN_ROOT="${RUN_ROOT:-$ROOT/results/$RUN_TAG}"
NOMINAL_OUT="${NOMINAL_OUT:-$RUN_ROOT/nominal}"
ROBUST_OUT="${ROBUST_OUT:-$RUN_ROOT/robust}"
NOMINAL_ARTIFACT_POINTER="${NOMINAL_ARTIFACT_POINTER:-$RUN_ROOT/nominal_artifact_path.txt}"

mkdir -p "$RUN_ROOT" "$NOMINAL_OUT" "$ROBUST_OUT" /mnt/beegfs/ydang/project1/logs

nominal_submit=$(
  sbatch --parsable \
    --export=ALL,ELOSPAT_ROOT="$ROOT",MATRAD_ROOT="$MATRAD",RUN_ROOT="$RUN_ROOT",NOMINAL_OUT="$NOMINAL_OUT",NOMINAL_ARTIFACT_POINTER="$NOMINAL_ARTIFACT_POINTER",ARC_SPACING="$ARC_SPACING" \
    "$ROOT/slurm/nominal_only.slurm"
)
nominal_job_id="${nominal_submit%%;*}"

robust_submit=$(
  sbatch --parsable \
    --dependency="afterok:$nominal_job_id" \
    --export=ALL,ELOSPAT_ROOT="$ROOT",MATRAD_ROOT="$MATRAD",RUN_ROOT="$RUN_ROOT",ROBUST_OUT="$ROBUST_OUT",NOMINAL_ARTIFACT_POINTER="$NOMINAL_ARTIFACT_POINTER",SAVE_FULL_DIJ="$SAVE_FULL_DIJ" \
    "$ROOT/slurm/robust_from_nominal.slurm"
)
robust_job_id="${robust_submit%%;*}"

cat <<EOF
Submitted nominal then robust workflow.
RUN_ROOT=$RUN_ROOT
NOMINAL_JOB_ID=$nominal_job_id
ROBUST_JOB_ID=$robust_job_id
DEPENDENCY=afterok:$nominal_job_id
NOMINAL_OUT=$NOMINAL_OUT
ROBUST_OUT=$ROBUST_OUT
NOMINAL_ARTIFACT_POINTER=$NOMINAL_ARTIFACT_POINTER

Check queue:
  squeue -j $nominal_job_id,$robust_job_id

Check logs:
  tail -f /mnt/beegfs/ydang/project1/logs/elospat_nominal_${nominal_job_id}.out
  tail -f /mnt/beegfs/ydang/project1/logs/elospat_robust_${robust_job_id}.out
EOF

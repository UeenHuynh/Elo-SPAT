#!/usr/bin/env bash
set -euo pipefail

ROOT="${ELOSPAT_ROOT:-/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time}"
MATRAD="${MATRAD_ROOT:-/mnt/beegfs/ydang/project1/data/matRad}"
ARC_SPACING="${ARC_SPACING:-20}"
SAVE_FULL_DIJ="${SAVE_FULL_DIJ:-false}"
NOMINAL_TIME="${NOMINAL_TIME:-168:00:00}"
ROBUST_TIME="${ROBUST_TIME:-168:00:00}"
RUN_TAG="${RUN_TAG:-full_nominal_then_robust_$(date +%Y-%m-%d_%H%M%S)}"
RUN_ROOT="${RUN_ROOT:-$ROOT/results/$RUN_TAG}"
NOMINAL_ROOT="${NOMINAL_ROOT:-$RUN_ROOT/nominal}"
ROBUST_ROOT="${ROBUST_ROOT:-$RUN_ROOT/robust}"
TASK39_NOMINAL_ARTIFACT_POINTER="${TASK39_NOMINAL_ARTIFACT_POINTER:-$RUN_ROOT/task39_nominal_artifact_path.txt}"

mkdir -p "$RUN_ROOT" "$NOMINAL_ROOT" "$ROBUST_ROOT" /mnt/beegfs/ydang/project1/logs

nominal_submit=$(
  sbatch --parsable \
    --time="$NOMINAL_TIME" \
    --export=ALL,ELOSPAT_ROOT="$ROOT",MATRAD_ROOT="$MATRAD",RUN_ROOT="$RUN_ROOT",NOMINAL_ROOT="$NOMINAL_ROOT",TASK39_NOMINAL_ARTIFACT_POINTER="$TASK39_NOMINAL_ARTIFACT_POINTER",ARC_SPACING="$ARC_SPACING" \
    "$ROOT/slurm/full_nominal_39_42.slurm"
)
nominal_job_id="${nominal_submit%%;*}"

robust_submit=$(
  sbatch --parsable \
    --dependency="afterok:$nominal_job_id" \
    --time="$ROBUST_TIME" \
    --export=ALL,ELOSPAT_ROOT="$ROOT",MATRAD_ROOT="$MATRAD",RUN_ROOT="$RUN_ROOT",NOMINAL_ROOT="$NOMINAL_ROOT",ROBUST_ROOT="$ROBUST_ROOT",TASK39_NOMINAL_ARTIFACT_POINTER="$TASK39_NOMINAL_ARTIFACT_POINTER",SAVE_FULL_DIJ="$SAVE_FULL_DIJ" \
    "$ROOT/slurm/full_robust_39_42.slurm"
)
robust_job_id="${robust_submit%%;*}"

cat <<EOF
Submitted FULL nominal Tasks 39-42, then FULL robust Tasks 39-42.
RUN_ROOT=$RUN_ROOT
NOMINAL_JOB_ID=$nominal_job_id
ROBUST_JOB_ID=$robust_job_id
DEPENDENCY=afterok:$nominal_job_id
NOMINAL_ROOT=$NOMINAL_ROOT
ROBUST_ROOT=$ROBUST_ROOT
TASK39_NOMINAL_ARTIFACT_POINTER=$TASK39_NOMINAL_ARTIFACT_POINTER
NOMINAL_TIME=$NOMINAL_TIME
ROBUST_TIME=$ROBUST_TIME

Check queue:
  squeue -j $nominal_job_id,$robust_job_id

Check nominal log:
  tail -f /mnt/beegfs/ydang/project1/logs/elospat_nom39_42_${nominal_job_id}.out

Check robust log:
  tail -f /mnt/beegfs/ydang/project1/logs/elospat_rob39_42_${robust_job_id}.out
EOF

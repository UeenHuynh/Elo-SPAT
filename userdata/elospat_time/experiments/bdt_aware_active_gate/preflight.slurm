#!/bin/bash
#SBATCH --job-name=bdt_gate_preflight
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/bdt_gate_preflight_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/bdt_gate_preflight_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
set -euo pipefail
PROJECT=/mnt/beegfs/ydang/project1/data/Elo-SPAT
EXP="$PROJECT/userdata/elospat_time/experiments/bdt_aware_active_gate"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "files=dir(fullfile('$EXP','*.m'));bad=false;for i=1:numel(files),m=checkcode(fullfile(files(i).folder,files(i).name),'-id');disp(files(i).name);disp(struct2table(m));bad=bad||any(strcmp({m.id},'SYNER'));end;assert(~bad);"

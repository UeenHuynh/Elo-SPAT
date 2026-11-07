#!/bin/bash
#SBATCH --job-name=gatefix_pre
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/gate_constrained_preflight_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/gate_constrained_preflight_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
set -euo pipefail
EXP=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/experiments/gate_constrained_fixed_support
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "files=dir(fullfile('$EXP','*.m'));bad=false;for i=1:numel(files),m=checkcode(fullfile(files(i).folder,files(i).name),'-id');fprintf('CHECKCODE %s: %d\n',files(i).name,numel(m));for j=1:numel(m),fprintf('%s L%d: %s\n',m(j).id,m(j).line,m(j).message);end;bad=bad||any(strcmp({m.id},'SYNER'));end;assert(~bad);disp('GATE_CONSTRAINED_PREFLIGHT_PASS');"

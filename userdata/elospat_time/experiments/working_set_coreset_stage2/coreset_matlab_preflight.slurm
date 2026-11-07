#!/bin/bash
#SBATCH --job-name=rms_core_check
#SBATCH --partition=LocalQ
#SBATCH --nodelist=mantis-01
#SBATCH --output=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/coreset_preflight_%j.out
#SBATCH --error=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/coreset_preflight_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail
ROOT=/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time
EXP="$ROOT/experiments/working_set_coreset_stage2"
/opt/mathworks/R2021b/bin/matlab -softwareopengl -batch "files=dir(fullfile('$EXP','*.m')); bad=false; for i=1:numel(files), p=fullfile(files(i).folder,files(i).name); m=checkcode(p,'-id'); fprintf('%s: %d messages\n',files(i).name,numel(m)); for j=1:numel(m), fprintf('  %s L%d C%d: %s\n',m(j).id,m(j).line,m(j).column,m(j).message); if strcmp(m(j).id,'SYNER'),bad=true;end; end; end; assert(~bad,'MATLAB syntax error detected.'); disp('CORESET_MATLAB_PREFLIGHT_PASS');"

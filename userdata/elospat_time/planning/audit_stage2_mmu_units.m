function audit = audit_stage2_mmu_units(outputDir, artifactFile)
%AUDIT_STAGE2_MMU_UNITS Trace matRad optimization weights to physical MU.
%
% This is a read-only audit. It does not optimize, alter an artifact, or
% infer a deliverability threshold from cfg.mmuThreshold. The result is the
% prerequisite for choosing the m_i values in Stage-2 ADMM Variant H.

narginchk(1, 2);
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 2 || isempty(artifactFile)
    artifactFile = fullfile(root, 'results', 'task25_27_pipeline_20deg.mat');
end
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

fprintf('Loading MU metadata from preserved artifact: %s\n', artifactFile);
source = load(artifactFile, 'cfg', 'dij', 'denseResult', 'machine');
dij = source.dij;
n = double(dij.totalNumOfBixels);

assert(isfield(dij, 'numParticlesPerMU') && ~isempty(dij.numParticlesPerMU), ...
    'audit_stage2_mmu_units:missingParticlesPerMU', ...
    'dij.numParticlesPerMU is required to convert w to physical MU.');
nppm = expand_vector(dij.numParticlesPerMU, n, 'numParticlesPerMU');
assert(all(isfinite(nppm) & nppm > 0), ...
    'audit_stage2_mmu_units:invalidParticlesPerMU', ...
    'numParticlesPerMU must be finite and positive.');

if isfield(dij, 'minMU') && ~isempty(dij.minMU)
    machineMinMU = expand_vector(dij.minMU, n, 'minMU');
    minMuSource = 'dij.minMU';
else
    machineMinMU = nan(n, 1);
    minMuSource = 'missing';
end
machineMinW = nppm .* machineMinMU ./ 1e6;

wDense = reshape(double(source.denseResult.w), [], 1);
assert(numel(wDense) == n, 'Dense result size does not match dij.');
physicalMU = wDense .* 1e6 ./ nppm;

planningWeightThreshold = NaN;
if isfield(source.cfg, 'mmuThreshold')
    planningWeightThreshold = double(source.cfg.mmuThreshold);
end
planningMUPerSpot = planningWeightThreshold .* 1e6 ./ nppm;

boundMU = false;
if isfield(source, 'cfg') && isfield(source.cfg, 'pln') && ...
        isfield(source.cfg.pln, 'propOpt') && ...
        isfield(source.cfg.pln.propOpt, 'boundMU')
    boundMU = logical(source.cfg.pln.propOpt.boundMU);
end

audit = struct();
audit.artifactFile = artifactFile;
audit.nSpots = n;
audit.formulaWeightToMU = 'MU_i = 1e6*w_i/numParticlesPerMU_i';
audit.formulaMUToWeight = 'wmin_i = numParticlesPerMU_i*minMU_i/1e6';
audit.machineMinMUSource = minMuSource;
audit.machineMinMU = summarize_vector(machineMinMU);
audit.numParticlesPerMU = summarize_vector(nppm);
audit.machineMinimumWeight = summarize_vector(machineMinW);
audit.planningWeightThreshold = planningWeightThreshold;
audit.planningThresholdPhysicalMU = summarize_vector(planningMUPerSpot);
audit.densePhysicalMU = summarize_vector(physicalMU);
audit.densePositivePhysicalMU = summarize_vector(physicalMU(physicalMU > 0));
audit.configuredBoundMU = boundMU;
audit.machineMMUVerified = all(isfinite(machineMinMU)) && any(machineMinMU > 0);
audit.planningThresholdIsPhysicalMMU = audit.machineMMUVerified && ...
    max(abs(planningMUPerSpot - machineMinMU)) <= ...
    1e-10 * max(1, max(abs(machineMinMU)));

save(fullfile(outputDir, 'stage2_mmu_units_audit.mat'), 'audit', ...
    'machineMinMU', 'machineMinW', 'nppm', 'physicalMU', ...
    'planningMUPerSpot', '-v7.3');
writetable(summary_table(audit), ...
    fullfile(outputDir, 'stage2_mmu_units_summary.csv'));
write_report(fullfile(outputDir, 'STAGE2_MMU_UNITS_AUDIT.md'), audit);

fprintf('Machine MMU verified: %d\n', audit.machineMMUVerified);
fprintf('Planning threshold equals machine MMU: %d\n', ...
    audit.planningThresholdIsPhysicalMMU);
fprintf('dij.minMU [min median max unique] = %.12g %.12g %.12g %d\n', ...
    audit.machineMinMU.min, audit.machineMinMU.median, ...
    audit.machineMinMU.max, audit.machineMinMU.nUnique);
fprintf('numParticlesPerMU [min median max unique] = %.12g %.12g %.12g %d\n', ...
    audit.numParticlesPerMU.min, audit.numParticlesPerMU.median, ...
    audit.numParticlesPerMU.max, audit.numParticlesPerMU.nUnique);
fprintf('cfg.mmuThreshold (weight units) = %.12g\n', planningWeightThreshold);
fprintf('Its physical-MU image [min median max] = %.12g %.12g %.12g\n', ...
    audit.planningThresholdPhysicalMU.min, ...
    audit.planningThresholdPhysicalMU.median, ...
    audit.planningThresholdPhysicalMU.max);
end

function x = expand_vector(x, n, name)
x = reshape(double(x), [], 1);
if isscalar(x)
    x = repmat(x, n, 1);
end
assert(numel(x) == n, 'Field %s has %d entries; expected %d.', ...
    name, numel(x), n);
end

function s = summarize_vector(x)
x = reshape(double(x), [], 1);
xf = x(isfinite(x));
s = struct('n', numel(x), 'nFinite', numel(xf), 'nUnique', 0, ...
    'min', NaN, 'median', NaN, 'max', NaN, 'mean', NaN);
if isempty(xf)
    return;
end
s.nUnique = numel(unique(xf));
s.min = min(xf);
s.median = median(xf);
s.max = max(xf);
s.mean = mean(xf);
end

function t = summary_table(a)
names = {'machine_min_MU'; 'num_particles_per_MU'; ...
    'machine_min_weight'; 'planning_threshold_physical_MU'; ...
    'dense_physical_MU'; 'dense_positive_physical_MU'};
v = {a.machineMinMU; a.numParticlesPerMU; a.machineMinimumWeight; ...
    a.planningThresholdPhysicalMU; a.densePhysicalMU; ...
    a.densePositivePhysicalMU};
n = cellfun(@(x) x.n, v);
nFinite = cellfun(@(x) x.nFinite, v);
nUnique = cellfun(@(x) x.nUnique, v);
minimum = cellfun(@(x) x.min, v);
medianValue = cellfun(@(x) x.median, v);
maximum = cellfun(@(x) x.max, v);
meanValue = cellfun(@(x) x.mean, v);
t = table(names, n, nFinite, nUnique, minimum, medianValue, maximum, meanValue);
end

function write_report(filename, a)
fid = fopen(filename, 'w');
assert(fid >= 0, 'Could not open report file: %s', filename);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '# Stage-2 MMU unit audit\n\n');
fprintf(fid, '- Artifact: `%s`\n', a.artifactFile);
fprintf(fid, '- Number of spots: %d\n', a.nSpots);
fprintf(fid, '- Conversion: $MU_i=10^6w_i/N_i$.\n');
fprintf(fid, '- Inverse lower bound: $w_i^{min}=N_i MU_i^{min}/10^6$.\n');
fprintf(fid, '- Machine MMU source: `%s`.\n', a.machineMinMUSource);
fprintf(fid, '- Machine MMU verified nonzero: **%d**.\n', a.machineMMUVerified);
fprintf(fid, '- `cfg.mmuThreshold` in weight units: `%.12g`.\n', ...
    a.planningWeightThreshold);
fprintf(fid, '- Planning threshold equals recorded machine MMU: **%d**.\n\n', ...
    a.planningThresholdIsPhysicalMMU);
fprintf(fid, '| Quantity | min | median | max | unique |\n');
fprintf(fid, '|---|---:|---:|---:|---:|\n');
report_row(fid, 'dij.minMU (physical MU)', a.machineMinMU);
report_row(fid, 'numParticlesPerMU', a.numParticlesPerMU);
report_row(fid, 'machine lower bound (weight)', a.machineMinimumWeight);
report_row(fid, 'cfg threshold mapped to MU', a.planningThresholdPhysicalMU);
report_row(fid, 'Dense positive MU', a.densePositivePhysicalMU);
fprintf(fid, ['\nA nonzero ADMM MMU threshold is authorized only if the recorded machine ' ...
    'MMU is nonzero and its provenance is confirmed, or if a distinct planning ' ...
    'MMU in physical MU is explicitly adopted and documented.\n']);
end

function report_row(fid, name, s)
fprintf(fid, '| %s | %.12g | %.12g | %.12g | %d |\n', ...
    name, s.min, s.median, s.max, s.nUnique);
end

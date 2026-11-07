% test_integrated_surrogate_gradient  Finite-difference ES/time gradients.

rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(rootDir, 'integrated'));
addpath(fullfile(rootDir, 'config'));

groups = repmat(struct('ix', [], 'cpIndex', 0, 'energy', 0), 4, 1);
groups(1) = struct('ix', [1;2], 'cpIndex', 1, 'energy', 100);
groups(2) = struct('ix', 3, 'cpIndex', 1, 'energy', 120);
groups(3) = struct('ix', [4;5], 'cpIndex', 2, 'energy', 110);
groups(4) = struct('ix', 6, 'cpIndex', 2, 'energy', 130);
cpGroups = repmat(struct('ix', [], 'cpIndex', 0, 'gantryAngle', 0), 2, 1);
cpGroups(1) = struct('ix', [1;2;3], 'cpIndex', 1, 'gantryAngle', 0);
cpGroups(2) = struct('ix', [4;5;6], 'cpIndex', 2, 'gantryAngle', 20);
map = repmat(struct('globalIx', 0), 6, 1);
for i = 1:6, map(i).globalIx = i; end
cfg = default_config();
machine = default_machine_model();
context = build_integrated_context(groups, cpGroups, map, machine, cfg);
w = [0.4; 0.8; 0.2; 0.7; 0.3; 0.5];
[~, gradES, ~, gradT] = evaluate_delivery_surrogates(w, context);

h = 1e-6;
fdES = zeros(size(w));
fdT = zeros(size(w));
for i = 1:numel(w)
    wp = w; wm = w;
    wp(i) = wp(i) + h;
    wm(i) = wm(i) - h;
    [fpES, ~, fpT] = evaluate_delivery_surrogates(wp, context);
    [fmES, ~, fmT] = evaluate_delivery_surrogates(wm, context);
    fdES(i) = (fpES - fmES) / (2*h);
    fdT(i) = (fpT - fmT) / (2*h);
end
assert(norm(fdES - gradES) / max(1, norm(fdES)) < 2e-5, ...
    'Energy-sequencing gradient failed finite differences.');
assert(norm(fdT - gradT) / max(1, norm(fdT)) < 2e-5, ...
    'Delivery-surrogate gradient failed finite differences.');
fprintf('[OK] Integrated ES/time surrogate gradients validated.\n');

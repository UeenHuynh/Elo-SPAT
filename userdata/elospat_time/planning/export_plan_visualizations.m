function plotFiles = export_plan_visualizations(ct, cst, dij, resultGUI, meta, outputDir, plotPrefix)
% export_plan_visualizations  Save the main dose plot (DVH) for plan results.
%
% INPUTS
%   ct        - matRad CT struct (used for metadata only)
%   cst       - matRad structure table
%   dij       - matRad dose influence matrix
%   resultGUI - optimization result with weights and/or dose cube
%   meta      - metadata struct, optionally with meta.outputPaths.runDir
%   outputDir - optional destination folder
%   plotPrefix - optional filename prefix
%
% OUTPUT
%   plotFiles - struct listing saved PNG file paths

narginchk(5, 7);

if nargin < 6 || isempty(outputDir)
    outputDir = derive_output_dir(meta);
end

if nargin < 7 || isempty(plotPrefix)
    plotPrefix = 'plan';
end

if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

plotFiles = struct('dvh', '', 'doseAxial', '', 'doseCoronal', '');

[doseCube, cubeLabel] = resolve_dose_cube(dij, resultGUI);

plotFiles.dvh = save_dvh_plot(cst, doseCube, outputDir, plotPrefix);
plotFiles.doseAxial = '';
plotFiles.doseCoronal = '';

end

function outputDir = derive_output_dir(meta)
outputDir = fullfile('userdata', 'elospat_time', 'results');

if isstruct(meta) && isfield(meta, 'outputPaths') && isstruct(meta.outputPaths)
    if isfield(meta.outputPaths, 'runDir') && ~isempty(meta.outputPaths.runDir)
        outputDir = meta.outputPaths.runDir;
    end
end
end

function [doseCube, cubeLabel] = resolve_dose_cube(dij, resultGUI)
cubeLabel = 'Dose';
doseCube = [];

if isstruct(resultGUI)
    if isfield(resultGUI, 'RBExDose') && ~isempty(resultGUI.RBExDose)
        doseCube = resultGUI.RBExDose;
        cubeLabel = 'RBEx dose';
        return;
    end

    if isfield(resultGUI, 'physicalDose') && ~isempty(resultGUI.physicalDose)
        doseCube = resultGUI.physicalDose;
        cubeLabel = 'Physical dose';
        return;
    end
end

if ~isstruct(resultGUI) || ~isfield(resultGUI, 'w') || isempty(resultGUI.w)
    error('export_plan_visualizations:missingDose', ...
        'resultGUI must contain a dose cube or optimized weights.');
end

resultWithCubes = matRad_calcCubes(resultGUI.w, dij);

if isfield(resultWithCubes, 'RBExDose') && ~isempty(resultWithCubes.RBExDose)
    doseCube = resultWithCubes.RBExDose;
    cubeLabel = 'RBEx dose';
elseif isfield(resultWithCubes, 'physicalDose') && ~isempty(resultWithCubes.physicalDose)
    doseCube = resultWithCubes.physicalDose;
    cubeLabel = 'Physical dose';
else
    error('export_plan_visualizations:missingDose', ...
        'Could not derive a dose cube from resultGUI.w and dij.');
end

end

function plotPath = save_dvh_plot(cst, doseCube, outputDir, plotPrefix)
dvh = matRad_calcDVH(cst, doseCube, 'cum');
fig = figure('Visible', 'off', 'Color', 'w');
hold on;

for i = 1:numel(dvh)
    plot(dvh(i).doseGrid, dvh(i).volumePoints, 'LineWidth', 1.5);
end

legendLabels = cell(size(cst, 1), 1);
for i = 1:size(cst, 1)
    legendLabels{i} = stringify_label(cst{i, 2}, sprintf('Structure %d', i));
end

legend(legendLabels, 'Interpreter', 'none', 'Location', 'bestoutside');
xlabel('Dose');
ylabel('Volume (%)');
title('Dose-volume histogram');
grid on;
box on;

plotPath = fullfile(outputDir, sprintf('%s_dvh.png', plotPrefix));
save_figure_png(fig, plotPath);
close(fig);
end

function value = stringify_label(rawValue, fallbackValue)
if isstring(rawValue)
    value = char(rawValue);
elseif ischar(rawValue)
    value = rawValue;
else
    value = fallbackValue;
end

value = strtrim(value);
if isempty(value)
    value = fallbackValue;
end
end

function save_figure_png(fig, plotPath)
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotPath, '-dpng', '-r150');
end

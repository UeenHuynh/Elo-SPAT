function plotFiles = export_arc_visualizations(pln, stf, dij, outputDir, plotPrefix)
% export_arc_visualizations  Save the main arc-geometry plot for Task 6 outputs.
%
% INPUTS
%   pln       - matRad plan struct with pln.propStf.gantryAngles
%   stf       - matRad steering information struct array
%   dij       - matRad dose influence matrix
%   outputDir - optional destination folder
%   plotPrefix - optional filename prefix
%
% OUTPUT
%   plotFiles - struct listing saved PNG file paths

narginchk(3, 5);

if nargin < 4 || isempty(outputDir)
    outputDir = derive_output_dir(pln);
end

if nargin < 5 || isempty(plotPrefix)
    plotPrefix = 'arc';
end

if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

plotFiles = struct('angles', '', 'beamBixels', '');

plotFiles.angles = save_angle_plot(pln, dij, outputDir, plotPrefix);
plotFiles.beamBixels = '';

end

function outputDir = derive_output_dir(pln)
outputDir = fullfile('userdata', 'elospat_time', 'results');

if isstruct(pln) && isfield(pln, 'rmsPersist') && isstruct(pln.rmsPersist)
    if isfield(pln.rmsPersist, 'runDir') && ~isempty(pln.rmsPersist.runDir)
        outputDir = pln.rmsPersist.runDir;
    end
end
end

function plotPath = save_angle_plot(pln, dij, outputDir, plotPrefix)
angles = pln.propStf.gantryAngles(:);
beamIdx = 1:numel(angles);

fig = figure('Visible', 'off', 'Color', 'w');
plot(beamIdx, angles, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on;
box on;
xlabel('Beam index');
ylabel('Gantry angle (degrees)');
title(sprintf('Discrete arc geometry | %d beams | %d bixels', ...
    numel(angles), safe_total_bixels(dij)));

plotPath = fullfile(outputDir, sprintf('%s_gantry_angles.png', plotPrefix));
save_figure_png(fig, plotPath);
close(fig);
end

function value = safe_total_bixels(dij)
value = 0;
if isstruct(dij) && isfield(dij, 'totalNumOfBixels') && isnumeric(dij.totalNumOfBixels)
    value = dij.totalNumOfBixels;
end
end

function save_figure_png(fig, plotPath)
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotPath, '-dpng', '-r150');
end

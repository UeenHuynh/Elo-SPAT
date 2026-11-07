function plotPath = export_group_index_visualization(groups, cpGroups, outputDir, plotPrefix)
% export_group_index_visualization  Save a control-point by energy-layer occupancy heatmap.
%
% INPUTS
%   groups    - struct array from build_group_index()
%   cpGroups  - struct array from build_group_index()
%   outputDir - directory where PNG should be saved
%   plotPrefix - optional filename prefix
%
% OUTPUT
%   plotPath  - saved PNG path

narginchk(3, 4);

if nargin < 4 || isempty(plotPrefix)
    plotPrefix = 'task9_groups';
end

if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

validate_groups(groups, cpGroups);

cpValues = reshape(double([cpGroups.cpIndex]), [], 1);
groupCp = reshape(double([groups.cpIndex]), [], 1);
groupEnergyIndex = reshape(double([groups.energyIndex]), [], 1);
groupSpotCount = reshape(double([groups.nSpots]), [], 1);

nCp = numel(cpValues);
nEnergy = max(groupEnergyIndex);
occupancy = zeros(nCp, nEnergy);

for i = 1:numel(groups)
    cpRow = find(cpValues == groupCp(i), 1, 'first');
    occupancy(cpRow, groupEnergyIndex(i)) = groupSpotCount(i);
end

fig = figure('Visible', 'off', 'Color', 'w');
imagesc(occupancy);
axis tight;
colorbar;
colormap(parula);
xlabel('Energy-layer index');
ylabel('Control-point index');
title('Control-point by energy-layer occupancy');

plotPath = fullfile(outputDir, sprintf('%s_cp_energy_occupancy.png', plotPrefix));
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotPath, '-dpng', '-r150');
close(fig);

end

function validate_groups(groups, cpGroups)
if ~isstruct(groups) || isempty(groups)
    error('export_group_index_visualization:invalidInput', ...
        'groups must be a non-empty struct array.');
end

if ~isstruct(cpGroups) || isempty(cpGroups)
    error('export_group_index_visualization:invalidInput', ...
        'cpGroups must be a non-empty struct array.');
end

requiredGroupFields = {'cpIndex', 'energyIndex', 'nSpots'};
requiredCpFields = {'cpIndex'};

for i = 1:numel(requiredGroupFields)
    if ~isfield(groups, requiredGroupFields{i})
        error('export_group_index_visualization:missingField', ...
            'groups.%s is required.', requiredGroupFields{i});
    end
end

for i = 1:numel(requiredCpFields)
    if ~isfield(cpGroups, requiredCpFields{i})
        error('export_group_index_visualization:missingField', ...
            'cpGroups.%s is required.', requiredCpFields{i});
    end
end
end

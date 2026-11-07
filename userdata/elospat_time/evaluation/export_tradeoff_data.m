function outputs = export_tradeoff_data(inputData, outputPrefix, outputDir, doseMetricColumn)
% export_tradeoff_data  Save dose-vs-estimated-BDT pairs for Pareto plotting.
%
% USAGE
%   outputs = export_tradeoff_data(methodTable)
%   outputs = export_tradeoff_data(methodEntries, outputPrefix, outputDir)
%
% OUTPUTS
%   outputs.tradeoffTable
%   outputs.csvFile
%   outputs.matFile

narginchk(1, 4);

methodTable = normalize_input_table(inputData);

if nargin < 2 || isempty(outputPrefix)
    outputPrefix = 'tradeoff_data';
end
if nargin < 3 || isempty(outputDir)
    outputDir = pwd;
end
if nargin < 4 || isempty(doseMetricColumn)
    doseMetricColumn = 'dose_targetD95Min';
end

if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

requiredColumns = {'Method', doseMetricColumn, 'delivery_T_total'};
missing = requiredColumns(~ismember(requiredColumns, methodTable.Properties.VariableNames));
if ~isempty(missing)
    error('export_tradeoff_data:missingColumns', ...
        'Missing required columns: %s', strjoin(missing, ', '));
end

tradeoffTable = table();
tradeoffTable.Method = methodTable.Method;
tradeoffTable.doseMetricName = repmat(string(doseMetricColumn), height(methodTable), 1);
tradeoffTable.doseMetricValue = methodTable.(doseMetricColumn);
tradeoffTable.estimatedBDT = methodTable.delivery_T_total;

optionalColumns = {'dose_oarDmaxMax', 'delivery_nActiveLayers', 'delivery_nActiveSpots', 'solver_iterations'};
for i = 1:numel(optionalColumns)
    if ismember(optionalColumns{i}, methodTable.Properties.VariableNames)
        tradeoffTable.(optionalColumns{i}) = methodTable.(optionalColumns{i});
    end
end

csvFile = fullfile(outputDir, sprintf('%s.csv', outputPrefix));
matFile = fullfile(outputDir, sprintf('%s.mat', outputPrefix));

writetable(tradeoffTable, csvFile);
save(matFile, 'tradeoffTable', 'doseMetricColumn', '-v7.3');

outputs = struct();
outputs.tradeoffTable = tradeoffTable;
outputs.csvFile = csvFile;
outputs.matFile = matFile;
end

function methodTable = normalize_input_table(inputData)
if istable(inputData)
    methodTable = inputData;
elseif isstruct(inputData)
    methodTable = compare_methods_table(inputData);
else
    error('export_tradeoff_data:invalidInput', ...
        'Input must be a table or method-entry struct array.');
end
end

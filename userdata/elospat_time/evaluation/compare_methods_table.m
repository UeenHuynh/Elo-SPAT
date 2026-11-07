function methodTable = compare_methods_table(methodEntries)
% compare_methods_table  Build a flat comparison table across methods.
%
% INPUT
%   methodEntries - struct array with fields:
%                   .name or .methodName
%                   .doseMetrics
%                   .deliveryMetrics
%                   .solverMetrics
%
% OUTPUT
%   methodTable - MATLAB table with one row per method and prefixed metric
%                 columns:
%                 dose_*, delivery_*, solver_*

narginchk(1, 1);

if ~isstruct(methodEntries) || isempty(methodEntries)
    error('compare_methods_table:invalidInput', ...
        'methodEntries must be a non-empty struct array.');
end

nMethods = numel(methodEntries);
methodNames = strings(nMethods, 1);
doseRows = cell(nMethods, 1);
deliveryRows = cell(nMethods, 1);
solverRows = cell(nMethods, 1);

allDoseFields = strings(0, 1);
allDeliveryFields = strings(0, 1);
allSolverFields = strings(0, 1);

for i = 1:nMethods
    methodNames(i) = string(resolve_method_name(methodEntries(i), i));

    doseRows{i} = flatten_metric_struct(resolve_metric_substruct(methodEntries(i), 'doseMetrics'));
    deliveryRows{i} = flatten_metric_struct(resolve_metric_substruct(methodEntries(i), 'deliveryMetrics'));
    solverRows{i} = flatten_metric_struct(resolve_metric_substruct(methodEntries(i), 'solverMetrics'));

    allDoseFields = union(allDoseFields, string(fieldnames(doseRows{i})));
    allDeliveryFields = union(allDeliveryFields, string(fieldnames(deliveryRows{i})));
    allSolverFields = union(allSolverFields, string(fieldnames(solverRows{i})));
end

tableStruct = struct();
tableStruct.Method = methodNames;
tableStruct = append_metric_block(tableStruct, doseRows, allDoseFields, 'dose_');
tableStruct = append_metric_block(tableStruct, deliveryRows, allDeliveryFields, 'delivery_');
tableStruct = append_metric_block(tableStruct, solverRows, allSolverFields, 'solver_');

methodTable = struct2table(tableStruct);
end

function name = resolve_method_name(entry, defaultIndex)
if isfield(entry, 'name') && ~isempty(entry.name)
    name = char(string(entry.name));
elseif isfield(entry, 'methodName') && ~isempty(entry.methodName)
    name = char(string(entry.methodName));
else
    name = sprintf('Method %d', defaultIndex);
end
end

function metricStruct = resolve_metric_substruct(entry, fieldName)
if ~isfield(entry, fieldName) || isempty(entry.(fieldName))
    metricStruct = struct();
    return;
end
metricStruct = entry.(fieldName);
end

function flat = flatten_metric_struct(metricStruct)
flat = struct();
if ~isstruct(metricStruct) || isempty(fieldnames(metricStruct))
    return;
end

if isfield(metricStruct, 'summary') && isstruct(metricStruct.summary)
    metricStruct = metricStruct.summary;
end

fields = fieldnames(metricStruct);
for i = 1:numel(fields)
    value = metricStruct.(fields{i});
    if isscalar(value) && (isnumeric(value) || islogical(value))
        flat.(fields{i}) = value;
    elseif ischar(value)
        flat.(fields{i}) = string(value);
    elseif isstring(value) && isscalar(value)
        flat.(fields{i}) = value;
    end
end
end

function tableStruct = append_metric_block(tableStruct, rows, fieldNames, prefix)
nRows = numel(rows);
for i = 1:numel(fieldNames)
    rawName = char(fieldNames(i));
    columnName = matlab.lang.makeValidName([prefix rawName]);
    sampleValue = find_first_value(rows, rawName);

    if isstring(sampleValue) || ischar(sampleValue)
        column = strings(nRows, 1);
        for r = 1:nRows
            if isfield(rows{r}, rawName)
                column(r) = string(rows{r}.(rawName));
            else
                column(r) = "";
            end
        end
    else
        column = NaN(nRows, 1);
        for r = 1:nRows
            if isfield(rows{r}, rawName)
                value = rows{r}.(rawName);
                if islogical(value)
                    column(r) = double(value);
                else
                    column(r) = double(value);
                end
            end
        end
    end

    tableStruct.(columnName) = column;
end
end

function value = find_first_value(rows, rawName)
value = NaN;
for i = 1:numel(rows)
    if isfield(rows{i}, rawName)
        value = rows{i}.(rawName);
        return;
    end
end
end

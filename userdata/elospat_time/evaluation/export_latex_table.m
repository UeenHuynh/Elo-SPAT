function [latexText, outputFile] = export_latex_table(inputData, outputFile)
% export_latex_table  Export a compact comparison table as LaTeX tabular.
%
% USAGE
%   latexText = export_latex_table(methodTable)
%   latexText = export_latex_table(methodEntries)
%   [latexText, outputFile] = export_latex_table(..., outputFile)
%
% Columns:
%   Method, CTV D95, CTV D98, OAR Dmax, Active EL, Active Spots,
%   Switch-up, T_energy, T_spot, T_gantry, T_total

narginchk(1, 2);

methodTable = normalize_input_table(inputData);
researchMeta = resolve_research_meta();

requiredColumns = { ...
    'Method', ...
    'dose_targetD95Min', ...
    'dose_targetD98Min', ...
    'dose_oarDmaxMax', ...
    'delivery_nActiveLayers', ...
    'delivery_nActiveSpots', ...
    'delivery_nSwitchUp', ...
    'delivery_T_energy', ...
    'delivery_T_spot', ...
    'delivery_T_gantry', ...
    'delivery_T_total'};

missing = requiredColumns(~ismember(requiredColumns, methodTable.Properties.VariableNames));
if ~isempty(missing)
    error('export_latex_table:missingColumns', ...
        'Missing required columns: %s', strjoin(missing, ', '));
end

header = [ ...
    '%% ' researchMeta.disclaimer newline ...
    '%% Machine model note: ' researchMeta.machineModelNote newline ...
    '%% Citation placeholder: ' researchMeta.citationPlaceholder newline ...
    '\begin{tabular}{lrrrrrrrrrr}' newline ...
    '\hline' newline ...
    'Method & CTV D95 & CTV D98 & OAR Dmax & Active EL & Active Spots & Switch-up & T\_energy & T\_spot & T\_gantry & T\_total \\' newline ...
    '\hline' newline];

bodyLines = strings(height(methodTable), 1);
for i = 1:height(methodTable)
    bodyLines(i) = sprintf([ ...
        '%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\'], ...
        escape_latex_text(string(methodTable.Method(i))), ...
        fmt_num(methodTable.dose_targetD95Min(i), 2), ...
        fmt_num(methodTable.dose_targetD98Min(i), 2), ...
        fmt_num(methodTable.dose_oarDmaxMax(i), 2), ...
        fmt_num(methodTable.delivery_nActiveLayers(i), 0), ...
        fmt_num(methodTable.delivery_nActiveSpots(i), 0), ...
        fmt_num(methodTable.delivery_nSwitchUp(i), 0), ...
        fmt_num(methodTable.delivery_T_energy(i), 2), ...
        fmt_num(methodTable.delivery_T_spot(i), 2), ...
        fmt_num(methodTable.delivery_T_gantry(i), 2), ...
        fmt_num(methodTable.delivery_T_total(i), 2));
end

footer = [newline '\hline' newline '\end{tabular}' newline];
latexText = [header strjoin(cellstr(bodyLines), newline) footer];

if nargin >= 2 && ~isempty(outputFile)
    fid = fopen(outputFile, 'w');
    if fid == -1
        error('export_latex_table:fileOpenFailed', ...
            'Could not open %s for writing.', outputFile);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s', latexText);
else
    outputFile = '';
end
end

function methodTable = normalize_input_table(inputData)
if istable(inputData)
    methodTable = inputData;
elseif isstruct(inputData)
    methodTable = compare_methods_table(inputData);
else
    error('export_latex_table:invalidInput', ...
        'Input must be a table or method-entry struct array.');
end
end

function researchMeta = resolve_research_meta()
if exist('research_only_metadata', 'file') == 2
    researchMeta = research_only_metadata();
else
    researchMeta = struct( ...
        'disclaimer', 'Research-only output.', ...
        'machineModelNote', 'Machine timing estimates are configurable proxies.', ...
        'citationPlaceholder', 'TODO: add citation.');
end
end

function text = fmt_num(value, decimals)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    text = '--';
    return;
end

if decimals == 0
    text = sprintf('%.0f', value);
else
    text = sprintf(['%.' num2str(decimals) 'f'], value);
end
end

function text = escape_latex_text(text)
text = replace(text, "\", '\textbackslash{}');
text = replace(text, "_", '\_');
text = replace(text, "%", '\%');
text = replace(text, "&", '\&');
text = char(text);
end

function plotFiles = replot_saved_result(resultMatPath, patientFile, outputDir)
% replot_saved_result  Rebuild PNG plots from a saved Task 5/6/7 MAT file.
%
% INPUTS
%   resultMatPath - path to a MAT artifact containing plan data
%   patientFile   - optional patient MAT file with ct and cst
%   outputDir     - optional output folder for regenerated plots
%
% OUTPUT
%   plotFiles - struct of saved plot paths

narginchk(1, 3);

if nargin < 2
    patientFile = '';
end

if nargin < 3 || isempty(outputDir)
    outputDir = fileparts(resultMatPath);
end

if exist(resultMatPath, 'file') ~= 2
    error('replot_saved_result:fileNotFound', ...
        'Result MAT file not found: %s', resultMatPath);
end

data = load(resultMatPath);
[~, plotPrefix] = fileparts(resultMatPath);

if isfield(data, 'resultGUI') || isfield(data, 'result')
    [ct, cst] = load_patient_context(data, patientFile);
    resultStruct = choose_result_struct(data);
    meta = choose_meta_struct(data);
    
    if ~isfield(data, 'dij')
        error('replot_saved_result:missingField', ...
            'Saved result is missing dij, which is required for dose plotting.');
    end
    
    plotFiles = export_plan_visualizations( ...
        ct, cst, data.dij, resultStruct, meta, outputDir, plotPrefix);
    print_plot_summary(plotFiles, outputDir);
    return;
end

if isfield(data, 'pln') && isfield(data, 'stf') && isfield(data, 'dij')
    plotFiles = export_arc_visualizations(data.pln, data.stf, data.dij, outputDir, plotPrefix);
    print_plot_summary(plotFiles, outputDir);
    return;
end

error('replot_saved_result:unsupportedArtifact', ...
    'Saved MAT file does not contain a supported Task 5/6/7 artifact.');

end

function resultStruct = choose_result_struct(data)
if isfield(data, 'resultGUI')
    resultStruct = data.resultGUI;
elseif isfield(data, 'result')
    resultStruct = data.result;
else
    error('replot_saved_result:missingField', ...
        'Saved result file does not contain resultGUI or result.');
end
end

function meta = choose_meta_struct(data)
if isfield(data, 'meta') && isstruct(data.meta)
    meta = data.meta;
else
    meta = struct();
end
end

function [ct, cst] = load_patient_context(data, patientFile)
if isfield(data, 'ct') && isfield(data, 'cst')
    ct = data.ct;
    cst = data.cst;
    return;
end

resolvedPatientFile = patientFile;
if isempty(resolvedPatientFile)
    resolvedPatientFile = which('PROSTATE.mat');
end

if isempty(resolvedPatientFile) || exist(resolvedPatientFile, 'file') ~= 2
    error('replot_saved_result:missingPatientContext', ...
        'Patient MAT file with ct and cst is required to redraw dose plots.');
end

patientData = load(resolvedPatientFile, 'ct', 'cst');

if ~isfield(patientData, 'ct') || ~isfield(patientData, 'cst')
    error('replot_saved_result:invalidPatientContext', ...
        'Patient file does not contain ct and cst variables.');
end

ct = patientData.ct;
cst = patientData.cst;
end

function print_plot_summary(plotFiles, outputDir)
fprintf('Saved plots to: %s\n', outputDir);

plotNames = fieldnames(plotFiles);
for i = 1:numel(plotNames)
    fieldName = plotNames{i};
    plotPath = plotFiles.(fieldName);
    if ischar(plotPath) && ~isempty(plotPath)
        fprintf('  %s: %s\n', fieldName, plotPath);
    end
end
end

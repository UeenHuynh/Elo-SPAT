function metrics = compute_solver_metrics(resultLike, meta, optimizer)
% compute_solver_metrics  Extract solver summary metrics from matRad outputs.
%
% USAGE
%   metrics = compute_solver_metrics(resultLike)
%   metrics = compute_solver_metrics(resultLike, meta)
%   metrics = compute_solver_metrics(resultLike, meta, optimizer)
%
% OUTPUT
%   metrics - struct with fields:
%             .finalObjective
%             .doseComponent
%             .regularizationComponent
%             .iterations
%             .runtimeSeconds
%             .converged
%             .exitStatus
%             .statusMessage
%             .optimizerName
%             .acceptableButNotConverged

narginchk(1, 3);

if nargin < 2 || isempty(meta)
    meta = struct();
end
if nargin < 3
    optimizer = [];
end

infoStruct = extract_info_struct(resultLike, meta, optimizer);
resultInfoStruct = extract_result_info_struct(resultLike, meta, optimizer);

metrics = struct();
metrics.finalObjective = first_finite_scalar( ...
    read_any(infoStruct, {'objective'}), ...
    read_any(resultLike, {'fval', 'objective'}), ...
    read_any(resultInfoStruct, {'objective'}), ...
    NaN);
metrics.doseComponent = first_finite_scalar( ...
    read_any(infoStruct, {'doseComponent', 'objectiveDose', 'doseTerm'}), ...
    read_any(resultInfoStruct, {'doseComponent', 'objectiveDose', 'doseTerm'}), ...
    NaN);
metrics.regularizationComponent = first_finite_scalar( ...
    read_any(infoStruct, {'regularizationComponent', 'regComponent', 'penaltyTerm'}), ...
    read_any(resultInfoStruct, {'regularizationComponent', 'regComponent', 'penaltyTerm'}), ...
    NaN);
metrics.iterations = first_finite_scalar( ...
    read_any(resultInfoStruct, {'iterations', 'iter'}), ...
    read_any(infoStruct, {'iterations', 'iter'}), ...
    NaN);
metrics.runtimeSeconds = first_finite_scalar( ...
    read_any(meta, {'runtime', 'runtimeSeconds'}), ...
    read_any(resultInfoStruct, {'cpu', 'runtime', 'runtimeSeconds'}), ...
    read_any(infoStruct, {'cpu', 'runtime', 'runtimeSeconds'}), ...
    NaN);
metrics.exitStatus = first_finite_scalar( ...
    read_any(meta, {'exitStatus', 'exitFlag'}), ...
    read_any(resultLike, {'exitStatus', 'exitFlag'}), ...
    read_any(resultInfoStruct, {'status', 'exitStatus', 'exitFlag'}), ...
    read_any(infoStruct, {'status', 'exitStatus', 'exitFlag'}), ...
    NaN);
metrics.converged = first_logical_value( ...
    read_any(meta, {'converged'}), ...
    read_any(resultLike, {'converged'}), ...
    infer_convergence_from_status(metrics.exitStatus), ...
    false);
metrics.statusMessage = first_text_value( ...
    read_any(meta, {'statusMessage', 'exitMessage', 'exitMsg'}), ...
    read_any(resultLike, {'statusMessage', 'exitMessage', 'exitMsg'}), ...
    status_text_from_code(metrics.exitStatus), ...
    '');
metrics.optimizerName = first_text_value( ...
    read_any(meta, {'optimizerName'}), ...
    read_any(resultLike, {'usedOptimizer'}), ...
    class_name_or_empty(optimizer), ...
    '');
metrics.acceptableButNotConverged = (~metrics.converged) && isequal(metrics.exitStatus, -1);

end

function infoStruct = extract_info_struct(resultLike, meta, optimizer)
infoStruct = struct();

candidate = read_any(resultLike, {'info'});
if isstruct(candidate) && ~isempty(fieldnames(candidate))
    infoStruct = candidate;
    return;
end

candidate = read_any(meta, {'info'});
if isstruct(candidate) && ~isempty(fieldnames(candidate))
    infoStruct = candidate;
    return;
end

candidate = read_any(optimizer, {'info'});
if isstruct(candidate) && ~isempty(fieldnames(candidate))
    infoStruct = candidate;
end
end

function resultInfoStruct = extract_result_info_struct(resultLike, meta, optimizer)
resultInfoStruct = struct();

candidate = read_any(meta, {'optimizerResultInfo', 'resultInfo'});
if isstruct(candidate) && ~isempty(fieldnames(candidate))
    resultInfoStruct = candidate;
    return;
end

candidate = read_any(optimizer, {'resultInfo'});
if isstruct(candidate) && ~isempty(fieldnames(candidate))
    resultInfoStruct = candidate;
    return;
end

candidate = read_any(resultLike, {'resultInfo'});
if isstruct(candidate) && ~isempty(fieldnames(candidate))
    resultInfoStruct = candidate;
end
end

function value = read_any(source, names)
value = [];
for i = 1:numel(names)
    if has_field_or_prop(source, names{i})
        value = read_field_or_prop(source, names{i});
        if ~isempty(value)
            return;
        end
    end
end
end

function tf = has_field_or_prop(source, name)
tf = false;
if isstruct(source)
    tf = isfield(source, name);
elseif isobject(source)
    tf = isprop(source, name);
end
end

function value = read_field_or_prop(source, name)
if isstruct(source)
    value = source.(name);
elseif isobject(source) && isprop(source, name)
    value = source.(name);
else
    value = [];
end
end

function value = first_finite_scalar(varargin)
fallback = varargin{end};
for i = 1:(nargin - 1)
    candidate = varargin{i};
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
        value = double(candidate);
        return;
    end
end
value = fallback;
end

function value = first_logical_value(varargin)
fallback = varargin{end};
for i = 1:(nargin - 1)
    candidate = varargin{i};
    if islogical(candidate) && isscalar(candidate)
        value = candidate;
        return;
    end
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
        value = logical(candidate);
        return;
    end
end
value = fallback;
end

function value = first_text_value(varargin)
fallback = varargin{end};
for i = 1:(nargin - 1)
    candidate = varargin{i};
    if ischar(candidate)
        candidate = strtrim(candidate);
        if ~isempty(candidate)
            value = candidate;
            return;
        end
    elseif isstring(candidate) && isscalar(candidate)
        candidate = strtrim(char(candidate));
        if ~isempty(candidate)
            value = candidate;
            return;
        end
    end
end
value = fallback;
end

function converged = infer_convergence_from_status(status)
converged = false;
if isnumeric(status) && isscalar(status) && isfinite(status)
    converged = any(status == [0, 1]);
end
end

function text = status_text_from_code(status)
if ~isnumeric(status) || ~isscalar(status) || ~isfinite(status)
    text = '';
    return;
end

switch double(status)
    case 0
        text = 'Solved';
    case 1
        text = 'Solved To Acceptable Level';
    case -1
        text = 'Maximum Number of Iterations Exceeded';
    otherwise
        text = sprintf('Status %g', double(status));
end
end

function value = class_name_or_empty(source)
if isempty(source)
    value = '';
elseif isobject(source)
    value = class(source);
else
    value = '';
end
end

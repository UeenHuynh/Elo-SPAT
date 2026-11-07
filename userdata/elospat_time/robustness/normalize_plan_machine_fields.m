function pln = normalize_plan_machine_fields(pln)
% normalize_plan_machine_fields  Normalize pln text fields used by matRad machine lookup.

narginchk(1, 1);

if ~isstruct(pln)
    error('normalize_plan_machine_fields:invalidInput', ...
        'pln must be a struct.');
end

if ~isfield(pln, 'radiationMode') || isempty(pln.radiationMode)
    pln.radiationMode = 'protons';
else
    pln.radiationMode = normalize_text_scalar(pln.radiationMode, 'pln.radiationMode');
end

if ~isfield(pln, 'machine') || isempty(pln.machine)
    pln.machine = 'Generic';
else
    pln.machine = normalize_text_scalar(pln.machine, 'pln.machine');
end
end

function value = normalize_text_scalar(value, fieldName)
if isstring(value)
    if ~isscalar(value)
        error('normalize_plan_machine_fields:invalidText', ...
            '%s must be a scalar string or character vector.', fieldName);
    end
    value = char(value);
elseif iscellstr(value)
    if ~isscalar(value)
        error('normalize_plan_machine_fields:invalidText', ...
            '%s must be a scalar string or character vector.', fieldName);
    end
    value = value{1};
elseif ~ischar(value)
    error('normalize_plan_machine_fields:invalidText', ...
        '%s must be a scalar string or character vector.', fieldName);
end

if isempty(value)
    error('normalize_plan_machine_fields:invalidText', ...
        '%s must not be empty.', fieldName);
end

if size(value, 1) > 1
    rows = strtrim(cellstr(value));
    rows = rows(~cellfun(@isempty, rows));
    if isempty(rows)
        error('normalize_plan_machine_fields:invalidText', ...
            '%s must not be empty.', fieldName);
    end

    uniqueRows = unique(rows, 'stable');
    if numel(uniqueRows) ~= 1
        error('normalize_plan_machine_fields:invalidText', ...
            '%s must resolve to a single scalar text value.', fieldName);
    end

    value = uniqueRows{1};
else
    value = strtrim(value);
end

if isempty(value)
    error('normalize_plan_machine_fields:invalidText', ...
        '%s must not be empty.', fieldName);
end
end

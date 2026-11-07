function orderedEnergies = order_layers_within_control_point(energies, mode, customOrder)
% order_layers_within_control_point  Order active energy layers within one CP.
%
% INPUTS
%   energies    - vector of active energy values at one control point
%   mode        - 'descending' (default) | 'ascending' | 'custom'
%   customOrder - optional vector defining the desired custom energy order
%
% OUTPUT
%   orderedEnergies - ordered energy vector with the same shape as energies

narginchk(1, 3);

if nargin < 2 || isempty(mode)
    mode = 'descending';
end

if nargin < 3
    customOrder = [];
end

if ~isnumeric(energies) || isempty(energies) || ~isvector(energies) || any(~isfinite(energies))
    error('order_layers_within_control_point:invalidInput', ...
        'energies must be a non-empty finite numeric vector.');
end

inputShape = size(energies);
energies = reshape(double(energies), [], 1);

switch lower(char(mode))
    case 'descending'
        orderedEnergies = sort(energies, 'descend');
    case 'ascending'
        orderedEnergies = sort(energies, 'ascend');
    case 'custom'
        if ~isnumeric(customOrder) || isempty(customOrder) || ~isvector(customOrder) || ...
                any(~isfinite(customOrder))
            error('order_layers_within_control_point:missingCustomOrder', ...
                'customOrder must be a non-empty finite numeric vector for mode ''custom''.');
        end

        customOrder = reshape(double(customOrder), [], 1);
        orderIx = zeros(numel(energies), 1);

        for i = 1:numel(energies)
            tolerance = max(1e-6, abs(energies(i)) * 1e-8);
            match = find(abs(customOrder - energies(i)) <= tolerance, 1, 'first');
            if isempty(match)
                error('order_layers_within_control_point:customOrderMismatch', ...
                    'Energy %.8g is not present in the supplied custom order.', ...
                    energies(i));
            end
            orderIx(i) = match;
        end

        [~, sortIx] = sort(orderIx, 'ascend');
        orderedEnergies = energies(sortIx);
    otherwise
        error('order_layers_within_control_point:invalidMode', ...
            'Unsupported mode ''%s''.', char(mode));
end

orderedEnergies = reshape(orderedEnergies, inputShape);

end

function [resultGUI, optimizer] = matRad_fluenceOptimization(dij, ~, ~, varargin)
% Test stub for wrapper-stage unit tests.

if nargin >= 4 && ~isempty(varargin{1})
    w0 = reshape(double(varargin{1}), [], 1);
else
    w0 = ones(double(dij.totalNumOfBixels), 1);
end

resultGUI = struct();
resultGUI.w = max(w0, 0);
resultGUI.fval = sum(resultGUI.w);

optimizer = struct();
optimizer.resultInfo = struct('status', 0, 'message', 'stub solved');
end

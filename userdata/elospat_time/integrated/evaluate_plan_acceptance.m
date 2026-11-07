function report = evaluate_plan_acceptance(doseQualityPass, ...
    maxConstraintViolation, supportPreserved, finalIpoptAcceptable, rule)
% evaluate_plan_acceptance  Clinical/feasibility gate, separate from stopping.

narginchk(4, 5);
if nargin < 5 || isempty(rule), rule = default_stabilization_rule(); end
report = struct();
report.doseQualityPass = logical(doseQualityPass);
report.constraintPass = isfinite(maxConstraintViolation) && ...
    maxConstraintViolation <= rule.constraintTolerance;
report.supportPreserved = logical(supportPreserved);
report.finalIpoptAcceptable = logical(finalIpoptAcceptable);
report.plan_accepted = report.doseQualityPass && report.constraintPass && ...
    report.supportPreserved && report.finalIpoptAcceptable;
end

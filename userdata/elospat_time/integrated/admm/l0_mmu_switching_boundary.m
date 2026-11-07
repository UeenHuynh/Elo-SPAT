function boundary = l0_mmu_switching_boundary(lambda0, rho, minimumWeight)
%L0_MMU_SWITCHING_BOUNDARY Exact a-value where off/on costs are equal.

m = reshape(double(minimumWeight), [], 1);
n = numel(m);
lambda0 = expand(lambda0, n);
assert(isscalar(rho) && isfinite(rho) && rho > 0);
threshold = sqrt(2 * lambda0 / rho);
boundary = threshold;
useSubMmuBranch = m > 0 & threshold < m;
boundary(useSubMmuBranch) = m(useSubMmuBranch) / 2 + ...
    lambda0(useSubMmuBranch) ./ (rho * m(useSubMmuBranch));
end

function x = expand(x, n)
x = reshape(double(x), [], 1);
if isscalar(x), x = repmat(x, n, 1); end
assert(numel(x) == n && all(isfinite(x)) && all(x >= 0));
end

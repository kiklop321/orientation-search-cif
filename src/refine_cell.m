function [cell, info] = refine_cell(cell0, hkl, tth_meas, wavelength)
%REFINE_CELL  Lattice parameters fitted to the 2theta of indexed spots.
%
%   [cell, info] = refine_cell(cell0, hkl, tth_meas, wavelength)
%
%   cell0      : starting cell [a b c alpha beta gamma] (A, deg), e.g. from the CIF
%   hkl        : N x 3 Miller indices of the matched spots
%   tth_meas   : N x 1 measured 2theta (deg)
%   wavelength : (A)
%
%   2theta of a reflection depends only on its d-spacing, i.e. only on the
%   cell and not on the orientation, so the cell can be fitted to the
%   measured 2theta alone:  minimise sum (2theta_calc(hkl, cell) - 2theta_meas)^2.
%   The free parameters follow from the cell of the CIF (angles are kept):
%     a = b = c                -> 1 parameter (cubic, rhombohedral setting)
%     a = b ~= c or a = c, ... -> 2 parameters (tetragonal, hexagonal)
%     a, b, c all different    -> 3 parameters (orthorhombic and lower)
%
%   info.n_par, info.rms_before, info.rms_after (deg), info.esd (A)
%   (standard deviations of a, b, c from the residuals and the Jacobian)

L = cell0(1:3);
tol = 1e-4 * max(L);
% which lengths are tied together
if abs(L(1) - L(2)) < tol && abs(L(2) - L(3)) < tol
    groups = {1:3};
elseif abs(L(1) - L(2)) < tol
    groups = {[1 2], 3};
elseif abs(L(1) - L(3)) < tol
    groups = {[1 3], 2};
elseif abs(L(2) - L(3)) < tol
    groups = {1, [2 3]};
else
    groups = {1, 2, 3};
end
np = numel(groups);
expand = @(p) assign(L, groups, p);
p0 = cellfun(@(g) L(g(1)), groups);

resid = @(p) tth_of(hkl, [expand(p), cell0(4:6)], wavelength) - tth_meas;
f = @(p) sum(resid(p).^2);
opts = optimset('TolX', 1e-8, 'TolFun', 1e-12, 'MaxFunEvals', 2000, 'Display', 'off');
p = fminsearch(f, p0, opts);
cell = [expand(p), cell0(4:6)];

r = resid(p);
info.n_par = np;
info.rms_before = sqrt(mean(resid(p0).^2));
info.rms_after  = sqrt(mean(r.^2));
% standard deviations from the numerical Jacobian
J = zeros(numel(r), np);
for k = 1:np
    dp = zeros(1, np);  dp(k) = 1e-6 * p(k);
    J(:, k) = (resid(p + dp) - r) / dp(k);
end
s2 = sum(r.^2) / max(numel(r) - np, 1);
esd_p = sqrt(diag(s2 * inv(J' * J)))';
info.esd = expand(esd_p);
end

function L = assign(L, groups, p)
for k = 1:numel(groups), L(groups{k}) = p(k); end
end

function tth = tth_of(hkl, cell, wavelength)
g = calculate_reciprocal_lattice(struct('a', cell(1:3), 'u', cell(4:6)), hkl, [1 0 0 0]);
s = vecnorm(g, 2, 2) * wavelength / (4 * pi);          % sin(theta)
tth = 2 * asind(min(s, 1));
end

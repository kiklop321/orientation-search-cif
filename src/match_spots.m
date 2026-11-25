function [T, calc] = match_spots(q, phi0, g, hkl, variant, expt, S)
%MATCH_SPOTS  Predicted spots of an orientation and their measured partners.
%
%   [T, calc] = match_spots(q, phi0, g, hkl, variant, expt, S)
%
%   q, phi0      : orientation from SEARCH_ORIENTATION
%   g, hkl, variant : from PHASE_REFLECTIONS
%   expt         : measured spots [2theta, beta] (deg)
%   S            : settings (wavelength, rotation_range, tth_range,
%                  tolerance(2))
%
%   calc : all predicted spots recorded during the rotation, columns
%          [2theta, beta, phi, h, k, l, variant]
%   T    : table with one row per predicted spot whose NEAREST measured spot
%          lies within the tolerance (S.tolerance(2) in both 2theta and beta):
%          hkl, variant, rotation angle phi - phi0,
%          predicted and measured 2theta and beta, their differences and
%          the index of the measured spot (row of expt)

[tth, b1, b2, p1, p2, missed] = ewald_angles(quatrotate(quatconj(q), g), S.wavelength);
h = hkl(~missed, :);  v = variant(~missed);
calc = [tth, b1, p1, h, v; tth, b2, p2, h, v];

width = S.rotation_range(2) - S.rotation_range(1);
keep = mod(calc(:,3) - phi0 - S.rotation_range(1), 360) <= width & ...
       calc(:,1) >= S.tth_range(1) & calc(:,1) <= S.tth_range(2);
calc = calc(keep, :);

% nearest measured spot (largest of |d2theta| and |dbeta|), accepted within the tolerance.
% beta runs from -90 to 270 deg and -90 = 270 (straight below the beam): spots
% within 5 deg of this seam are added once more with beta -+ 360
near_top = find(expt(:,2) > 265);  near_bottom = find(expt(:,2) < -85);
expt_w = [expt; expt(near_top, :) - [0 360]; expt(near_bottom, :) + [0 360]];
row_w  = [(1:size(expt, 1))'; near_top; near_bottom];
[iw, dist] = knnsearch(expt_w, calc(:,1:2), 'Distance', 'chebychev');
idx = row_w(iw);
hit = dist <= S.tolerance(2);
c = calc(hit, :);
e = expt(idx(hit), :);
rot = mod(c(:,3) - phi0 + 180, 360) - 180;
T = table(c(:,4), c(:,5), c(:,6), c(:,7), rot, c(:,1), e(:,1), c(:,1) - e(:,1), ...
    c(:,2), e(:,2), mod(c(:,2) - e(:,2) + 180, 360) - 180, idx(hit), 'VariableNames', {'h', 'k', 'l', 'variant', ...
    'phi_minus_phi0', 'tth_calc', 'tth_meas', 'd_tth', 'beta_calc', 'beta_meas', 'd_beta', 'spot'});
T = sortrows(T, 'tth_calc');
end

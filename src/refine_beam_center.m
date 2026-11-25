function [beam_center, info] = refine_beam_center(beam_center, geom, peaks_px, T, g, hkl, variant, q, wavelength)
%REFINE_BEAM_CENTER  Direct-beam position fitted to the matched spots of one crystal.
%
%   [beam_center, info] = refine_beam_center(beam_center, geom, peaks_px, T, g, hkl, variant, q, wavelength)
%
%   beam_center : starting direct-beam pixel [row, column]
%   geom        : geometry of the frame (RAPIDutils.getDiffractionGeometry);
%                 the distance (radius of the image plate) is kept
%   peaks_px    : measured spots in pixels, row k = spot k of expt (FIND_BRAGG_PEAKS)
%   T           : matched spots (MATCH_SPOTS), T.spot = row of peaks_px
%   g, hkl, variant : reflections (PHASE_REFLECTIONS)
%   q           : orientation (quaternion) of the crystal
%   wavelength  : (A)
%
%   No calibrant is needed: the matched spots of a crystal whose cell and
%   orientation are known fix the beam centre. The pixel of every measured
%   spot is converted to (2theta, beta) with a trial beam centre and
%   compared with the predicted position of its reflection. Free parameters:
%     beam centre (row, column)                                     2
%     small tilt of the crystal about x and y, as nuisance           2
%   (a rotation about z, the rotation axis of the sample, changes only the
%   rotation angle phi at which a reflection diffracts, not its 2theta or
%   beta, so it cannot be determined from the spot positions and is not fitted)
%   The residuals are d2theta and dbeta*sin(2theta) (both deg, i.e. arc on
%   the Ewald sphere, so that beta of spots near the beam does not dominate).
%   The cell is kept: it is refined in turn by REFINE_CELL (see
%   ORIENTATION_PIPELINE). A shift of the beam centre moves the spots in
%   2theta differently for different beta, a change of the cell moves all
%   in the same direction, so the two can be separated if the matched spots
%   cover a range of beta. The distance cannot be separated from the cell
%   and is not refined.
%
%   info.esd         standard deviation of [row, column] (px)
%   info.rms_before, info.rms_after   rms of the residuals (deg)
%   info.n           number of matched spots used
%   info.rotation    nuisance tilt [about x, about y] (deg) found

n = height(T);
px = peaks_px(T.spot, 1:2);                              % [row, column]
[~, ig] = ismember([T.h T.k T.l T.variant], [hkl variant], 'rows');
g_lab = quatrotate(quatconj(q), g(ig, :));

% which of the two diffraction solutions each spot is
[~, b1, b2] = ewald_full(g_lab, wavelength);
use2 = abs(mod(b2 - T.beta_calc + 180, 360) - 180) < abs(mod(b1 - T.beta_calc + 180, 360) - 180);

resid = @(p) residuals(p, beam_center, geom, px, g_lab, use2, wavelength);
p = zeros(4, 1);                                         % [d_row d_col w_x w_y]
r = resid(p);
rms0 = sqrt(mean(r.^2));
% Gauss-Newton (the problem is nearly linear for small shifts)
step = [0.1; 0.1; 1e-3; 1e-3];
for it = 1:30
    J = jacobian(resid, p, r, step);
    dp = -(J' * J) \ (J' * r);
    for half = 1:10                                      % shorten the step until it helps
        r_new = resid(p + dp);
        if sum(r_new.^2) < sum(r.^2), break; end
        dp = dp / 2;
    end
    if sum(r_new.^2) >= sum(r.^2), break; end
    p_new = p + dp;
    p = p_new;  r = r_new;
    if max(abs(dp(1:2))) < 1e-3, break; end
end
J = jacobian(resid, p, r, step);
s2 = sum(r.^2) / max(numel(r) - numel(p), 1);
C = s2 * inv(J' * J);

beam_center = beam_center + p(1:2)';
info.esd = sqrt(diag(C(1:2, 1:2)))';
info.rms_before = rms0;
info.rms_after = sqrt(mean(r.^2));
info.n = n;
info.rotation = p(3:4)';
end

%% ------------------------------------------------------------------------
function r = residuals(p, bc0, geom, px, g_lab, use2, wavelength)
geom.zero = fliplr(bc0 + p(1:2)');                      % geom.zero = [column, row]
[tth_m, beta_m] = RAPIDutils.ipPixels2diff(geom, px(:,2), px(:,1));
tth_m = tth_m * 180/pi;  beta_m = beta_m * 180/pi;

w = [p(3); p(4); 0] * pi/180;                           % small tilt (Rodrigues)
a = norm(w);
if a > 0
    k = w / a;  Kx = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    R = eye(3) + sin(a) * Kx + (1 - cos(a)) * Kx^2;
else
    R = eye(3);
end
[tth_c, b1, b2] = ewald_full((R * g_lab')', wavelength);
beta_c = b1;  beta_c(use2) = b2(use2);

d_tth = tth_c - tth_m;
d_beta = mod(beta_c - beta_m + 180, 360) - 180;
r = [d_tth; d_beta .* sind(tth_m)];
r(isnan(r)) = 10;                                       % reflection left the Ewald sphere
end

function J = jacobian(f, p, r, step)
J = zeros(numel(r), numel(p));
for k = 1:numel(p)
    dp = zeros(size(p));  dp(k) = step(k);
    J(:, k) = (f(p + dp) - f(p - dp)) / (2 * step(k));
end
end

function [tth, b1, b2] = ewald_full(g, wavelength)
% EWALD_ANGLES with NaN for vectors that do not reach the Ewald sphere
tth = nan(size(g, 1), 1);  b1 = tth;  b2 = tth;
[t, x1, x2, ~, ~, missed] = ewald_angles(g, wavelength);
tth(~missed) = t;  b1(~missed) = x1;  b2(~missed) = x2;
end

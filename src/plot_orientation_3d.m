function fig = plot_orientation_3d(q, phi0, phase, geom, T, S, name, n_rays)
%PLOT_ORIENTATION_3D  What an orientation (q, phi0) means in the experiment, in 3D.
%
%   fig = plot_orientation_3d(q, phi0, phase, geom, T, S, name, n_rays)
%
%   q, phi0  orientation from SEARCH_ORIENTATION (or orientation.mat)
%   phase    struct from PHASE_FROM_CIF (cell used for the drawing and indices)
%   geom     detector geometry (RAPIDutils.getDiffractionGeometry with the beam centre)
%   T        matched spots of this orientation (MATCH_SPOTS)
%   S        settings (S.rotation_range, S.wavelength)
%   name     label for the title
%   n_rays   number of diffracted beams drawn in panel 4 (lowest 2theta first)
%
%   Four panels: 1. the unit cell in the reference position (a || x = beam,
%   b || y, c || z = rotation axis) inside the cylindrical image plate;
%   2. after the quaternion q (a tilt about a horizontal axis that brings the
%   crystal direction n to z); 3. after phi0 (turn about z) = the orientation in
%   the middle of the exposure; 4. the rotation during the exposure with the
%   diffracted beams of matched reflections ending at their spots. Prints the
%   axis and angle of q, the crystal directions along x, y, z and the Euler
%   angles (ORIENTATION_TO_EULER), and checks that the drawn beams hit the
%   plate at the predicted spots.

if nargin < 8, n_rays = 12; end
[~, B, A] = calculate_reciprocal_lattice(struct('a', phase.cell(1:3), 'u', phase.cell(4:6)), [0 0 0], [1 0 0 0]);
rot = S.rotation_range;

Rz = @(a) [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
L0 = quat2dcm(q)';                      % crystal -> lab after the quaternion (phi = 0)
L  = Rz(phi0) * L0;                     % ... and after phi0: the found orientation

% the quaternion as an angle about an axis
q = q / norm(q);
q_angle = 2 * acosd(min(1, abs(q(1))));
if q_angle > 1e-6, q_axis = sign(q(1)) * q(2:4) / norm(q(2:4)); else, q_axis = [1 0 0]; end
n_cart = (L0' * [0; 0; 1])';            % crystal direction that q brings to z
[n_uvw, n_dev] = nearest_indices(n_cart, A, 6);

% crystal directions along the lab axes in the middle of the exposure
[ux, dx] = nearest_indices((L' * [1; 0; 0])', A, 6);
[uy, dy] = nearest_indices((L' * [0; 1; 0])', A, 6);
[uz, dz] = nearest_indices((L' * [0; 0; 1])', A, 6);

txt = {
    sprintf('%s: q = [%.4f %.4f %.4f %.4f], phi0 = %.2f deg, rotation during the exposure %g to %g deg', ...
        name, q, phi0, rot)
    sprintf('Quaternion q = rotation by %.2f deg about the axis [%.3f %.3f %.3f] (horizontal, q(4) = 0)', q_angle, q_axis)
    sprintf('  it brings the crystal direction [%.3f %.3f %.3f] ~ [%s] (%.1f deg off) to the rotation axis z', ...
        n_cart / max(abs(n_cart)), idx2str(n_uvw), n_dev)
    sprintf('phi0 = %.2f deg: then the crystal is turned about z by phi0 - the orientation in the middle of the exposure', phi0)
    'In the middle of the exposure:'
    sprintf('  along the beam (x):          [%s]  (%.1f deg off)', idx2str(ux), dx)
    sprintf('  horizontal, across (y):      [%s]  (%.1f deg off)', idx2str(uy), dy)
    sprintf('  along the rotation axis (z): [%s]  (%.1f deg off)  - the same as n, phi0 does not move it', idx2str(uz), dz)
    };
fprintf('%s\n', txt{:});

% Euler angles (Bunge, as in EBSD) - see src/orientation_to_euler.m for the steps
euler  = orientation_to_euler(q, phi0);
euler0 = orientation_to_euler(q, 0);              % the quaternion alone
txt_euler = {
    sprintf('Euler angles (Bunge z-x''-z'', sample frame X = beam, Z = rotation axis): phi1 = %.2f, Phi = %.2f, phi2 = %.2f deg', euler)
    sprintf('  from q alone (phi0 = 0): phi1 = %.2f, Phi = %.2f, phi2 = %.2f deg; phi0 = %.2f deg just adds to phi1', euler0, phi0)
    sprintf('  Phi = angle between the crystal axis c [001] and the rotation axis Z (= the tilt of q here, %.2f deg)', q_angle)
    };
fprintf('%s\n', txt_euler{:});

% detector
Rp = 1;                                              % plate radius in the drawing
mm = Rp / geom.distance;                             % drawing units per mm
% columns < beam column = towards +y (long side), rows < beam row = up
ang_pos = (geom.zero(1) - 1) * geom.pixSize(1) / geom.distance * 180 / pi;           % towards +y (deg)
ang_neg = (geom.pixNumber(1) - geom.zero(1)) * geom.pixSize(1) / geom.distance * 180 / pi;  % towards -y
h_up    = (geom.zero(2) - 1) * geom.pixSize(2) * mm;
h_down  = (geom.pixNumber(2) - geom.zero(2)) * geom.pixSize(2) * mm;
plate = struct('a1', -ang_neg, 'a2', ang_pos, 'h1', -h_down, 'h2', h_up, 'R', Rp);
fprintf('Image plate: from %.0f deg (towards -y) to %.0f deg (towards +y) around z, height %+.0f to %+.0f mm\n', ...
    -ang_neg, ang_pos, -h_down / mm, h_up / mm);

% diffracted beams of the matched reflections: k_out = k_in + g at the angle phi
K = 2 * pi / S.wavelength;
T = sortrows(T, 'tth_calc');
T = T(1:min(n_rays, height(T)), :);
rays = zeros(height(T), 3);  hits = rays;  phis = zeros(height(T), 1);
for k = 1:height(T)
    phis(k) = phi0 + T.phi_minus_phi0(k);
    g_lab = Rz(phis(k)) * L0 * ([T.h(k) T.k(k) T.l(k)] * B)';
    k_out = [K 0 0] + g_lab';
    rays(k, :) = k_out / norm(k_out);
    hits(k, :) = rays(k, :) * Rp / norm(rays(k, 1:2));    % where it meets the cylinder
end
% check: the drawn hit points against the predicted spots on the detector
[px, pz] = RAPIDutils.diff2ipPixels(geom, T.tth_calc * pi / 180, T.beta_calc * pi / 180);
col3d = geom.zero(1) - atan2d(hits(:,2), hits(:,1)) * pi / 180 * geom.distance / geom.pixSize(1);
row3d = geom.zero(2) - hits(:,3) / mm / geom.pixSize(2);
fprintf('Check of the drawing: diffracted beams hit the plate within %.1f px of the predicted spots\n', ...
    max(hypot(col3d - px, row3d - pz)));

% figure
fig = figure('Color', 'w', 'Position', [40 40 1500 820], 'Name', 'orientation');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
As = A / max(vecnorm(A, 2, 2)) * 0.5;           % cell edge in the drawing

% 1. start
ax = nexttile(tl);  setup_axes(ax, plate, true);
draw_cell(ax, eye(3), As, 1);
title(ax, '1. START: a || x (beam), b || y, c || z (rotation axis)');

% 2. quaternion
ax = nexttile(tl);  setup_axes(ax, plate, false);
draw_cell(ax, eye(3), As, 0.12);                 % ghost of the start
draw_cell(ax, L0, As, 1);
plot3(ax, [-0.75 0.75] * q_axis(1), [-0.75 0.75] * q_axis(2), [0 0], 'm-', 'LineWidth', 2.5);
text(ax, 0.78 * q_axis(1), 0.78 * q_axis(2), 0.04, sprintf('axis of q\n(%.1f%s)', q_angle, char(176)), 'Color', 'm', 'FontWeight', 'bold');
quiver3(ax, 0, 0, 0, 0, 0, 0.8, 0, 'Color', [0.85 0.5 0], 'LineWidth', 3, 'MaxHeadSize', 0.4);
text(ax, 0.03, 0.03, 0.86, sprintf('[%s] -> z', idx2str(n_uvw)), 'Color', [0.85 0.5 0], 'FontWeight', 'bold');
title(ax, sprintf('2. QUATERNION: tilt by %.1f%s about a horizontal axis, [%s] goes to z', q_angle, char(176), idx2str(n_uvw)));

% 3. phi0
ax = nexttile(tl);  setup_axes(ax, plate, false);
draw_cell(ax, L0, As, 0.12);                     % ghost of step 2
draw_cell(ax, L, As, 1);
draw_arc(ax, 0, phi0, 0.55, 0.62, [0 0.45 0.9], sprintf('\\phi_0 = %.1f%s', phi0, char(176)));
title(ax, sprintf('3. PHI0: turn about z by %.1f%s = orientation in the middle of the exposure', phi0, char(176)));

% 4. exposure
ax = nexttile(tl);  setup_axes(ax, plate, false);
draw_cell(ax, Rz(rot(1)) * L, As, 0.12);
draw_cell(ax, Rz(rot(2)) * L, As, 0.12);
draw_cell(ax, L, As, 1);
draw_arc(ax, phi0 + rot(1), phi0 + rot(2), 0.55, 0.62, [0 0.45 0.9], ...
    sprintf('\\phi_0 %+g ... %+g%s', rot(1), rot(2), char(176)));
for k = 1:height(T)
    plot3(ax, [0 hits(k,1)], [0 hits(k,2)], [0 hits(k,3)], '-', 'Color', [1 0.55 0 0.7], 'LineWidth', 1);
    plot3(ax, hits(k,1), hits(k,2), hits(k,3), 'o', 'MarkerFaceColor', [1 0.55 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
    text(ax, hits(k,1), hits(k,2), hits(k,3) + 0.04, sprintf('%d%d%d', abs([T.h(k) T.k(k) T.l(k)])), 'FontSize', 7);
end
title(ax, sprintf('4. EXPOSURE: sample turns %g...%g%s about z; %d strongest matched reflections -> spots', ...
    rot, char(176), height(T)));

title(tl, sprintf(['%s: q = [%.3f %.3f %.3f %.3f], \\phi_0 = %.2f%s   ->   Euler (Bunge) ' ...
    '\\phi_1 = %.2f%s, \\Phi = %.2f%s, \\phi_2 = %.2f%s'], strrep(name, '_', '\_'), q, phi0, char(176), ...
    euler(1), char(176), euler(2), char(176), euler(3), char(176)), 'FontWeight', 'bold');
describe_figure(fig, { ...
    ['HOW THE ORIENTATION IS DESCRIBED: the search does not try all 3D rotations directly. It first chooses which crystal direction n ' ...
     'lies along the rotation axis z (2 angles) - that is the quaternion q, a tilt about a horizontal axis (magenta) - and then how ' ...
     'far the crystal is turned about z (1 angle) - that is phi0. Together 3 angles = one orientation.'], ...
    ['DRAWING: grey cylinder = image plate (sample in its axis; it reaches further towards +y than towards -y, the beam centre ' ...
     'is the red dot where the direct beam hits it, the black x is the middle of the plate, 42.6 deg further towards +y); red arrow = X-ray beam (+x); dashed = rotation axis z. Cell: a = [100] red, ' ...
     'b = [010] green, c = [001] blue; faint cell = previous step. Panel 4: faint cells = start and end of the rotation during the ' ...
     'exposure; orange lines = diffracted beams of matched reflections, each leaving at its own angle phi in that range, ending at ' ...
     'its spot on the plate. The numbers (axis and angle of q, [uvw] along x, y, z) are printed in the command window. ' ...
     'In MATLAB the panels can be turned with the mouse (Rotate 3D in the figure toolbar).'], ...
    sprintf(['EULER ANGLES (Bunge, as in EBSD; sample frame X = beam, Z = rotation axis; src/orientation_to_euler.m): q -> rotation matrix ' ...
     'quat2dcm(q); with phi0: g = (Rz(phi0) * quat2dcm(q)'')'' (sample -> crystal); Phi = acos(g33), phi1 = atan2(g31, -g32), ' ...
     'phi2 = atan2(g13, g23). Here phi1 = %.2f, Phi = %.2f, phi2 = %.2f deg. q alone gives Phi = %.2f (= the tilt of q, the angle ' ...
     'between c [001] and Z) and phi2 = %.2f deg; phi0 = %.2f deg only adds to phi1 (%.2f + %.2f). Crystal symmetry gives equivalent ' ...
     'sets (24 for cubic).'], euler, euler0(2), euler0(3), phi0, euler0(1), phi0)});
end

%% ------------------------------------------------------------------------
function setup_axes(ax, plate, show_labels)
hold(ax, 'on');
% image plate: part of a cylinder around z
a = linspace(plate.a1, plate.a2, 120);
[Aa, Hh] = meshgrid(a, [plate.h1 plate.h2]);
surf(ax, plate.R * cosd(Aa), plate.R * sind(Aa), Hh, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.12, ...
    'EdgeColor', 'none');
for h = [plate.h1 plate.h2]
    plot3(ax, plate.R * cosd(a), plate.R * sind(a), h + 0 * a, '-', 'Color', [0.5 0.5 0.5]);
end
for a_edge = [plate.a1 plate.a2]
    plot3(ax, plate.R * cosd(a_edge) * [1 1], plate.R * sind(a_edge) * [1 1], [plate.h1 plate.h2], '-', 'Color', [0.5 0.5 0.5]);
end
% beam, beam centre, rotation axis
quiver3(ax, -1.3, 0, 0, 1.15, 0, 0, 0, 'Color', 'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.25);
plot3(ax, [0 plate.R], [0 0], [0 0], 'r:', 'LineWidth', 1);
plot3(ax, plate.R, 0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
am = (plate.a1 + plate.a2) / 2;  hm = (plate.h1 + plate.h2) / 2;           % middle of the plate
plot3(ax, plate.R * cosd(am), plate.R * sind(am), hm, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
plot3(ax, [0 0], [0 0], [plate.h1 plate.h2] * 1.05, 'k--', 'LineWidth', 1);
if show_labels
    text(ax, -1.3, 0, 0.08, 'X-ray beam (+x)', 'Color', 'r', 'FontWeight', 'bold');
    text(ax, plate.R, 0, -0.1, 'beam centre', 'Color', 'r');
    text(ax, 0, 0, plate.h2 * 1.1, 'rotation axis z', 'HorizontalAlignment', 'center');
    text(ax, plate.R * cosd(60), plate.R * sind(60), plate.h1 + 0.1, 'image plate', 'Color', [0.4 0.4 0.4]);
    text(ax, 0, 0.45, -0.75, '+y', 'Color', [0.4 0.4 0.4]);
end
axis(ax, 'equal');  axis(ax, 'off');
view(ax, [-38 18]);
xlim(ax, [-1.35 1.05]);  ylim(ax, [-1.05 1.05]);  zlim(ax, [plate.h1 plate.h2] * 1.12);
end

function draw_cell(ax, L, As, alpha)
% unit cell (rows of As = a, b, c in crystal Cartesian) turned by L, centred at 0
[i, j, k] = ndgrid(0:1);
f = [i(:) j(:) k(:)];
V = (f - 0.5) * As * L';
F = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
patch(ax, 'Vertices', V, 'Faces', F, 'FaceColor', [0.75 0.82 1], 'FaceAlpha', 0.35 * alpha, ...
    'EdgeColor', [0.2 0.2 0.4], 'EdgeAlpha', alpha, 'LineWidth', 1.2);
cols = [0.85 0 0; 0 0.6 0; 0 0 0.9];
lab = {'a [100]', 'b [010]', 'c [001]'};
for n = 1:3
    v = As(n, :) * L' / norm(As(n, :)) * 0.55;
    h = quiver3(ax, 0, 0, 0, v(1), v(2), v(3), 0, 'Color', [cols(n, :) alpha], 'LineWidth', 2, 'MaxHeadSize', 0.4);
    if alpha < 1, h.LineStyle = ':'; end
    if alpha == 1
        text(ax, 1.15 * v(1), 1.15 * v(2), 1.15 * v(3), lab{n}, 'Color', cols(n, :), 'FontWeight', 'bold', 'FontSize', 9);
    end
end
end

function draw_arc(ax, a1, a2, radius, height, col, label)
% arrow-like arc about z from angle a1 to a2 (deg)
a = linspace(a1, a2, 60);
plot3(ax, radius * cosd(a), radius * sind(a), height + 0 * a, '-', 'Color', col, 'LineWidth', 2.5);
plot3(ax, radius * cosd(a(end)), radius * sind(a(end)), height, '>', 'Color', col, 'MarkerFaceColor', col);
plot3(ax, [0 radius * cosd(a1)], [0 radius * sind(a1)], [height height], ':', 'Color', col);
plot3(ax, [0 radius * cosd(a2)], [0 radius * sind(a2)], [height height], ':', 'Color', col);
am = (a1 + a2) / 2;
text(ax, 1.1 * radius * cosd(am), 1.1 * radius * sind(am), height + 0.06, label, 'Color', col, 'FontWeight', 'bold');
end

function s = idx2str(h)
s = strjoin(arrayfun(@(x) sprintf('%d', x), round(h), 'UniformOutput', false), ' ');
end

function [uvw, dev] = nearest_indices(v, A, nmax)
% integer direction [uvw] (|u|,|v|,|w| <= nmax) closest to the Cartesian vector v
[h, k, l] = ndgrid(-nmax:nmax);
H = [h(:) k(:) l(:)];
H = H(any(H, 2), :);
H = H(gcd(gcd(abs(H(:,1)), abs(H(:,2))), abs(H(:,3))) == 1, :);
D = H * A;
c = (D ./ vecnorm(D, 2, 2)) * (v(:) / norm(v));
[cmax, i] = max(c);
uvw = H(i, :);
dev = acosd(min(cmax, 1));
end

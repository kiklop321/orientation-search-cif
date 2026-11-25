function [fig, fig_anim] = plot_orientation_steps(q, phi0, phase, name, varargin)
%PLOT_ORIENTATION_STEPS  How the vector n, the quaternion q and phi0 turn the unit cell.
%
%   [fig, fig_anim] = plot_orientation_steps(q, phi0, phase, name)
%   [fig, fig_anim] = plot_orientation_steps(..., 'animate', true, 'gif', gif_file, 'mp4', mp4_file)
%
%   q, phi0  orientation from SEARCH_ORIENTATION (or orientation.mat)
%   phase    struct from PHASE_FROM_CIF
%   name     label for the title
%   options  'animate' true/false (default false): play the movement in a
%            second figure; 'gif', 'mp4': file names to save it (default
%            none); 'verbose' (default true): print n, the tilt and phi0
%
%   The search describes an orientation by a direction n FIXED IN THE
%   CRYSTAL (2 parameters: a unit vector), the quaternion q of the shortest
%   rotation that brings n to the rotation axis z (a tilt about the
%   horizontal axis n x z, q(4) = 0), and a turn phi0 about z (1 parameter).
%   fig shows six stages: start, half tilt, n on z, half phi0, final, final
%   seen along z. fig_anim (with 'animate') plays the same movement.

ip = inputParser;
ip.addParameter('animate', false);
ip.addParameter('gif', '');
ip.addParameter('mp4', '');
ip.addParameter('verbose', true);
ip.parse(varargin{:});
animate = ip.Results.animate;  gif_file = ip.Results.gif;  mp4_file = ip.Results.mp4;  verbose = ip.Results.verbose;
q = q(:)' / norm(q);
[~, ~, A] = calculate_reciprocal_lattice(struct('a', phase.cell(1:3), 'u', phase.cell(4:6)), [0 0 0], [1 0 0 0]);

Rz = @(a) [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
L0 = quat2dcm(q)';                                  % crystal -> lab after q
n_c = (L0' * [0; 0; 1])';                           % n in the crystal frame (Cartesian)
n_uvw = nearest_uvw(n_c, A);

% q as a tilt: angle and horizontal axis, and the tilt by a fraction t of it
tilt = acosd(max(-1, min(1, n_c(3) / norm(n_c))));   % angle between n and z at the start
ax_q = cross(n_c, [0 0 1]);
if norm(ax_q) < 1e-12, ax_q = [1 0 0]; end
ax_q = ax_q / norm(ax_q);
rot_t = @(t) rodrigues(ax_q, t * tilt);
assert(norm(rot_t(1) - L0) < 1e-9, 'q is not the shortest rotation n -> z');

if verbose, fprintf('n = [%s] (in the crystal; %.1f deg from c [001]) is %.2f deg from z at the start\n', ...
    n_uvw, acosd(n_c(3) / norm(n_c)), tilt); end
if verbose, fprintf('q: tilt by %.2f deg about the horizontal axis [%.3f %.3f %.3f]; then phi0 = %.2f deg about z\n', tilt, ax_q, phi0); end

% the movement: tilt 0 -> 1, then phi0 0 -> 1
stage = @(s) stage_matrix(s, rot_t, Rz, phi0);       % s in [0, 2]

% figure with six stages
As = A / max(vecnorm(A, 2, 2)) * 0.9;               % cell edge in the drawing
fig = figure('Color', 'w', 'Position', [40 40 1500 820], 'Name', 'orientation steps');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
stages = [0 0.5 1 1.5 2];
names = {sprintf('1. START: a || x, b || y, c || z; n = [%s] is %.1f%s from z', n_uvw, tilt, char(176)), ...
         sprintf('2. tilt by q, halfway (%.1f%s)', tilt / 2, char(176)), ...
         sprintf('3. after q: n on z (tilt %.1f%s about the magenta axis)', tilt, char(176)), ...
         sprintf('4. turn about z by \\phi_0, halfway (%.1f%s); n stays on z', phi0 / 2, char(176)), ...
         sprintf('5. FINAL: after \\phi_0 = %.1f%s (middle of the exposure)', phi0, char(176))};
for k = 1:numel(stages)
    ax = nexttile(tl);
    draw_stage(ax, stages(k), stage, As, n_c, ax_q, tilt, phi0, n_uvw);
    title(ax, names{k}, 'FontSize', 10);
end
ax = nexttile(tl);                                   % final, seen from above
draw_stage(ax, 2, stage, As, n_c, ax_q, tilt, phi0, n_uvw);
view(ax, [0 90]);
title(ax, sprintf('6. FINAL seen from above (along z): a, b, c turned by \\phi_0 about n'), 'FontSize', 10);

title(tl, sprintf('%s: n = [%s] -> z by q (tilt %.1f%s), then \\phi_0 = %.1f%s about z', ...
    strrep(name, '_', '\_'), n_uvw, tilt, char(176), phi0, char(176)), 'FontWeight', 'bold');
describe_figure(fig, { ...
    sprintf(['n (orange) is a direction of the crystal, here [%s]: it is fixed in the cell and turns with it. The quaternion q tilts ' ...
     'the whole cell by %.1f deg about the horizontal axis n x z (magenta) until n points along z, the rotation axis of the sample ' ...
     '(panels 1-3; the orange dotted arc is the path of the tip of n). This fixes 2 of the 3 angles of the orientation: which crystal ' ...
     'direction is vertical. phi0 then turns the whole cell about z by %.1f deg (panels 3-5; the red dotted arc is the path of a [100]); n ' ...
     'does not move, a, b and c go round it. Together: the orientation in the middle of the exposure; during the exposure the ' ...
     'sample keeps turning about z, phi0 -30 ... +30 deg.'], n_uvw, tilt, phi0), ...
    ['Black arrow = X-ray beam (+x), dashed = rotation axis z, faint cell = start position. Cell axes a [100] red, b [010] green, ' ...
     'c [001] blue. In MATLAB the panels can be turned with the mouse (Rotate 3D).']});

% animation, played in a figure; saved as GIF / MP4 if file names are given
fig_anim = [];
if animate
    if ~isempty(gif_file) && isfile(gif_file), delete(gif_file); end
    f2 = figure('Color', 'w', 'Position', [60 60 700 640], 'Name', 'orientation animation');  fig_anim = f2;
    ax = axes(f2);
    if ~isempty(mp4_file), vid = VideoWriter(mp4_file, 'MPEG-4');  vid.FrameRate = 20;  vid.Quality = 95;  open(vid); end
    s_all = [zeros(1, 8), linspace(0, 1, 40), ones(1, 10), linspace(1, 2, 50), 2 * ones(1, 15)];
    frame_size = [];
    for s = s_all
        cla(ax);
        draw_stage(ax, s, stage, As, n_c, ax_q, tilt, phi0, n_uvw);
        if s <= 1
            title(ax, sprintf('q: tilt %.1f%s of %.1f%s - n = [%s] goes to z', s * tilt, char(176), tilt, char(176), n_uvw));
        else
            title(ax, sprintf('\\phi_0: turn about z %.1f%s of %.1f%s - n stays on z', (s - 1) * phi0, char(176), phi0, char(176)));
        end
        drawnow;
        if ~isempty(gif_file), exportgraphics(f2, gif_file, 'Append', true, 'Resolution', 70); end
        if ~isempty(mp4_file)
            fr = getframe(f2);
            if isempty(frame_size), frame_size = size(fr.cdata, [1 2]); end
            writeVideo(vid, imresize(fr.cdata, frame_size));      % all frames the same size
        end
        if isempty(gif_file) && isempty(mp4_file), pause(0.04); end   % only played: real speed
    end
    if ~isempty(mp4_file), close(vid); end
end

end

%% ------------------------------------------------------------------------
function L = stage_matrix(s, rot_t, Rz, phi0)
% s = 0..1: tilt by the fraction s of q; s = 1..2: then turn about z by (s-1)*phi0
if s <= 1
    L = rot_t(s);
else
    L = Rz((s - 1) * phi0) * rot_t(1);
end
end

function draw_stage(ax, s, stage, As, n_c, ax_q, tilt, phi0, n_uvw)
hold(ax, 'on');
L = stage(s);
% laboratory: beam and rotation axis
quiver3(ax, -1.5, 0, 0, 0.9, 0, 0, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.4);
text(ax, -1.5, 0, 0.12, 'beam +x', 'FontSize', 8);
plot3(ax, [0 0], [0 0], [-1.1 1.25], 'k--');
text(ax, 0.03, 0, 1.3, 'z', 'FontWeight', 'bold');
% start position, faint
draw_cell(ax, eye(3), As, 0.15);
% path of the tip of n during the tilt, and of a during phi0
t = linspace(0, min(s, 1), 40);
tip = zeros(numel(t), 3);
for i = 1:numel(t), tip(i, :) = (stage(t(i)) * n_c')' / norm(n_c) * 1.05; end
plot3(ax, tip(:,1), tip(:,2), tip(:,3), ':', 'Color', [0.95 0.5 0], 'LineWidth', 2);
if s > 1
    p = linspace(1, s, 40);
    a_tip = zeros(numel(p), 3);
    for i = 1:numel(p), a_tip(i, :) = (stage(p(i)) * (As(1, :) / norm(As(1, :)))')' * 0.8; end
    plot3(ax, a_tip(:,1), a_tip(:,2), a_tip(:,3), ':', 'Color', [0.85 0 0], 'LineWidth', 2);
end
% tilt axis of q (horizontal)
if s < 1.001
    plot3(ax, [-0.9 0.9] * ax_q(1), [-0.9 0.9] * ax_q(2), [0 0], 'm-', 'LineWidth', 2);
end
% the cell and n
draw_cell(ax, L, As, 1);
nv = (L * n_c')' / norm(n_c) * 1.05;
quiver3(ax, 0, 0, 0, nv(1), nv(2), nv(3), 0, 'Color', [0.95 0.5 0], 'LineWidth', 3.5, 'MaxHeadSize', 0.3);
text(ax, nv(1) * 1.08, nv(2) * 1.08, nv(3) * 1.08 + 0.05, sprintf('n = [%s]', n_uvw), ...
    'Color', [0.85 0.4 0], 'FontWeight', 'bold');
axis(ax, 'equal');  axis(ax, 'off');
xlim(ax, [-1.55 1.1]);  ylim(ax, [-1.1 1.1]);  zlim(ax, [-1.1 1.4]);
view(ax, [-35 18]);
end

function draw_cell(ax, L, As, alpha)
[i, j, k] = ndgrid(0:1);
f = [i(:) j(:) k(:)];
V = (f - 0.5) * As * L';
F = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
patch(ax, 'Vertices', V, 'Faces', F, 'FaceColor', [0.75 0.82 1], 'FaceAlpha', 0.3 * alpha, ...
    'EdgeColor', [0.2 0.2 0.4], 'EdgeAlpha', alpha, 'LineWidth', 1.2);
if alpha < 1, return; end
cols = [0.85 0 0; 0 0.6 0; 0 0 0.9];
lab = {'a [100]', 'b [010]', 'c [001]'};
for n = 1:3
    v = As(n, :) * L' / norm(As(n, :)) * 0.8;
    quiver3(ax, 0, 0, 0, v(1), v(2), v(3), 0, 'Color', cols(n, :), 'LineWidth', 2, 'MaxHeadSize', 0.3);
    text(ax, 1.12 * v(1), 1.12 * v(2), 1.12 * v(3), lab{n}, 'Color', cols(n, :), 'FontSize', 8, 'FontWeight', 'bold');
end
end

function R = rodrigues(axis, angle)
k = axis(:) / norm(axis);
K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
R = eye(3) + sind(angle) * K + (1 - cosd(angle)) * K^2;
end

function s = nearest_uvw(v, A)
[h, k, l] = ndgrid(-6:6);
H = [h(:) k(:) l(:)];  H = H(any(H, 2), :);
D = H * A;
[~, i] = max((D ./ vecnorm(D, 2, 2)) * (v(:) / norm(v)));
s = strjoin(arrayfun(@(x) sprintf('%d', x), H(i, :), 'UniformOutput', false), ' ');
end

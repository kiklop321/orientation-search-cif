%% VISUALIZE_TWO_GRAINS  Orientation relationship of two grains found in one frame.
%
% Reads orientation.mat written by ORIENTATION_PIPELINE with S.n_grains > 1
% (e.g. example_two_grains.m) and shows, as in EBSD:
%   - the two unit cells in the sample frame, faces coloured by the inverse
%     pole figure (IPF) colour of the beam direction; the best common plane
%     and the best common direction are drawn in both cells (green)
%   - pole figures {100}, {110}, {111} of both grains (projection along the
%     beam, rotation axis vertical); poles of the two grains closer than
%     tol are ringed in green
%   - inverse pole figure of the beam direction (x) and of the rotation
%     axis (z) for both grains (cubic phases)
% and writes to relationship.txt:
%   - misorientation angle and axis (smallest over the crystal symmetry)
%   - nearest coincidence-site-lattice (CSL) relation and the Brandon
%     criterion (cubic phases)
%   - pairs of lattice planes and of lattice directions of the two grains
%     that are parallel within tol (indices up to max_index), and how many
%     such pairs two grains with a RANDOM misorientation have on average:
%     only if clearly more are found is the relation more than chance
%
% Note: diffraction gives the orientations of the two grains, not the
% position of the boundary between them. A common plane found here is a
% candidate for the boundary / habit plane; whether the boundary really lies
% on it has to be checked in the microstructure (e.g. trace in EBSD or
% light microscopy).
%
% Frame: x = beam, z = rotation axis of the sample, both cells shown at the
% centre phi0 of the exposure of grain 1 (the sample frame).

clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));

%% Settings ===============================================================
result_dir = fullfile(here, 'results', 'Ti24Nb0O_two_grains');   % output of example_two_grains.m
pair       = [1 2];         % which two grains
tol        = 3;             % planes / directions closer than this are parallel (deg)
max_index  = 2;             % Miller indices up to this are tested
max_listed = 12;            % at most this many pairs are listed
colours    = [0 0 0; 1 0 1];   % grain 1, grain 2 (as in orientation_overlay.png)
highlight  = [1 0.8 0];       % common plane and direction

%% Orientations ===========================================================
r = load(fullfile(result_dir, 'orientation.mat'), 'grains_saved', 'phase', 'image_file');
G = r.grains_saved(pair);
phase = r.phase;
Sym = phase.parent_rotations;
cubic = size(Sym, 3) == 24;
[~, B, A] = calculate_reciprocal_lattice(struct('a', phase.cell(1:3), 'u', phase.cell(4:6)), [0 0 0], [1 0 0 0]);
[~, frame_name] = fileparts(r.image_file);

% crystal (Cartesian) -> sample frame: g_sample = Rz(phi0) * quat2dcm(q)' * g_crystal
% (as in MATCH_SPOTS / EWALD_ANGLES)
Rz = @(a) [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
L = cell(1, 2);
for k = 1:2
    L{k} = Rz(G(k).phi0) * quat2dcm(G(k).q)';
end

%% Misorientation =========================================================
M = L{1}' * L{2};                          % crystal 2 -> crystal 1 coordinates
[mis_angle, mis_axis] = disorientation(M, Sym);
uvw_axis = mis_axis / A;                   % direction indices in crystal 1
if cubic, uvw_axis = sort(abs(uvw_axis), 'descend'); end
[axis_int, axis_dev] = nearest_indices(uvw_axis * A, A, 6);

txt = {};
txt{end+1} = sprintf('Frame %s, phase %s, a b c = %.4f %.4f %.4f A', frame_name, phase.name, phase.cell(1:3));
txt{end+1} = sprintf('Grain %d: %d spots, grain %d: %d spots', pair(1), G(1).n_spots, pair(2), G(2).n_spots);
txt{end+1} = '';
for k = 1:2
    [bx, dbx] = nearest_indices((L{k}' * [1; 0; 0])', A, 6);
    [bz, dbz] = nearest_indices((L{k}' * [0; 0; 1])', A, 6);
    txt{end+1} = sprintf('Grain %d: beam (x) along [%s] (%.1f deg off), rotation axis (z) along [%s] (%.1f deg off)', ...
        pair(k), idx2str(bx), dbx, idx2str(bz), dbz); %#ok<SAGROW>
end
txt{end+1} = '';
txt{end+1} = sprintf('Misorientation: %.2f deg about [%.3f %.3f %.3f] of grain %d  (~ <%s>, %.1f deg off)', ...
    mis_angle, uvw_axis / max(abs(uvw_axis)), pair(1), idx2str(axis_int), axis_dev);

csl_str = '';
if cubic
    [csl, csl_dev] = nearest_csl(M, Sym);
    brandon = 15 / sqrt(csl.sigma);
    if mis_angle < 15
        csl_str = 'low-angle boundary (Sigma 1)';
    elseif csl_dev <= brandon
        csl_str = sprintf('Sigma %d (%.2f deg <%s>), deviation %.1f deg <= Brandon %.1f deg', ...
            csl.sigma, csl.angle, idx2str(csl.axis), csl_dev, brandon);
    else
        csl_str = sprintf('no CSL relation; nearest Sigma %d (%.2f deg <%s>) is %.1f deg off (Brandon %.1f deg)', ...
            csl.sigma, csl.angle, idx2str(csl.axis), csl_dev, brandon);
    end
    txt{end+1} = ['CSL: ' csl_str];
end

%% Parallel planes and directions =========================================
H = index_list(max_index);
Nc = unit(H * B);                          % plane normals, crystal frame
Dc = unit(H * A);                          % directions, crystal frame
Pp = parallel_pairs(Nc * L{1}', Nc * L{2}', H, tol);
Pd = parallel_pairs(Dc * L{1}', Dc * L{2}', H, tol);

% the same count for random misorientations (reference for chance coincidences)
rng(1);
n_rand = zeros(500, 1);
for t = 1:numel(n_rand)
    qr = randn(1, 4);
    n_rand(t) = size(parallel_pairs(Nc, Nc * quat2dcm(qr / norm(qr))', H, tol), 1);
end

txt{end+1} = '';
txt{end+1} = sprintf(['Planes parallel within %g deg (indices <= %d): %d pairs; two grains with a random ' ...
    'misorientation have %.1f +- %.1f such pairs'], tol, max_index, size(Pp, 1), mean(n_rand), std(n_rand));
txt = [txt, pair_lines(Pp, H, '(%s)', pair, max_listed)];
txt{end+1} = '';
txt{end+1} = sprintf('Directions parallel within %g deg:', tol);
if cubic
    txt{end+1} = '  (cubic: direction [uvw] is the normal of plane (hkl) = (uvw), so the list is the same)';
end
txt = [txt, pair_lines(Pd, H, '[%s]', pair, max_listed)];
txt{end+1} = '';
txt{end+1} = ['Note: a common plane is a candidate for the boundary / habit plane; the boundary ' ...
    'itself is not seen by diffraction.'];

fprintf('%s\n', txt{:});
fid = fopen(fullfile(result_dir, 'relationship.txt'), 'w');
fprintf(fid, '%s\n', txt{:});
fclose(fid);

%% Figure =================================================================
fig = figure('Color', 'w', 'Position', [40 40 1650 880]);
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- unit cells
ax = nexttile(tl, 1, [2 2]);
hold(ax, 'on');
As = A / max(vecnorm(A, 2, 2));            % cell scaled to edge <= 1
centres = [0 -1.4 0; 0 1.4 0];                   % side by side along y
for k = 1:2
    if cubic
        face = ipf_colour((L{k}' * [1; 0; 0])');
    else
        face = colours(k, :);
    end
    draw_cell(ax, L{k}, centres(k, :), As, face, colours(k, :));
    text(ax, centres(k, 1), centres(k, 2), centres(k, 3) + 1.05, sprintf('grain %d (%d spots)', pair(k), G(k).n_spots), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', colours(k, :), 'FontSize', 11);
end
legend_txt = {};
if ~isempty(Pp)
    b = Pp(1, :);
    for k = 1:2
        n = Nc(b(k), :) * L{k}';
        draw_plane(ax, centres(k, :), n, As, L{k}, highlight);
    end
    legend_txt{end+1} = sprintf('common plane (%s)_%d || (%s)_%d, %.1f deg', ...
        idx2str(H(b(1), :)), pair(1), idx2str(b(4) * H(b(2), :)), pair(2), b(3));
end
if ~isempty(Pd)
    b = Pd(1, :);
    for k = 1:2
        dvec = Dc(b(k), :) * L{k}';
        if k == 2, dvec = b(4) * dvec; end
        p = centres(k, :) + [-1; 1] * 0.85 * dvec;
        plot3(ax, p(:, 1), p(:, 2), p(:, 3), '-', 'Color', highlight * 0.75, 'LineWidth', 4);
    end
    legend_txt{end+1} = sprintf('common direction [%s]_%d || [%s]_%d, %.1f deg', ...
        idx2str(H(b(1), :)), pair(1), idx2str(b(4) * H(b(2), :)), pair(2), b(3));
end
quiver3(ax, -2.2, 0, -1.0, 4.4, 0, 0, 0, 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'MaxHeadSize', 0.15);
text(ax, -2.2, 0, -1.15, 'beam (x)', 'Color', [0.85 0.33 0.1], 'FontWeight', 'bold');
quiver3(ax, 0, 0, -1.0, 0, 0, 2.0, 0, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'MaxHeadSize', 0.2);
text(ax, 0, 0, 1.12, 'rotation axis (z)', 'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'center');
axis(ax, 'equal');  axis(ax, 'off');  view(ax, [-62 14]);  axis(ax, 'tight');  camzoom(ax, 1.15);
title(ax, {sprintf('misorientation %.1f%s about <%s>', mis_angle, char(176), idx2str(axis_int)), ...
    strjoin(legend_txt, newline)}, 'FontWeight', 'normal', 'FontSize', 12);
if cubic
    text(ax, 0, 0, -1.35, 'faces: IPF colour of the beam direction; yellow: common plane and direction', ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

% --- pole figures
families = [1 0 0; 1 1 0; 1 1 1];
tiles = [3 4 7];
for f = 1:3
    axp = nexttile(tl, tiles(f));
    pole_figure(axp, families(f, :), B, Sym, L, colours, tol, pair, f == 1);
end

% --- inverse pole figure
axi = nexttile(tl, 8);
if cubic
    ipf_plot(axi, L, colours, pair);
else
    axis(axi, 'off');
    text(axi, 0.5, 0.5, 'IPF only for cubic phases', 'HorizontalAlignment', 'center');
end

t = sprintf('%s: grains %d and %d', frame_name, pair);
if ~isempty(csl_str), t = [t ' - ' csl_str]; end
title(tl, t, 'Interpreter', 'none', 'FontWeight', 'bold');
exportgraphics(fig, fullfile(result_dir, 'grain_relationship.png'), 'Resolution', 200);
fprintf('Saved to %s\n', result_dir);

%% ========================================================================
%  Local functions
%  ========================================================================
function v = unit(v)
v = v ./ vecnorm(v, 2, 2);
end

function s = idx2str(h)
s = strjoin(arrayfun(@(x) sprintf('%d', x), round(h), 'UniformOutput', false), ' ');
end

function H = index_list(n)
% all (h k l) with |indices| <= n, coprime, first non-zero index positive
[h, k, l] = ndgrid(-n:n);
H = [h(:) k(:) l(:)];
H = H(any(H, 2), :);
H = H(gcd(gcd(abs(H(:,1)), abs(H(:,2))), abs(H(:,3))) == 1, :);
first = zeros(size(H, 1), 1);
for i = 1:size(H, 1), first(i) = H(i, find(H(i, :), 1)); end
H = H(first > 0, :);
end

function [uvw, dev] = nearest_indices(v, A, nmax)
% integer direction [uvw] (|u|,|v|,|w| <= nmax) closest to the Cartesian vector v
H = index_list(nmax);  H = [H; -H];
c = unit(H * A) * (v(:) / norm(v));
[cmax, i] = max(c);
uvw = H(i, :);
dev = acosd(min(cmax, 1));
end

function R = rodrigues(axis, angle)
k = axis(:) / norm(axis);
K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
R = eye(3) + sind(angle) * K + (1 - cosd(angle)) * K^2;
end

function [ang, ax] = disorientation(M, Sym)
% smallest rotation angle of S_i * M * S_j and its axis
ang = inf;  ax = [0 0 1];
for i = 1:size(Sym, 3)
    for j = 1:size(Sym, 3)
        Mij = Sym(:,:,i) * M * Sym(:,:,j);
        a = acosd(max(-1, min(1, (trace(Mij) - 1) / 2)));
        if a < ang - 1e-9
            ang = a;
            [V, D] = eig(Mij);
            [~, m] = min(abs(diag(D) - 1));
            ax = real(V(:, m))';
            s = [Mij(3,2) - Mij(2,3), Mij(1,3) - Mij(3,1), Mij(2,1) - Mij(1,2)];
            if dot(ax, s) < 0, ax = -ax; end
        end
    end
end
end

function [best, dev] = nearest_csl(M, Sym)
% nearest CSL misorientation of a cubic crystal (Sigma <= 29)
T = [ 3 60.00 1 1 1;  5 36.87 1 0 0;  7 38.21 1 1 1;  9 38.94 1 1 0; 11 50.48 1 1 0;
     13 22.62 1 0 0; 13 27.80 1 1 1; 15 48.19 2 1 0; 17 28.07 1 0 0; 17 61.93 2 2 1;
     19 26.53 1 1 0; 19 46.83 1 1 1; 21 21.79 1 1 1; 21 44.42 2 1 1; 23 40.46 3 1 1;
     25 16.26 1 0 0; 25 51.68 3 3 1; 27 31.59 1 1 0; 27 35.43 2 1 0; 29 43.60 1 0 0;
     29 46.40 2 2 1];
dev = inf;
for t = 1:size(T, 1)
    Rc = rodrigues(T(t, 3:5), T(t, 2));
    % angle between Rc and every equivalent S_i * M * S_j; the symmetry must be
    % applied on both sides explicitly (a trace only sees the product S_j * S_i)
    d = inf;
    for i = 1:size(Sym, 3)
        d = min(d, disorientation(Rc' * Sym(:,:,i) * M, Sym));
    end
    if d < dev
        dev = d;
        best = struct('sigma', T(t, 1), 'angle', T(t, 2), 'axis', T(t, 3:5));
    end
end
end

function P = parallel_pairs(V1, V2, H, tol)
% rows [i j deviation sign] of vectors V1(i,:) || V2(j,:) within tol,
% simplest indices first
C = V1 * V2';
[i, j] = find(acosd(min(abs(C), 1)) < tol);
if isempty(i), P = zeros(0, 4); return; end
c = C(sub2ind(size(C), i, j));
P = [i, j, acosd(min(abs(c), 1)), sign(c)];
complexity = sum(H(i, :).^2, 2) + sum(H(j, :).^2, 2);
[~, o] = sortrows([complexity, P(:, 3)]);
P = P(o, :);
end

function lines = pair_lines(P, H, fmt, pair, max_listed)
lines = {};
if isempty(P)
    lines{1} = '  none';
    return
end
for r = 1:min(size(P, 1), max_listed)
    lines{end+1} = sprintf(['  ' fmt '_%d || ' fmt '_%d   %.2f deg'], ...
        idx2str(H(P(r, 1), :)), pair(1), idx2str(P(r, 4) * H(P(r, 2), :)), pair(2), P(r, 3)); %#ok<AGROW>
end
if size(P, 1) > max_listed
    lines{end+1} = sprintf('  ... %d more', size(P, 1) - max_listed);
end
end

function draw_cell(ax, L, c, As, face, edge)
% unit cell (rows of As = cell vectors) rotated by L, centred at c
[i, j, k] = ndgrid(0:1);
f = [i(:) j(:) k(:)];
V = (f - 0.5) * As * L' + c;
F = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
patch(ax, 'Vertices', V, 'Faces', F, 'FaceColor', face, 'FaceAlpha', 0.45, ...
    'EdgeColor', edge, 'LineWidth', 1.8);
lab = {'[100]', '[010]', '[001]'};
for n = 1:3
    v = As(n, :) * L' * 0.8;
    quiver3(ax, c(1), c(2), c(3), v(1), v(2), v(3), 0, 'Color', edge, 'LineWidth', 1.2, 'MaxHeadSize', 0.3);
    text(ax, c(1) + 1.15 * v(1), c(2) + 1.15 * v(2), c(3) + 1.15 * v(3), lab{n}, 'Color', edge, 'FontSize', 9);
end
end

function draw_plane(ax, c, n, As, L, col)
% section of the unit cell through its centre with the plane of normal n
[i, j, k] = ndgrid(0:1);
f = [i(:) j(:) k(:)];
V = (f - 0.5) * As * L';
p = zeros(0, 3);
for a = 1:8
    for b = a+1:8
        if sum(abs(f(a, :) - f(b, :))) ~= 1, continue; end      % cell edge
        da = dot(n, V(a, :));  db = dot(n, V(b, :));
        if da * db <= 0 && da ~= db
            t = da / (da - db);
            p(end+1, :) = V(a, :) + t * (V(b, :) - V(a, :)); %#ok<AGROW>
        end
    end
end
p = uniquetol(p, 1e-9, 'ByRows', true);
if size(p, 1) < 3, return; end
e1 = unit(p(1, :));  e2 = cross(n, e1);
[~, o] = sort(atan2(p * e2', p * e1'));
p = p(o, :) + c;
patch(ax, p(:, 1), p(:, 2), p(:, 3), col, 'FaceAlpha', 0.85, 'EdgeColor', col * 0.6, 'LineWidth', 1.5);
end

function pole_figure(ax, hkl, B, Sym, L, colours, tol, pair, show_legend)
% stereographic pole figure of the family {hkl}, projection along the beam (x)
n0 = unit(hkl * B);
P = zeros(size(Sym, 3), 3);
for s = 1:size(Sym, 3), P(s, :) = (Sym(:,:,s) * n0')'; end
P = uniquetol([P; -P], 1e-6, 'ByRows', true);
hold(ax, 'on');
t = linspace(0, 2*pi, 361);
plot(ax, cos(t), sin(t), 'k-', 'LineWidth', 0.8);
plot(ax, [-1 1], [0 0], ':', 'Color', [0.6 0.6 0.6]);
plot(ax, [0 0], [-1 1], ':', 'Color', [0.6 0.6 0.6]);
mk = {'o', 's'};
Q = cell(1, 2);
h = gobjects(1, 2);
for k = 1:2
    Q{k} = P * L{k}';
    Q{k} = Q{k}(Q{k}(:, 1) >= -1e-9, :);                   % hemisphere towards the beam
    xy = [Q{k}(:, 2), Q{k}(:, 3)] ./ (1 + Q{k}(:, 1));
    h(k) = plot(ax, xy(:, 1), xy(:, 2), mk{k}, 'Color', colours(k, :), 'MarkerSize', 7, 'LineWidth', 1.5);
end
% poles of the two grains closer than tol
ang = acosd(min(abs(Q{1} * Q{2}'), 1));
[i, ~] = find(ang < tol);
if ~isempty(i)
    xy = [Q{1}(i, 2), Q{1}(i, 3)] ./ (1 + Q{1}(i, 1));
    hr = plot(ax, xy(:, 1), xy(:, 2), 'o', 'Color', [0 0.7 0], 'MarkerSize', 16, 'LineWidth', 2);
else
    hr = gobjects(0);
end
axis(ax, 'equal');  axis(ax, [-1.1 1.1 -1.1 1.1]);  axis(ax, 'off');
text(ax, 1.08, 0, 'y', 'FontSize', 9);  text(ax, 0, 1.1, 'z (rot.)', 'FontSize', 9, 'HorizontalAlignment', 'center');
title(ax, sprintf('{%s} pole figure', sprintf('%d', hkl)), 'FontWeight', 'normal', 'Interpreter', 'none');
if show_legend
    lg = {sprintf('grain %d', pair(1)), sprintf('grain %d', pair(2))};
    if ~isempty(hr), h = [h, hr(1)]; lg{end+1} = sprintf('common (< %g%s)', tol, char(176)); end
    legend(ax, h, lg, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
end
end

function v = ipf_reduce(v)
% cubic standard triangle 001-101-111: z >= x >= y >= 0
v = sort(abs(v(:)'), 'ascend');            % [small middle large]
v = [v(2) v(1) v(3)] / norm(v);
end

function rgb = ipf_colour(v)
% IPF colour: red 001, green 101, blue 111
v = ipf_reduce(v);
C = [0 0 1; 1/sqrt(2) 0 1/sqrt(2); 1/sqrt(3) 1/sqrt(3) 1/sqrt(3)]';
w = max(C \ v', 0)';
rgb = w / max(w);
end

function ipf_plot(ax, L, colours, pair)
% inverse pole figure (cubic) of the beam (x) and the rotation axis (z)
[X, Y] = meshgrid(linspace(0, 0.43, 300), linspace(0, 0.38, 270));
d = 1 + X.^2 + Y.^2;
vx = 2 * X ./ d;  vy = 2 * Y ./ d;  vz = (1 - X.^2 - Y.^2) ./ d;
inside = vx >= vy & vy >= 0 & vz >= vx;
img = ones(numel(X), 3);
for i = find(inside)'
    img(i, :) = ipf_colour([vx(i) vy(i) vz(i)]);
end
img = reshape(img, [size(X), 3]);
image(ax, X(1, :), Y(:, 1), img, 'AlphaData', double(inside));
set(ax, 'YDir', 'normal');  hold(ax, 'on');
proj = @(v) [v(:, 1) ./ (1 + v(:, 3)), v(:, 2) ./ (1 + v(:, 3))];
s = linspace(0, 1, 50)';
edge = [proj(unit([s, 0 * s, ones(50, 1)])); ...
        proj(unit([ones(50, 1), s, ones(50, 1)])); ...
        proj(unit([1 - s, 1 - s, ones(50, 1)]))];
plot(ax, edge(:, 1), edge(:, 2), 'k-', 'LineWidth', 1);
text(ax, 0, -0.025, '001', 'HorizontalAlignment', 'center');
text(ax, 0.414, -0.025, '101', 'HorizontalAlignment', 'center');
text(ax, 0.366, 0.39, '111', 'HorizontalAlignment', 'center');
dirs = {[1 0 0], 'x (beam)', 'o'; [0 0 1], 'z (rot.)', '^'};
for k = 1:2
    for m = 1:2
        v = ipf_reduce((L{k}' * dirs{m, 1}')');
        p = proj(v);
        plot(ax, p(1), p(2), dirs{m, 3}, 'MarkerSize', 10, 'LineWidth', 2, ...
            'MarkerEdgeColor', colours(k, :), 'MarkerFaceColor', ipf_colour(v));
        text(ax, p(1) + 0.012, p(2) + 0.012, sprintf('%d: %s', pair(k), dirs{m, 2}), ...
            'FontSize', 8, 'Color', colours(k, :));
    end
end
axis(ax, 'equal');  axis(ax, [-0.03 0.45 -0.05 0.41]);  axis(ax, 'off');
title(ax, 'IPF of beam (o) and rotation axis (\Delta)', 'FontWeight', 'normal');
end

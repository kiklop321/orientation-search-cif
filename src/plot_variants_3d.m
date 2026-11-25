function [fig1, fig2, fig3] = plot_variants_3d(V, parent_cif, child_cif, OR, parent_name, child_name)
%PLOT_VARIANTS_3D  The child variants in the parent crystal, in 3D, with the parallel planes and directions.
%
%   [fig1, fig2, fig3] = plot_variants_3d(V, parent_cif, child_cif, OR, parent_name, child_name)
%
%   V          : 3 x 3 x N variants (child Cartesian -> parent Cartesian),
%                from VARIANTS_MTEX / VARIANTS_FROM_OR
%   parent_cif, child_cif : the CIF files (or structs from READ_CIF)
%   OR         : struct with parent_plane, child_plane, parent_dir, child_dir
%                (3 or 4 indices), or with parent_dirs, child_dirs (two
%                direction pairs); [] or an MTEX orientation = no OR drawn
%   names      : labels
%
%   fig1  overview: left, the parent unit cell (grey) with the c axes of all
%         variants as coloured arrows from its centre; right, variant 1 to
%         scale inside the parent cell with the parallel planes and directions
%   fig2  every variant: its unit cell (to scale) inside the parent cell, and
%         the pair of parallel PLANES (parent: grey section of the parent
%         cell, child: coloured section of the child cell) and the pair of
%         parallel DIRECTIONS (parent: black line, child: coloured dashed
%         line). Each variant realises the OR with its own members of the
%         families, e.g. (1-10)beta || (001)alpha''; their indices are in the
%         title (hexagonal phases in 4-index notation).
%   fig3  viewer: one variant, large; choose it in the menu, turn it with
%         the mouse (rotation is on)
%
%   Hexagonal phases (e.g. omega, alpha') are drawn as the hexagonal prism
%   (three unit cells around c), the others as the unit cell.
%
%   Frames: crystal X || a, Z || c*; the picture frame is the parent crystal.

pc = parent_cif;  cc = child_cif;
if ~isstruct(pc), pc = read_cif(pc); end
if ~isstruct(cc), cc = read_cif(cc); end
[~, Bp, Ap] = calculate_reciprocal_lattice(struct('a', pc.cell(1:3), 'u', pc.cell(4:6)), [0 0 0], [1 0 0 0]);
[~, Bc, Ac] = calculate_reciprocal_lattice(struct('a', cc.cell(1:3), 'u', cc.cell(4:6)), [0 0 0], [1 0 0 0]);
P = struct('A', Ap, 'B', Bp, 'R', point_group_rotations(pc), 'hex', is_hex(pc));
C = struct('A', Ac, 'B', Bc, 'R', point_group_rotations(cc), 'hex', is_hex(cc));
N = size(V, 3);
s = 1 / max(vecnorm(Ap, 2, 2));                       % drawing scale: parent edge = 1
Aps = Ap * s;  Acs = Ac * s;
cols = lines(max(N, 7));
% half-size of the plots: the larger of the two drawn cells (hexagonal: the prism)
lim = 1.05 * max([max(vecnorm(body(Aps, P.hex), 2, 2)), max(vecnorm(body(Acs, C.hex), 2, 2)), 0.6]);

% the OR as two pairs: {type, parent indices, child indices}
pairs = {};
if isstruct(OR) && isfield(OR, 'parent_plane')
    pairs = {'plane', OR.parent_plane, OR.child_plane; 'dir', OR.parent_dir, OR.child_dir};
elseif isstruct(OR) && isfield(OR, 'parent_dirs')
    pairs = {'dir', OR.parent_dirs(1,:), OR.child_dirs(1,:); 'dir', OR.parent_dirs(2,:), OR.child_dirs(2,:)};
end

%% figure 1: overview and variant 1 =======================================
fig1 = figure('Color', 'w', 'Position', [40 40 1400 640], 'Name', 'variants in the parent crystal');
tl = tiledlayout(fig1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);  hold(ax, 'on');
draw_cell(ax, eye(3), Aps, [0.6 0.6 0.6], 0.08, true, {'a', 'b', 'c'}, P.hex);
for j = 1:N
    c = (V(:,:,j) * Acs(3, :)')' / norm(Acs(3, :)) * 0.85;
    quiver3(ax, 0, 0, 0, c(1), c(2), c(3), 0, 'Color', cols(j, :), 'LineWidth', 2.5, 'MaxHeadSize', 0.25);
    text(ax, 1.08 * c(1), 1.08 * c(2), 1.08 * c(3), sprintf('%d: [%s]', j, idx_str(c / Ap, P.hex, 'dir')), ...
        'Color', cols(j, :) * 0.8, 'FontWeight', 'bold', 'FontSize', 8);
end
finish(ax, 0.75);
title(ax, sprintf('c axes of the %d %s variants in the %s cell', N, child_name, parent_name), 'Interpreter', 'none');

ax = nexttile(tl);  hold(ax, 'on');
lab = draw_variant(ax, V(:,:,1), Aps, Acs, cols(1, :), pairs, P, C, true, parent_name, child_name);
finish(ax, lim);
title(ax, [{sprintf('variant 1 in the %s cell', parent_name)}, lab], 'Interpreter', 'none', 'FontSize', 10);
describe_figure(fig1, {['Left: the ' parent_name ' unit cell (grey, axes a b c) and, from its centre, the c axis of every ' ...
    child_name ' variant (colour, with the ' parent_name ' direction it lies along). Right: variant 1 to scale inside the ' ...
    parent_name ' cell (axes of the variant: a'' b'' c''). PARALLEL PLANES: ' parent_name ' = grey section through the ' ...
    'centre of its cell, ' child_name ' = coloured section of the variant cell; PARALLEL DIRECTIONS: ' parent_name ...
    ' = black line, ' child_name ' = coloured dashed line. They coincide.']});

%% figure 2: every variant ================================================
ncol = ceil(sqrt(N * 1.6));  nr = ceil(N / ncol);
fig2 = figure('Color', 'w', 'Position', [60 40 1400 min(330 * nr + 60, 1000)], 'Name', 'all variants');
tl = tiledlayout(fig2, nr, ncol, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:N
    ax = nexttile(tl);  hold(ax, 'on');
    lab = draw_variant(ax, V(:,:,j), Aps, Acs, cols(j, :), pairs, P, C, false, parent_name, child_name);
    finish(ax, lim);
    title(ax, [{sprintf('variant %d', j)}, lab], 'FontSize', 8, 'Interpreter', 'none');
end
describe_figure(fig2, {['Every ' child_name ' variant (coloured cell - hexagonal phases as the hexagonal prism of three cells; axes a red, b green, c blue) to scale inside the ' ...
    parent_name ' unit cell (grey), all in the same view of the ' parent_name ' crystal. Title: the planes and directions that ' ...
    'are parallel in this variant - each variant uses other members of the same families. Planes: ' parent_name ' grey, ' ...
    child_name ' coloured; directions: ' parent_name ' black line, ' child_name ' coloured dashed. ' ...
    'In MATLAB the panels can be turned with the mouse (Rotate 3D).']});

%% figure 3: one variant, large and interactive ===========================
fig3 = figure('Color', 'w', 'Position', [80 40 1000 900], 'Name', 'variant viewer (turn with the mouse)');
ax3 = axes(fig3, 'Position', [0.05 0.12 0.9 0.78]);
uicontrol(fig3, 'Style', 'text', 'Units', 'normalized', 'Position', [0.02 0.95 0.08 0.03], ...
    'String', 'variant:', 'BackgroundColor', 'w', 'FontSize', 11, 'HorizontalAlignment', 'right');
uicontrol(fig3, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [0.11 0.95 0.12 0.035], ...
    'String', compose('%d', 1:N), 'FontSize', 11, 'Callback', @(src, ~) show(src.Value));
uicontrol(fig3, 'Style', 'text', 'Units', 'normalized', 'Position', [0.02 0.01 0.96 0.09], 'BackgroundColor', 'w', ...
    'FontSize', 9, 'HorizontalAlignment', 'left', 'String', ['Turn the picture with the mouse (rotation is on); ' ...
    'choose the variant above. Grey cell: ' parent_name ' (axes a b c); coloured cell: the ' child_name ' variant (axes a'' b'' c''). ' ...
    'Parallel PLANES: ' parent_name ' grey section, ' child_name ' coloured section. Parallel DIRECTIONS: ' parent_name ...
    ' black line, ' child_name ' coloured dashed line. The indices are in the title.']);
show(1);
rotate3d(fig3, 'on');

    function show(j)
        % draw variant j large in the viewer (keeps the current view direction)
        if isempty(ax3.Children), vw = [-35 20]; else, [az, el] = view(ax3); vw = [az el]; end
        cla(ax3);  hold(ax3, 'on');
        lab3 = draw_variant(ax3, V(:,:,j), Aps, Acs, cols(j, :), pairs, P, C, true, parent_name, child_name);
        finish(ax3, lim);  view(ax3, vw);
        title(ax3, [{sprintf('variant %d of %d: %s in %s', j, N, child_name, parent_name)}, lab3], ...
            'Interpreter', 'none', 'FontSize', 12);
    end
end

%% ------------------------------------------------------------------------
function lab = draw_variant(ax, Vj, Aps, Acs, col, pairs, P, C, big, pn, cn)
% parent cell, variant cell and the realised parallel planes / directions
draw_cell(ax, eye(3), Aps, [0.55 0.55 0.55], 0.05, big, {'a', 'b', 'c'}, P.hex);
if big, labels = {'a''', 'b''', 'c'''}; else, labels = {'a', 'b', 'c'}; end
draw_cell(ax, Vj, Acs, col, 0.15, true, labels, C.hex);
lab = {};
for k = 1:size(pairs, 1)
    [vp, vc, ip, ic] = realised_pair(Vj, pairs{k, 1}, pairs{k, 2}, pairs{k, 3}, P, C);
    if strcmp(pairs{k, 1}, 'plane')
        draw_plane(ax, vp, Aps, [0.25 0.25 0.25], 0.3, P.hex);
        draw_plane(ax, vc, Acs * Vj', col, 0.35, C.hex);
        lab{end+1} = sprintf('(%s) %s || (%s) %s   planes', ip, pn, ic, cn); %#ok<AGROW>
    else
        d1 = vp / norm(vp) * 0.75;  d2 = vc / norm(vc) * 0.75;
        plot3(ax, [-1 1] * d1(1), [-1 1] * d1(2), [-1 1] * d1(3), 'k-', 'LineWidth', 4);
        plot3(ax, [-1 1] * d2(1), [-1 1] * d2(2), [-1 1] * d2(3), '--', 'Color', col * 0.8, 'LineWidth', 2.5);
        lab{end+1} = sprintf('[%s] %s || [%s] %s   directions', ip, pn, ic, cn); %#ok<AGROW>
    end
end
end

function [vp, vc, ip, ic] = realised_pair(Vj, type, hp, hc, P, C)
% the members of the parent and child families that are parallel in this variant
% (vp, vc: Cartesian in the parent frame; ip, ic: their indices as text)
if strcmp(type, 'plane')
    np = three(hp, 'hkl') * P.B;  nc = three(hc, 'hkl') * C.B;
else
    np = three(hp, 'uvw') * P.A;  nc = three(hc, 'uvw') * C.A;
end
best = inf;
for m = 1:size(C.R, 3)
    w = (Vj * C.R(:,:,m) * nc')';                     % an equivalent child vector, parent frame
    for k = 1:size(P.R, 3)
        u = (P.R(:,:,k) * np')';
        e = acosd(min(1, max(-1, dot(u, w) / norm(u) / norm(w))));
        if e < best, best = e;  vp = u;  vc = w;  cm = (C.R(:,:,m) * nc')'; end
    end
end
if strcmp(type, 'plane')
    ip = idx_str(vp / P.B, P.hex, 'plane');  ic = idx_str(cm / C.B, C.hex, 'plane');
else
    ip = idx_str(vp / P.A, P.hex, 'dir');    ic = idx_str(cm / C.A, C.hex, 'dir');
end
end

function draw_cell(ax, V, A, col, alpha, axes_on, labels, hex)
% unit cell with basis rows A, turned by V, centred at 0; hexagonal phases as
% the hexagonal prism (three cells around c)
if nargin < 8, hex = false; end
[Pt, F] = body(A, hex);
Pt = Pt * V';
patch(ax, 'Vertices', Pt, 'Faces', F, 'FaceColor', col, 'FaceAlpha', alpha, 'EdgeColor', col * 0.7, 'LineWidth', 1.2);
if ~axes_on, return; end
ac = [0.85 0 0; 0 0.6 0; 0 0 0.9];
for n = 1:3
    v = (V * A(n, :)')' * 0.5;
    quiver3(ax, 0, 0, 0, v(1), v(2), v(3), 0, 'Color', ac(n, :), 'LineWidth', 1.5, 'MaxHeadSize', 0.4);
    if ~isempty(labels)
        text(ax, 1.15 * v(1), 1.15 * v(2), 1.15 * v(3), labels{n}, 'Color', ac(n, :), 'FontWeight', 'bold', 'FontSize', 8);
    end
end
end

function [Pt, F, E] = body(A, hex)
% vertices Pt, faces F (NaN-padded) and edges E of the drawn cell, centred at 0
% A: basis rows a, b, c (already turned, if needed)
if hex
    a1 = A(1, :);  a2 = A(2, :);  c = A(3, :);
    h = [a1; a1 + a2; a2; -a1; -a1 - a2; -a2];       % regular hexagon of edge a
    Pt = [h - c / 2; h + c / 2];
    side = [(1:6)', [2:6 1]', [8:12 7]', (7:12)'];
    F = [[1:6; 7:12], nan(2, 0)];
    F = [F; side, nan(6, 2)];
    E = [(1:6)', [2:6 1]'; (7:12)', [8:12 7]'; (1:6)', (7:12)'];
else
    [i, j, k] = ndgrid(0:1);
    f = [i(:) j(:) k(:)];
    Pt = (f - 0.5) * A;
    F = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
    E = zeros(0, 2);
    for a = 1:8
        for b = a+1:8
            if sum(abs(f(a, :) - f(b, :))) == 1, E(end+1, :) = [a b]; end %#ok<AGROW>
        end
    end
end
end

function draw_plane(ax, n, A, col, alpha, hex)
% section through the centre of the drawn cell (basis rows A, already turned) with the plane of normal n
if nargin < 6, hex = false; end
n = n / norm(n);
[Pt, ~, E] = body(A, hex);
p = zeros(0, 3);
for e = 1:size(E, 1)
    a = E(e, 1);  b = E(e, 2);
    da = dot(n, Pt(a, :));  db = dot(n, Pt(b, :));
    if da * db <= 0 && da ~= db
        t = da / (da - db);
        p(end+1, :) = Pt(a, :) + t * (Pt(b, :) - Pt(a, :)); %#ok<AGROW>
    end
end
p = uniquetol(p, 1e-9, 'ByRows', true);
if size(p, 1) < 3, return; end
e1 = p(1, :) / norm(p(1, :));  e2 = cross(n, e1);
[~, o] = sort(atan2(p * e2', p * e1'));
patch(ax, p(o, 1), p(o, 2), p(o, 3), col, 'FaceAlpha', alpha, 'EdgeColor', col, 'LineWidth', 1.5);
end

function finish(ax, lim)
axis(ax, 'equal');  axis(ax, 'off');
xlim(ax, [-1 1] * lim);  ylim(ax, [-1 1] * lim);  zlim(ax, [-1 1] * lim);
view(ax, [-35 20]);
end

function h = is_hex(cif)
h = abs(cif.cell(6) - 120) < 0.5 && abs(cif.cell(1) - cif.cell(2)) < 1e-3 * cif.cell(1);
end

function v = three(v, type)
% 4-index (hkil / UVTW) -> 3-index (hkl / uvw)
if numel(v) == 4
    if strcmp(type, 'hkl'), v = v([1 2 4]); else, v = [v(1) - v(3), v(2) - v(3), v(4)]; end
end
end

function s = idx_str(v, hex, type)
% integer indices of v (3-index, may be scaled), as text; 4 indices for hexagonal phases
v = v / max(abs(v));
for m = 1:12                                        % smallest multiple giving integers
    w = v * m;
    if all(abs(w - round(w)) < 1e-3), break; end
end
w = round(w);
if hex
    if strcmp(type, 'plane')
        w = [w(1), w(2), -(w(1) + w(2)), w(3)];
    else
        w = [2*w(1) - w(2), 2*w(2) - w(1), -(w(1) + w(2)), 3*w(3)];
    end
end
g = 0;  for x = abs(w), g = gcd(g, x); end
w = w / max(g, 1);
s = strjoin(arrayfun(@(x) sprintf('%d', x), w, 'UniformOutput', false), ' ');
end

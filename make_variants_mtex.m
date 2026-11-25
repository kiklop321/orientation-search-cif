%% MAKE_VARIANTS_MTEX  Orientations of all child variants in a parent crystal, with MTEX.
%
% Generic: give a PARENT phase and a CHILD phase (two CIF files) and the
% orientation relationship (OR) between them. The script computes with
% MTEX (https://mtex-toolbox.github.io) all variants of the child that the
% OR allows in one parent crystal and gives their orientations:
%
%   - in the frame of the parent crystal (always): Euler angles, which
%     parent directions the child axes lie along, misorientations between
%     the variants;
%   - in the sample frame (optional, if the orientation of the parent grain
%     is given - e.g. found by the orientation search, or as Euler angles):
%     Euler angles of every variant as it lies in the sample.
%
% Output (folder variants/, name set below):
%   <name>.mat               V (3x3xN, child -> parent frame) and the tables;
%                            read by the orientation search:
%                            phase = phase_from_cif(child_cif, 'variants', '<name>.mat')
%   <name>_variants.csv      the table of the variants (e.g. for Excel)
%   <name>_polefigures.png   pole figures of the variants in the parent frame (MTEX)
%   <name>_3d_overview.png, <name>_3d_variants.png, <name>_3d_variant1.png
%                            the variants drawn in 3D in the parent cell; the viewer
%                            window (one variant large, choose it, turn it with the mouse)
%                            stays open
%
% Examples of ORs (set in the settings):
%   beta -> alpha'' (orthorhombic):  {110}_b || (001)_a'',   <111>_b || <110>_a''
%   beta -> alpha'  (hexagonal):     {110}_b || (0001)_a',   <111>_b || <2-1-10>_a'   (Burgers)
%   fcc -> bcc:  orientation.KurdjumovSachs(csParent, csChild), orientation.NishiyamaWassermann(...)
% The plane and the direction of each phase must be perpendicular (the
% direction lies in the plane), e.g. [1-11] lies in (110).
%
% MTEX gives the symmetry, the cell and the OR; atoms and structure factors
% of the child (for the search) come from its CIF (phase_from_cif).
% Frames: crystal X || a, Z || c* (as everywhere in this tool); sample frame
% of the search: X = beam, Z = rotation axis.

clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here) || ~isfile(fullfile(here, 'make_variants_mtex.m')), here = fileparts(which('make_variants_mtex.m')); end
addpath(fullfile(here, 'src'));

%% Settings ===============================================================
mtex_path  = '';           % folder with startup_mtex.m, e.g. 'C:\tools\mtex-6.1.1';
                           % '' = MTEX already on the path or started (startup_mtex)

% --- the two phases
parent_cif = fullfile(here, 'cif', 'beta_Ti_COD9012924.cif');           % parent: beta-Ti (Im-3m)
child_cif  = fullfile(here, 'cif', 'omega_Ti_COD9012925.cif');   % child:  omega-Ti (P6/mmm)

% --- the orientation relationship: (plane) || (plane) and [direction] || [direction]
%     3 indices, or 4 for hexagonal phases (hkil / UVTW)
OR.parent_plane = [1 1 1];    OR.child_plane = [0 0 0 1];        % (111)_beta  || (0001)_omega
OR.parent_dir   = [1 -1 0];   OR.child_dir   = [1 1 -2 0];       % [1-10]_beta || [11-20]_omega   (4 variants)
% orthorhombic martensite alpha'' instead (OR of the article, 12 variants):
%   child_cif = fullfile(here, 'cif', 'alpha2_Ti20Nb_Cmcm.cif');             % alpha'' Ti-20Nb, Cmcm
%   OR.parent_plane = [1 1 0];   OR.child_plane = [0 0 1];          % (110)_beta  || (001)_alpha''
%   OR.parent_dir   = [1 -1 1];  OR.child_dir   = [1 1 0];          % [1-11]_beta || [110]_alpha''
%   child_label = 'alpha''''';  name = 'alpha2_in_beta';
% hexagonal martensite alpha' instead (Burgers):
%   child_cif = fullfile(here, 'cif', 'alpha_Ti_hcp_COD9008517.cif');            % alpha-Ti, P6_3/mmc
%   OR.parent_plane = [1 1 0];   OR.child_plane = [0 0 0 1];        % (110)_beta  || (0001)_alpha'
%   OR.parent_dir   = [1 -1 1];  OR.child_dir   = [2 -1 -1 0];      % [1-11]_beta || [2-1-10]_alpha'
%   child_label = 'alpha''';  name = 'alpha_hcp_in_beta';
% or two pairs of parallel directions, or an OR predefined in MTEX:
%   OR = struct('parent_dirs', [1 1 0; 1 -1 1], 'child_dirs', [0 0 1; 1 1 0]);
%   OR = 'Burgers';   (also 'KurdjumovSachs', 'NishiyamaWassermann', 'Pitsch', 'Bain', 'GreningerTrojano',
%                      'PitschSchrader', 'ShojiNishiyama')

% --- orientation of the parent grain (optional), for the variants in the sample frame
parent_orientation = [];
%   []                                    only the parent crystal frame
%   [phi1 Phi phi2]                       Bunge Euler angles of the parent (deg)
%   fullfile(here, 'results', 'X', 'orientation.mat')   result of the orientation search
%                                         (grain 1; a search with these variants or of the parent phase)

parent_label = 'beta';                % names in the tables and pictures ('' = from the CIF)
child_label  = 'omega';

name = 'omega_in_beta';                % output: variants/<name>.mat, _variants.csv, _polefigures.png

%% Start MTEX =============================================================
if ~exist('crystalSymmetry', 'class')
    if isempty(mtex_path)
        error('Set mtex_path (the folder with startup_mtex.m) or start MTEX first.');
    end
    addpath(mtex_path);
    startup_mtex;
end

%% Variants (MTEX) ========================================================
[V, p2c, info] = variants_mtex(parent_cif, child_cif, OR);
csP = info.csParent;  csC = info.csChild;
N = size(V, 3);
fprintf('\n%s -> %s: %d variants\n', csP.mineral, csC.mineral, N);
fprintf('parent -> child misorientation: %.2f deg (Euler %.2f %.2f %.2f deg)\n', ...
    angle(p2c) / degree, p2c.phi1 / degree, p2c.Phi / degree, p2c.phi2 / degree);

%% Orientations of the variants ===========================================
pc = read_cif(parent_cif);  cc = read_cif(child_cif);
[~, Bp, Ap] = calculate_reciprocal_lattice(struct('a', pc.cell(1:3), 'u', pc.cell(4:6)), [0 0 0], [1 0 0 0]);
[~, Bc, Ac] = calculate_reciprocal_lattice(struct('a', cc.cell(1:3), 'u', cc.cell(4:6)), [0 0 0], [1 0 0 0]);
Rp = point_group_rotations(pc);  Rc = point_group_rotations(cc);

% parent orientation in the sample frame (crystal -> sample), if given
Lp = [];
if ~isempty(parent_orientation)
    if ischar(parent_orientation) || isstring(parent_orientation)
        r = load(parent_orientation, 'grains_saved');
        Rz = @(a) [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
        Lp = Rz(r.grains_saved(1).phi0) * quat2dcm(r.grains_saved(1).q)';
    else
        e = parent_orientation;
        Rz = @(a) [cosd(a) sind(a) 0; -sind(a) cosd(a) 0; 0 0 1];      % passive
        Rx = @(a) [1 0 0; 0 cosd(a) sind(a); 0 -sind(a) cosd(a)];
        Lp = (Rz(e(3)) * Rx(e(2)) * Rz(e(1)))';                         % g' (g: sample -> crystal)
    end
    fprintf('parent orientation in the sample frame: Euler (Bunge) %.2f %.2f %.2f deg\n', euler_bunge(Lp'));
end

c001 = strings(N, 1);  c100 = c001;  err = zeros(N, 1);
eul_parent = zeros(N, 3);  eul_sample = nan(N, 3);  mis1 = zeros(N, 1);
for j = 1:N
    c001(j) = uvw_str((V(:,:,j) * ([0 0 1] * Ac)')', Ap);
    c100(j) = uvw_str((V(:,:,j) * ([1 0 0] * Ac)')', Ap);
    if isstruct(OR) && isfield(OR, 'parent_plane')
        err(j) = or_error(V(:,:,j), Rp, Rc, three(OR.parent_plane, 'hkl') * Bp, three(OR.child_plane, 'hkl') * Bc, ...
            three(OR.parent_dir, 'uvw') * Ap, three(OR.child_dir, 'uvw') * Ac);
    else
        err(j) = NaN;
    end
    eul_parent(j, :) = euler_bunge(V(:,:,j)');                 % sample = parent crystal frame
    if ~isempty(Lp), eul_sample(j, :) = euler_bunge((Lp * V(:,:,j))'); end
end
% misorientations between the variants (smallest over the child symmetry)
M = zeros(N);
for i = 1:N
    for j = 1:N
        M(i, j) = disorientation_angle(V(:,:,i), V(:,:,j), Rc);
    end
end
mis1 = M(:, 1);

T = table((1:N)', c001, c100, eul_parent(:,1), eul_parent(:,2), eul_parent(:,3), mis1, err, ...
    'VariableNames', {'variant', 'child_001_along_parent', 'child_100_along_parent', ...
    'phi1_parent_frame', 'Phi_parent_frame', 'phi2_parent_frame', 'misorientation_to_variant_1', 'OR_check_deg'});
if ~isempty(Lp)
    T.phi1_sample = eul_sample(:,1);  T.Phi_sample = eul_sample(:,2);  T.phi2_sample = eul_sample(:,3);
end
fprintf('\n');  disp(T);
angs = uniquetol(M(triu(true(N), 1)), 1e-3);
fprintf('Misorientations between variants (deg): %s\n', strjoin(compose('%.2f', angs), ', '));

%% Pole figures (MTEX), in the frame of the parent crystal =================
hC = Miller({0,0,1}, {1,0,0}, {0,1,0}, csC, 'uvw');
cols = lines(N);
mf = newMtexFigure('layout', [1 3]);
for j = 1:N
    opts = {'MarkerFaceColor', cols(j,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'label', {num2str(j)}, 'upper'};
    if j > 1, opts = [opts, {'add2all'}]; end %#ok<AGROW>
    plotPDF(info.oriChild(j), hC, opts{:});
end
hold all
pv = vector3d(symmetrise(Miller({1,1,0}, csP, 'uvw'), 'unique'));     % parent frame = plot frame
plot(pv, 'Marker', 's', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', 14, 'add2all', 'upper');
hold off
drawNow(mf, 'figSize', 'large');                          % lay out the pole figures
mtexTitle(mf.children(1), '[001]');
set(mf.parent, 'Name', sprintf('%s variants in %s, pole figures in the parent frame (squares = parent <110>)', ...
    csC.mineral, csP.mineral));

%% 3D pictures of the variants in the parent cell ========================
if isempty(parent_label), parent_label = sprintf('%s (%s)', pc.name, pc.space_group); end
if isempty(child_label),  child_label  = sprintf('%s (%s)', cc.name, cc.space_group); end
[fig3d1, fig3d2, fig3d3] = plot_variants_3d(V, pc, cc, OR, parent_label, child_label);   % fig3d3: viewer, turn with the mouse

%% Save ===================================================================
out_dir = fullfile(here, 'variants');
if ~isfolder(out_dir), mkdir(out_dir); end
out_file = fullfile(out_dir, [name '.mat']);
p2c_euler = [p2c.phi1 p2c.Phi p2c.phi2] / degree;
OR_saved = OR;
if isa(OR, 'orientation'), OR_saved = struct('mtex_euler', p2c_euler); end
variants_table = T;  misorientations = M;  parent_L = Lp;
parent_cif = relpath(parent_cif, here);  child_cif = relpath(child_cif, here);   % portable
save(out_file, 'V', 'parent_cif', 'child_cif', 'OR_saved', 'p2c_euler', 'variants_table', 'misorientations', 'parent_L', ...
    'parent_label', 'child_label');
writetable(T, fullfile(out_dir, [name '_variants.csv']));
saveFigure(fullfile(out_dir, [name '_polefigures.png']), mf.parent);
exportgraphics(fig3d1, fullfile(out_dir, [name '_3d_overview.png']), 'Resolution', 150);
exportgraphics(fig3d2, fullfile(out_dir, [name '_3d_variants.png']), 'Resolution', 150);
exportgraphics(fig3d3, fullfile(out_dir, [name '_3d_variant1.png']), 'Resolution', 150);
fprintf('\nSaved %s (+ _variants.csv, _polefigures.png, _3d_*.png)\nuse: phase = phase_from_cif(child_cif, ''variants'', ''%s'');\n', ...
    out_file, out_file);

%% ------------------------------------------------------------------------
function p = relpath(p, base)
% path relative to base (e.g. 'cif/x.cif'), if p lies in base
if startsWith(p, [base filesep]), p = p(numel(base) + 2:end); end
end

function v = three(v, type)
% 4-index (hkil / UVTW) -> 3-index (hkl / uvw)
if numel(v) == 4
    if strcmp(type, 'hkl'), v = v([1 2 4]); else, v = [v(1) - v(3), v(2) - v(3), v(4)]; end
end
end

function s = uvw_str(v, A)
% closest integer direction [uvw] (indices <= 4) of the Cartesian vector v, A = direct basis (rows)
[h, k, l] = ndgrid(-4:4);
H = [h(:) k(:) l(:)];  H = H(any(H, 2), :);
H = H(gcd(gcd(abs(H(:,1)), abs(H(:,2))), abs(H(:,3))) == 1, :);
D = H * A;
[~, i] = max((D ./ vecnorm(D, 2, 2)) * (v(:) / norm(v)));
s = sprintf('[%d %d %d]', H(i, :));
end

function e = or_error(V, Rp, Rc, np, nc, dp, dc)
% largest angle between the given parallel pairs for variant V, allowing for the
% symmetry of both phases (a variant realises an equivalent of the given pair)
e = inf;
for k = 1:size(Rp, 3)
    for m = 1:size(Rc, 3)
        W = V * Rc(:,:,m);
        a1 = acosd(min(1, abs(dot((Rp(:,:,k) * np') / norm(np), W * nc' / norm(nc)))));
        a2 = acosd(min(1, abs(dot((Rp(:,:,k) * dp') / norm(dp), W * dc' / norm(dc)))));
        e = min(e, max(a1, a2));
    end
end
end

function a = disorientation_angle(Vi, Vj, Rc)
% smallest rotation angle between two child orientations (child -> parent
% matrices), over the child symmetry on both sides
a = inf;
for p = 1:size(Rc, 3)
    for q = 1:size(Rc, 3)
        D = (Vi * Rc(:,:,p))' * (Vj * Rc(:,:,q));
        a = min(a, acosd(max(-1, min(1, (trace(D) - 1) / 2))));
    end
end
end

function [V, p2c, info] = variants_mtex(parent_cif, child_cif, OR)
%VARIANTS_MTEX  Variants of a product phase in a parent crystal, computed with MTEX.
%
%   [V, p2c, info] = variants_mtex(parent_cif, child_cif, OR)
%
%   parent_cif, child_cif : CIF files of the parent (e.g. beta) and the child
%                           (e.g. alpha'') phase
%   OR : the orientation relationship, one of
%        - two pairs of parallel PLANE and DIRECTION (as usually written,
%          e.g. (110)_beta || (001)_alpha'', [1-11]_beta || [110]_alpha''):
%            OR.parent_plane = [1 1 0];   OR.child_plane = [0 0 1];
%            OR.parent_dir   = [1 -1 1];  OR.child_dir   = [1 1 0];
%        - two pairs of parallel DIRECTIONS (as in VARIANTS_FROM_OR):
%            OR.parent_dirs = [1 1 0; 1 -1 1];  OR.child_dirs = [0 0 1; 1 1 0];
%        - the name of an OR predefined in MTEX, e.g. 'Burgers',
%          'KurdjumovSachs', 'NishiyamaWassermann', 'Pitsch', 'Bain',
%          'GreningerTrojano', 'PitschSchrader', 'ShojiNishiyama'
%        - an MTEX parent-to-child misorientation, e.g.
%            OR = orientation.Burgers(csParent, csChild)
%          (its crystal symmetries must be those of the two CIF files)
%        Indices: 3 (hkl / uvw); for hexagonal phases also 4 (hkil / UVTW),
%        e.g. hexagonal martensite alpha' in beta (Burgers):
%            OR.parent_plane = [1 1 0];   OR.child_plane = [0 0 0 1];
%            OR.parent_dir   = [1 -1 1];  OR.child_dir   = [2 -1 -1 0];
%
%   V    : 3 x 3 x N; V(:,:,j) turns the Cartesian frame of the child into
%          the Cartesian frame of the parent for variant j - the same as
%          VARIANTS_FROM_OR and phase.variants (frames: x || a, z || c*)
%   p2c  : the MTEX parent-to-child misorientation
%   info : csParent, csChild (MTEX crystalSymmetry), oriChild (the variants
%          as MTEX orientations for the parent orientation = identity)
%
%   MTEX (https://mtex-toolbox.github.io) must be on the MATLAB path and
%   started (startup_mtex). The crystal symmetries are built from the space
%   group and cell read by READ_CIF (MTEX's own CIF reader does not know the
%   newer CIF keys such as _space_group_name_H-M_alt), with the alignment
%   X || a, Z || c* used by this tool. MTEX only needs the symmetry and the
%   cell; atoms and structure factors still come from PHASE_FROM_CIF.
%
%   The variants are variants(p2c, parent orientation = identity): all
%   child orientations p * inv(p2c) for the rotations p of the parent,
%   with those equivalent by the child symmetry removed.

if ~exist('crystalSymmetry', 'class')
    error('variants_mtex:noMTEX', ['MTEX is not on the MATLAB path. Add it and start it, e.g.\n' ...
        '  addpath(''C:\\...\\mtex-6.1.1''); startup_mtex']);
end

pc = read_cif(parent_cif);
cc = read_cif(child_cif);
csP = cs_from_cif(pc, 'parent');
csC = cs_from_cif(cc, 'child');

if ischar(OR) || isstring(OR)
    % an OR predefined in MTEX, by name, e.g. 'Burgers', 'KurdjumovSachs'
    try
        p2c = orientation.(char(OR))(csP, csC);
    catch err
        error('variants_mtex:ORname', ['MTEX has no predefined OR "%s" for these phases (%s). Known: Bain, ' ...
            'KurdjumovSachs, NishiyamaWassermann, Pitsch, GreningerTrojano, Burgers, PitschSchrader, ' ...
            'ShojiNishiyama.'], char(OR), err.message);
    end
elseif isa(OR, 'orientation')
    p2c = OR;
    if ~(p2c.CS.Laue == csP.Laue && p2c.SS.Laue == csC.Laue)
        error('variants_mtex:symmetry', 'The symmetries of the given OR differ from those of the CIF files.');
    end
    p2c = orientation(rotation(p2c), csP, csC);          % use the symmetries built from the CIFs
elseif all(isfield(OR, {'parent_plane', 'child_plane', 'parent_dir', 'child_dir'}))
    p2c = orientation.map(miller(OR.parent_plane, csP, 'hkl'), miller(OR.child_plane, csC, 'hkl'), ...
                          miller(OR.parent_dir, csP, 'uvw'), miller(OR.child_dir, csC, 'uvw'));
elseif all(isfield(OR, {'parent_dirs', 'child_dirs'}))
    p2c = orientation.map(miller(OR.parent_dirs(1,:), csP, 'uvw'), miller(OR.child_dirs(1,:), csC, 'uvw'), ...
                          miller(OR.parent_dirs(2,:), csP, 'uvw'), miller(OR.child_dirs(2,:), csC, 'uvw'));
else
    error('variants_mtex:OR', ['OR must have the fields parent_plane, child_plane, parent_dir, child_dir, ' ...
        'or parent_dirs, child_dirs, or be an MTEX orientation.']);
end

oriChild = variants(p2c, orientation.id(csP));
V = zeros(3, 3, length(oriChild));
for j = 1:length(oriChild)
    V(:,:,j) = matrix(oriChild(j));                      % child Cartesian -> parent Cartesian
end
info = struct('csParent', csP, 'csChild', csC, 'oriChild', oriChild);
end

%% ------------------------------------------------------------------------
function cs = cs_from_cif(cif, label)
% MTEX crystalSymmetry from the space group and the cell of READ_CIF
sg = strrep(strtrim(cif.space_group), ' ', '');
sg = regexprep(sg, ':.*$', '');                          % origin / setting suffix
try
    cs = crystalSymmetry(sg, cif.cell(1:3), cif.cell(4:6) * pi / 180, 'X||a', 'Z||c*', 'mineral', cif.name);
catch err
    error('variants_mtex:spacegroup', 'MTEX does not know the space group "%s" of the %s CIF (%s).', ...
        cif.space_group, label, err.message);
end
end

function m = miller(v, cs, type)
% Miller indices as an MTEX Miller object: 3 indices (hkl / uvw) or, for
% hexagonal and trigonal phases, 4 indices (hkil / UVTW, with i = -(h+k), T = -(U+V))
if numel(v) == 4
    if abs(v(3) + v(1) + v(2)) > 1e-9
        error('variants_mtex:index', 'Four-index notation needs the 3rd index = -(1st + 2nd): [%s].', num2str(v));
    end
    if strcmp(type, 'hkl'), type = 'hkil'; else, type = 'UVTW'; end
end
v = num2cell(v);
m = Miller(v{:}, cs, type);
end

function phase = phase_from_cif(cif_file, varargin)
%PHASE_FROM_CIF  Phase for the orientation search defined by a CIF file.
%
%   phase = phase_from_cif(cif_file)
%       single crystal of the phase (no variants)
%
%   phase = phase_from_cif(cif_file, 'parent', parent_cif, 'OR', OR)
%       all variants of the phase formed in a parent crystal with the given
%       orientation relationship (see VARIANTS_FROM_OR); the orientation of
%       the PARENT is searched
%
%   phase = phase_from_cif(cif_file, 'variants', variants_file)
%   phase = phase_from_cif(cif_file, 'parent', parent_cif, 'variants', V)
%       the same with variants computed elsewhere, e.g. with MTEX by
%       make_variants_mtex.m: variants_file is the .mat file it writes
%       (V and the parent CIF are taken from it), or V a 3x3xN array
%       (V(:,:,j): child Cartesian -> parent Cartesian, as VARIANTS_FROM_OR)
%
%   Further options (name, value):
%     'hkl_max'  reflections -hkl_max..hkl_max are used (default 6);
%                'auto' = all reflections that can diffract at the wavelength
%                (d >= lambda/2); more reflections at high angles, slower search
%     'cell'     [a b c alpha beta gamma] overriding the cell of the CIF
%                (e.g. lattice parameters measured on the sample)
%     'name'     label used in the output
%
%   The returned struct has the fields used by PHASE_REFLECTIONS and
%   SEARCH_ORIENTATION:
%     name, cell, elements, Z (atomic numbers), positions [x y z occ type],
%     hkl_max, variants (3x3xN), parent_rotations (3x3xM, proper rotations
%     of the crystal whose orientation is searched), cif (the parsed file)

ip = inputParser;
ip.addParameter('parent', '');
ip.addParameter('OR', []);
ip.addParameter('variants', []);
ip.addParameter('hkl_max', 6);
ip.addParameter('cell', []);
ip.addParameter('name', '');
ip.parse(varargin{:});
o = ip.Results;

cif = read_cif(cif_file);
if ~isempty(o.cell), cif.cell = o.cell; end
[positions, elements] = expand_positions(cif);
Z = zeros(1, numel(elements));
for k = 1:numel(elements)
    zk = element_Z(elements{k});
    if isempty(zk)
        warning('phase_from_cif:element', ['Unknown element "%s" in %s; these atoms are ' ...
            'ignored (scattering factor 0).'], elements{k}, cif_file);
        zk = 0;
    end
    Z(k) = zk;
end

phase.cif       = cif;
phase.cell      = cif.cell;
phase.elements  = elements;
phase.Z         = Z;
phase.positions = positions;
phase.hkl_max   = o.hkl_max;

% variants given (e.g. from MTEX): a .mat file of make_variants_mtex or a 3x3xN array
V_given = [];
if ~isempty(o.variants)
    if ischar(o.variants) || isstring(o.variants)
        f = load(o.variants, 'V', 'parent_cif');
        V_given = f.V;
        if isempty(o.parent)
            o.parent = f.parent_cif;
            if ~isfile(o.parent)                % moved: look in cif/ next to the variants folder
                [~, nm, ext] = fileparts(f.parent_cif);
                o.parent = fullfile(fileparts(fileparts(char(o.variants))), 'cif', [nm ext]);
            end
        end
    else
        V_given = o.variants;
    end
    if isempty(o.parent)
        error('phase_from_cif:parent', 'With ''variants'' as an array the parent CIF is needed: ''parent'', parent_cif.');
    end
    for j = 1:size(V_given, 3)
        if norm(V_given(:,:,j)' * V_given(:,:,j) - eye(3)) > 1e-6 || det(V_given(:,:,j)) < 0
            error('phase_from_cif:variants', 'Variant %d is not a rotation matrix.', j);
        end
    end
end

if isempty(o.parent)
    phase.variants = eye(3);
    phase.parent_rotations = point_group_rotations(cif);
    label = sprintf('%s (%s), single crystal', cif.name, cif.space_group);
else
    if isempty(V_given) && (isempty(o.OR) || ~all(isfield(o.OR, {'parent_dirs', 'child_dirs'})))
        error('phase_from_cif:OR', ['With ''parent'' the orientation relationship ''OR'' is needed: ' ...
            'OR.parent_dirs = [u v w; u v w], OR.child_dirs = [u v w; u v w].']);
    end
    parent = read_cif(o.parent);
    if isempty(V_given)
        phase.variants = variants_from_OR(parent, cif, o.OR);
    else
        phase.variants = V_given;
    end
    phase.parent_rotations = point_group_rotations(parent);
    label = sprintf('%s (%s), %d variants in %s (%s)', cif.name, cif.space_group, ...
        size(phase.variants, 3), parent.name, parent.space_group);
end
if isempty(o.name), phase.name = label; else, phase.name = o.name; end
fprintf('Phase: %s\n  %d atoms in the unit cell (%s), %d symmetry operations\n', phase.name, ...
    size(positions, 1), strjoin(elements, ', '), size(cif.W, 3));
end

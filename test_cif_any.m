%% TEST_CIF_ANY  CIF files of all crystal systems (cif/test_set/, from COD).
%
% For every file: read it, expand the atoms, count the proper rotations and
% check that all reflections forbidden by the space group (centring, glide
% planes, screw axes) are absent from the reflection list. If the file lists
% its symmetry operations AND a space-group symbol, the operations taken
% from the symbol (syminfo.lib) are compared with those of the file.

clear;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
d = fullfile(here, 'cif', 'test_set');
odd = @(x) mod(x, 2) ~= 0;

% file, phase, expected proper rotations, reflections that must be absent
T = {
 'COD_1000041.cif', 'NaCl',          24, @(h,k,l) ~((~odd(h) & ~odd(k) & ~odd(l)) | (odd(h) & odd(k) & odd(l)))
 'COD_1100108.cif', 'alpha-Fe',      24, @(h,k,l) odd(h+k+l)
 'COD_1010944.cif', 'cristobalite',  24, @(h,k,l) ~((~odd(h) & ~odd(k) & ~odd(l)) | (odd(h) & odd(k) & odd(l))) | ...
                                             (h == 0 & mod(k+l, 4) ~= 0) | (k == 0 & mod(h+l, 4) ~= 0) | (l == 0 & mod(h+k, 4) ~= 0)
 'COD_1575875.cif', 'Mg',            12, @(h,k,l) (h == 0 & k == 0 & odd(l)) | ((h == k | k == -2*h | h == -2*k) & odd(l))
 'COD_1000017.cif', 'corundum',       6, @(h,k,l) mod(-h+k+l, 3) ~= 0 | ((h == -k | h == 0 | k == 0) & odd(l))
 'COD_1011097.cif', 'quartz',         6, @(h,k,l) h == 0 & k == 0 & mod(l, 3) ~= 0
 'COD_1530150.cif', 'rutile',         8, @(h,k,l) (h == 0 & odd(k+l)) | (k == 0 & odd(h+l)) | (k == 0 & l == 0 & odd(h))
 'COD_1010942.cif', 'anatase',        8, @(h,k,l) odd(h+k+l) | (l == 0 & odd(h)) | (h == 0 & k == 0 & mod(l, 4) ~= 0) | ...
                                             (abs(h) == abs(k) & mod(2*h+l, 4) ~= 0)
 'COD_1010912.cif', 'baddeleyite',    2, @(h,k,l) (h == 0 & l == 0 & odd(k)) | (k == 0 & odd(l))
 'COD_2107372.cif', 'albite',         1, @(h,k,l) odd(h+k)
 'synthetic_two_blocks.cif', 'NaCl 2 blocks', 24, @(h,k,l) ~((~odd(h) & ~odd(k) & ~odd(l)) | (odd(h) & odd(k) & odd(l)))
};

all_ok = true;
fprintf('\n%-16s %-13s %-15s %-30s %5s %4s %5s  %-10s %s\n', 'file', 'phase', 'space group', 'symmetry from', ...
    'atoms', 'rot', 'refl', 'extinction', 'ops from symbol = ops in file');
for t = 1:size(T, 1)
    f = fullfile(d, T{t,1});
    c = read_cif(f);
    pos = expand_positions(c);
    R = point_group_rotations(c);
    ph = phase_from_cif(f, 'hkl_max', 6);
    [~, hkl] = phase_reflections(ph, 0.7093);
    bad = T{t,4}(hkl(:,1), hkl(:,2), hkl(:,3));
    ext_ok = ~any(bad);
    rot_ok = size(R, 3) == T{t,3};

    % operations from the symbol vs. operations in the file
    cmp = 'n/a';
    if startsWith(c.symmetry_source, 'operations listed')
        try
            warning('off', 'spacegroup_ops:setting');
            [ops2, sg] = spacegroup_ops('hall', c.hall, 'hm', c.space_group, 'number', c.sg_number, 'cell', c.cell);
            if strcmp(sg.matched_by, 'number')
                cmp = "n/a, setting " + string(c.space_group) + " not in library";
            else
                cmp = string(same_ops(c, ops2)) + " (" + sg.matched_by + ": " + sg.xhm + ")";
            end
        catch err
            cmp = "not in library (" + string(err.identifier) + ")";
        end
    end
    ok = ext_ok && rot_ok;
    all_ok = all_ok && ok;
    src = c.symmetry_source; if numel(src) > 30, src = [src(1:27) '...']; end
    fprintf('%-16s %-13s %-15s %-30s %5d %4d %5d  %-10s %s\n', T{t,1}, T{t,2}, c.space_group, src, ...
        size(pos, 1), size(R, 3), size(hkl, 1), string(ext_ok) + "/" + string(rot_ok), cmp);
    if ~ext_ok
        fprintf('      forbidden reflections present, e.g. %s\n', mat2str(hkl(find(bad, 3), :)));
    end
end
fprintf('\nALL TESTS %s\n', string(all_ok));

function same = same_ops(c, ops2)
% compare two sets of operations (W, w mod 1)
key = @(W, w) round([W(:); mod(w(:), 1)]' * 1e6) / 1e6;
A = zeros(size(c.W, 3), 12);
for k = 1:size(c.W, 3), A(k, :) = key(c.W(:,:,k), c.w(:,k)); end
B = zeros(numel(ops2), 12);
c2 = struct('ops', {ops2});
for k = 1:numel(ops2)
    [W, w] = parse(ops2{k});
    B(k, :) = key(W, w);
end
A(:, 10:12) = mod(round(A(:, 10:12) * 1e6), 1e6) / 1e6;
B(:, 10:12) = mod(round(B(:, 10:12) * 1e6), 1e6) / 1e6;
same = isequal(sortrows(unique(A, 'rows')), sortrows(unique(B, 'rows')));
end

function [W, w] = parse(op)
tmp = [tempname '.cif'];
fid = fopen(tmp, 'w');
fprintf(fid, 'data_x\n_cell_length_a 1\n_cell_length_b 1\n_cell_length_c 1\nloop_\n_symmetry_equiv_pos_as_xyz\n''%s''\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nX 0 0 0\n', op);
fclose(fid);
c = read_cif(tmp);
delete(tmp);
W = c.W(:,:,1);  w = c.w(:,1);
end

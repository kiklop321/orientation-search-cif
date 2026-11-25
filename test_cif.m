%% TEST_CIF  Checks of the CIF reader, the symmetry expansion and the variants.
%
% Reads every CIF in cif/, expands the asymmetric unit and checks the
% number of atoms, of proper rotations and the systematic extinctions
% expected for the space group. Finally the 12 alpha'' variants are
% generated from the beta and alpha'' CIFs and checked to form a set that
% the rotations of beta map onto itself.

clear;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
cifdir = fullfile(here, 'cif');
lambda = 0.7093;
ok = true;

% file, expected atoms in the cell, expected proper rotations, extinction rule (true = must be absent)
tests = {
  'beta_Ti_COD9012924.cif',          2, 24, @(h,k,l) mod(h+k+l, 2) ~= 0                          % beta-Ti, Im-3m: I-centring
  'beta_TiNb_COD1538066.cif',          4, 24, @(h,k,l) mod(h+k+l, 2) ~= 0                          % beta Ti-Nb, Im-3m, Ti/Nb mixed
  'alpha_Ti_hcp_COD9008517.cif',          2, 12, @(h,k,l) h == 0 & k == 0 & mod(l, 2) ~= 0            % alpha-Ti, P6_3/mmc: 00l, l odd
  'omega_Ti_COD9012925.cif',          3, 12, @(h,k,l) false(size(h))                              % omega-Ti, P6/mmm: none
  'alpha2_Ti20Nb_Cmcm.cif',   8,  4, @(h,k,l) mod(h+k, 2) ~= 0 | (k == 0 & mod(l, 2) ~= 0) % alpha'', Cmcm: C, c-glide
};
fprintf('%-24s %-14s %-26s %5s %5s %5s %7s  %s\n', 'file', 'space group', 'cell', 'ops', 'atoms', 'rot', 'refl', 'result');
for t = 1:size(tests, 1)
    f = fullfile(cifdir, tests{t,1});
    c = read_cif(f);
    [pos, el] = expand_positions(c);
    R = point_group_rotations(c);
    ph = phase_from_cif(f, 'hkl_max', 4);
    [~, hkl] = phase_reflections(ph, lambda);
    absent = tests{t,4}(hkl(:,1), hkl(:,2), hkl(:,3));
    pass = size(pos, 1) == tests{t,2} && size(R, 3) == tests{t,3} && ~any(absent);
    ok = ok && pass;
    fprintf('%-24s %-14s %-26s %5d %5d %5d %7d  %s\n', tests{t,1}, c.space_group, ...
        sprintf('%.3f %.3f %.3f %g %g %g', c.cell), size(c.W, 3), size(pos, 1), size(R, 3), ...
        size(hkl, 1), string(pass) + " (" + strjoin(el, ',') + ")");
end

% ---- variants beta -> alpha'' from the CIFs
% The set of all variants must be closed under the rotations of the parent:
% a rotation of beta only permutes the variants, so the reciprocal vectors
% of all variants together must map onto themselves.
OR.parent_dirs = [1 1 0; 1 -1 1];      % [110]_beta || [001]_alpha'', [1-11]_beta || [110]_alpha''
OR.child_dirs  = [0 0 1; 1 1 0];
new = phase_from_cif(fullfile(cifdir, 'alpha2_Ti20Nb_Cmcm.cif'), 'parent', fullfile(cifdir, 'beta_Ti_COD9012924.cif'), 'OR', OR);
[g_new, hkl_new] = phase_reflections(new, lambda);
closed = true;
for k = 1:size(new.parent_rotations, 3)
    closed = closed && all(ismembertol((new.parent_rotations(:,:,k) * g_new')', g_new, 1e-9, 'ByRows', true));
end
pass = size(new.variants, 3) == 12 && closed;
ok = ok && pass;
fprintf('\nbeta -> alpha'''' from CIF: %d variants, %d reflections per variant; closed under the beta rotations: %s\n', ...
    size(new.variants, 3), size(hkl_new, 1) / 12, string(closed));
fprintf('\nALL TESTS %s\n', string(ok));

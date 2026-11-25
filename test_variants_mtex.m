%% TEST_VARIANTS_MTEX  Variants from MTEX against the variants of this tool.
%
% Needs MTEX (set mtex_path). For two orientation relationships the
% variants computed by MTEX (src/variants_mtex.m) are compared with those of
% VARIANTS_FROM_OR: every MTEX variant must be one of ours, up to the
% rotations of the child (same variant), and the numbers must agree.
%   1. beta -> alpha'' (Cmcm), the OR of the article, given as planes + directions
%   2. beta -> alpha (hexagonal, P6_3/mmc), Burgers OR as predefined in MTEX
%      (orientation.Burgers) - checks also the axis convention of a
%      hexagonal cell (X || a, Z || c*)
%   3. the same Burgers OR given in 4-index notation,
%      (110)_beta || (0001)_alpha, [1-11]_beta || [2-1-10]_alpha

clear;
here = fileparts(mfilename('fullpath'));
if isempty(here) || ~isfile(fullfile(here, 'test_variants_mtex.m')), here = fileparts(which('test_variants_mtex.m')); end
addpath(fullfile(here, 'src'));
mtex_path = '';                  % folder with startup_mtex.m ('' = MTEX already started)
if ~exist('crystalSymmetry', 'class')
    if isempty(mtex_path), error('Set mtex_path (the folder with startup_mtex.m) or start MTEX first.'); end
    addpath(mtex_path);  startup_mtex;
end

beta  = fullfile(here, 'cif', 'beta_Ti_COD9012924.cif');
cases = { ...
    'beta -> alpha'''' (article OR)', fullfile(here, 'cif', 'alpha2_Ti20Nb_Cmcm.cif'), ...
        struct('parent_plane', [1 1 0], 'child_plane', [0 0 1], 'parent_dir', [1 -1 1], 'child_dir', [1 1 0]), ...
        struct('parent_dirs', [1 1 0; 1 -1 1], 'child_dirs', [0 0 1; 1 1 0]);
    'beta -> alpha hcp (Burgers, MTEX)', fullfile(here, 'cif', 'alpha_Ti_hcp_COD9008517.cif'), ...
        'Burgers', ...
        struct('parent_dirs', [1 1 0; -1 1 1], 'child_dirs', [0 0 1; 1 0 0]);   % [2-1-10] = [100]
    'beta -> alpha hcp (Burgers, 4 indices)', fullfile(here, 'cif', 'alpha_Ti_hcp_COD9008517.cif'), ...
        struct('parent_plane', [1 1 0], 'child_plane', [0 0 0 1], 'parent_dir', [1 -1 1], 'child_dir', [2 -1 -1 0]), ...
        struct('parent_dirs', [1 1 0; 1 -1 1], 'child_dirs', [0 0 1; 1 0 0]) };
ok = true;
for c = 1:size(cases, 1)
    pc = read_cif(beta);  cc = read_cif(cases{c, 2});
    OR = cases{c, 3};
    if ischar(OR)                                     % MTEX predefined relationship
        [~, ~, inf0] = variants_mtex(beta, cases{c, 2}, cases{c, 4});
        OR = orientation.Burgers(inf0.csParent, inf0.csChild);
    end
    Vm = variants_mtex(beta, cases{c, 2}, OR);
    Vo = variants_from_OR(pc, cc, cases{c, 4});
    Rc = point_group_rotations(cc);
    found = zeros(1, size(Vm, 3));
    for j = 1:size(Vm, 3)
        for i = 1:size(Vo, 3)
            for m = 1:size(Rc, 3)
                if norm(Vm(:,:,j) - Vo(:,:,i) * Rc(:,:,m)) < 1e-6, found(j) = i; end
            end
        end
    end
    pass = size(Vm, 3) == size(Vo, 3) && all(found > 0) && numel(unique(found)) == numel(found);
    ok = ok && pass;
    fprintf('%-38s MTEX %2d variants, own %2d; every MTEX variant one of ours: %s\n', cases{c, 1}, ...
        size(Vm, 3), size(Vo, 3), string(pass));
end
fprintf('\nALL TESTS %s\n', string(ok));

%% ORIENTATION SEARCH from one 2D X-ray diffraction frame
%
% Finds the Bragg spots on an image-plate frame and the orientation of a
% crystal, defined by a CIF file, that reproduces them.
%
%   1. read the frame                     (Rigaku R-AXIS .img)
%   2. find the spots                     (src/find_bragg_peaks.m)
%   3. build the reflections of the phase (cif/*.cif, src/phase_from_cif.m, src/phase_reflections.m)
%   4. search the orientation             (src/search_orientation.m)
%   5. assign the spots, plot, save       (src/match_spots.m, src/plot_overlay.m)
%
% Everything that can be changed is in section 0. Output goes to
% results/<name>/ (see README.md).
%
% Method: Salata et al., J. Appl. Cryst. 59 (2026), Sections 3.1 and 4.3.

clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'), ...
        fullfile(here, 'external', 'RigakuRapidUtils'), fullfile(here, 'external', 'sphere'));

%% 0. Settings ============================================================

% --- input frame
name        = 'Ti20Nb0O_alpha2_single';    % name of the output folder (results/<name>/)
image_file  = fullfile(here, 'data', 'Ti20Nb0O.img');   % the frame (.img); example frames: data/README.md

% --- instrument
S.wavelength = 0.7093;              % (A), Mo K-alpha
beam_center  = [1274, 3247];        % detector pixel [row, column] of the direct beam, found
                                    % by making the 2theta peaks as narrow as possible; the
                                    % nominal value in the .img header, [1287, 3246], fits
                                    % the spots much worse. The sample-detector distance
                                    % (127.26 mm = radius of the image plate) is read from the header.
S.refine_beam_center = false;       % true = refine beam_center to the matched spots (no calibrant
                                    % needed; start within ~10 px), see README, step 6
S.beam_center_scan   = 20;          % before that, coarse scan +-this many px (0 = none)

% --- angles
S.rotation_range = [-30 30];        % rotation of the sample during the exposure,
                                    % relative to its centre phi0 (deg)
S.phi0_range     = [0 360];         % where the centre phi0 is searched (deg);
                                    % keep [0 360]: phi0 is not the goniometer angle; narrow
                                    % only around a phi0 found before for the same crystal
S.phi0_step      = 0.25;            % (deg)
S.tth_range      = [0 180];         % only predicted spots with 2theta in this
                                    % range are matched (deg)
plot_max_tth     = 180;             % predicted spots drawn up to this 2theta (deg); all = 180

% --- spot search (on log(log(I+1)) of the frame)
P.min_height           = 1.7;       % minimum peak height
P.min_prominence_map   = 0.09;      % minimum prominence in the (2theta, beta) map
P.min_prominence_pixel = 0.01;      % minimum prominence in detector pixels
P.min_distance_px      = 6;         % closer spots are merged (px)
P.refine_radius_px     = 3;         % map hits are moved to the brightest pixel within (px)
P.tth_beta_map         = [0.6, pi/2, 1.2, 2*pi];   % map centre and size [2theta beta d2theta dbeta] (rad)

% --- orientation search
S.mesh_degree   = [7 10];           % icosahedron subdivisions, stage 1 and 2
S.tolerance     = [0.25 0.3];       % max |d2theta| and |dbeta| for a match, stage 1 and 2 (deg)
S.unmatched_tol = 2;                % tolerance for counting unexplained spots (deg)
S.refine_cell   = 10;               % max. times the lattice parameters are refined to the matched
                                    % spots and the orientation re-optimised; stops when the
                                    % cell no longer changes (0 = keep the CIF cell)
S.n_grains      = 1;                % MAXIMUM number of grains: 1 = one crystal; > 1 = after
                                    % each grain the search is repeated on the unexplained spots
                                    % and stops at the first grain with < S.min_spots_per_grain
S.min_spots_per_grain = 15;         % a further grain must explain at least this many spots

% --- phase: a CIF file (cell, symmetry operations, atoms)
cif_file = fullfile(here, 'cif', 'alpha2_Ti20Nb_Cmcm.cif');          % alpha'' Ti-20Nb (Cmcm)
phase = phase_from_cif(cif_file, 'hkl_max', 6);
%   Options are given as name-value pairs, in any order and combination, e.g.
%     phase = phase_from_cif(cif_file);                          % defaults (hkl_max = 6)
%     phase = phase_from_cif(cif_file, 'hkl_max', 9);
%     phase = phase_from_cif(cif_file, 'hkl_max', 6, 'cell', [3.30 3.30 3.30 90 90 90]);
%
%   'cell', [a b c alpha beta gamma]  lattice parameters of the sample (A, deg),
%                                     replace those of the CIF
%   'hkl_max', n   which reflections are predicted (default 6):
%                  1. all (h k l) with EACH index from -n to n (not the sum);
%                  2. only those that can diffract, d >= wavelength/2;
%                  3. without the extinct ones (structure factor F = 0 from the
%                     atoms of the CIF, e.g. h+k+l odd for a body-centred cell).
%                  The list is COMPLETE up to
%                     2theta = 2*asin(wavelength*(n+1) / (2*max(a,b,c)))
%                  (a reflection needing index n+1 has d <= max(a,b,c)/(n+1));
%                  above it some weak high-index reflections are missing.
%                  For Mo K-alpha and n = 6:  beta-Ti (a = 3.30 A)  complete to ~98 deg,
%                                             alpha'' (b = 4.86 A)  complete to ~62 deg.
%                  6 is enough for the search: it relies on the strong, low-angle
%                  spots (n = 9 gives the same orientations). A larger n only adds
%                  high-angle predictions (mostly green circles) and is slower.
%   'hkl_max', 'auto'  all reflections that can diffract at this wavelength,
%                  n = ceil(2*max(a,b,c)/wavelength): 10 for beta-Ti, 14 for alpha''
% Special case - all variants of a product phase in a parent crystal
% (the orientation of the parent is then searched), e.g. alpha'' in beta:
%   OR.parent_dirs = [1 1 0; 1 -1 1];  OR.child_dirs = [0 0 1; 1 1 0];
%   phase = phase_from_cif(cif_file, 'parent', fullfile(here, 'cif', 'beta_Ti_COD9012924.cif'), 'OR', OR);
% or with the variants made by MTEX in make_variants_mtex.m (OR as planes + directions):
%   phase = phase_from_cif(cif_file, 'variants', fullfile(here, 'variants', 'alpha2_in_beta.mat'));
S.search_region = phase.parent_rotations;   % fundamental region of the crystal symmetry
                                            % ('octant' = as in the original scripts)

%% 1-5. Read the frame, find the spots and the orientation, save =========
out_dir = fullfile(here, 'results', name);
out = orientation_pipeline(image_file, beam_center, phase, S, P, out_dir, name, plot_max_tth);

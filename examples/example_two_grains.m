%% EXAMPLE_TWO_GRAINS  Orientation search for several grains in one frame.
%
% The frame Ti24Nb0O (beta Ti-24Nb) contains two beta grains. With
% S.n_grains > 1 the search is repeated on the spots that the grains found
% so far do not explain; it stops when a further grain explains fewer than
% S.min_spots_per_grain spots. Grains are marked in different colours in
% results/Ti24Nb0O_two_grains/orientation_overlay.png and numbered in the
% column 'grain' of spots_matched.csv.

clear; close all;
here = fileparts(mfilename('fullpath'));
here = fileparts(here);                  % the tool folder (this script is in examples/)
addpath(fullfile(here, 'src'), fullfile(here, 'external', 'RigakuRapidUtils'), fullfile(here, 'external', 'sphere'));

S.wavelength = 0.7093;  beam_center = [1274, 3247];
S.rotation_range = [-30 30];  S.phi0_range = [0 360];  S.phi0_step = 0.25;  S.tth_range = [0 180];
S.mesh_degree = [7 10];  S.tolerance = [0.25 0.3];  S.unmatched_tol = 2;
S.refine_cell = 10;
S.n_grains = 3;                    % look for AT MOST three grains (two are found, see above)
S.min_spots_per_grain = 15;
P.min_height = 1.7;  P.min_prominence_map = 0.09;  P.min_prominence_pixel = 0.01;
P.min_distance_px = 6;  P.refine_radius_px = 3;  P.tth_beta_map = [0.6, pi/2, 1.2, 2*pi];

phase = phase_from_cif(fullfile(here, 'cif', 'beta_Ti_COD9012924.cif'), 'hkl_max', 6);   % beta-Ti
S.search_region = phase.parent_rotations;

out = orientation_pipeline(fullfile(here, 'data', 'two_grains', 'Ti24Nb0O.img'), beam_center, phase, S, P, ...
    fullfile(here, 'results', 'Ti24Nb0O_two_grains'), 'Ti24Nb0O - beta-Ti', 180);   % last argument: plot up to 2theta = 180 deg

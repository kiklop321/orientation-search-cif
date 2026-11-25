%% TEST_BETA_TI  Orientation search with the beta-Ti CIF on four beta frames and a control.
%
% Ti-24Nb and Ti-28Nb consist of the beta phase only (article, Section 4.1),
% so a single beta crystal should explain most of their spots; four frames
% of these alloys with one grain each are tested. Ti-20Nb is
% fully transformed to alpha'' and serves as a negative control: beta should
% explain only few spots there.
%
% Phase: beta-Ti, Im-3m, a = 3.282 A (COD 9012924). After the first search
% the lattice parameter is refined to the matched spots (S.refine_cell).

clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'), ...
        fullfile(here, 'external', 'RigakuRapidUtils'), fullfile(here, 'external', 'sphere'));
data = fullfile(here, 'data');                                    % frames used here

% same settings as run_orientation_search.m
S.wavelength = 0.7093;  beam_center = [1274, 3247];
S.rotation_range = [-30 30];  S.phi0_range = [0 360];  S.phi0_step = 0.25;  S.tth_range = [0 180];
S.mesh_degree = [7 10];  S.tolerance = [0.25 0.3];  S.unmatched_tol = 2;
S.refine_cell = 10;                    % refine the lattice parameter until it converges (max. 10 times)
P.min_height = 1.7;  P.min_prominence_map = 0.09;  P.min_prominence_pixel = 0.01;
P.min_distance_px = 6;  P.refine_radius_px = 3;  P.tth_beta_map = [0.6, pi/2, 1.2, 2*pi];
plot_max_tth = 180;

phase = phase_from_cif(fullfile(here, 'cif', 'beta_Ti_COD9012924.cif'), 'hkl_max', 6);
S.search_region = phase.parent_rotations;

% data/beta_single_grain/: frames of Ti-24Nb and Ti-28Nb that contain one
% beta grain (checked: a second search on the unexplained spots finds no
% further grain); data/Ti20Nb0O.img: alpha'' control
frames = {'beta_single_grain/Ti28Nb0O', 'beta_single_grain/Ti28Nb_1', 'beta_single_grain/Ti24Nb_1', ...
          'beta_single_grain/Ti24Nb0O-50-40', 'Ti20Nb0O'};
summary = table('Size', [0 9], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'frame', 'spots', 'matched_spots', 'fraction', 'phi0', 'mean_d_tth', 'std_d_beta', 'a_refined', 'a_esd'});
for f = frames
    [~, fname] = fileparts(f{1});
    out = orientation_pipeline(fullfile(data, [f{1} '.img']), beam_center, phase, S, P, ...
        fullfile(here, 'results', [fname '_beta_Ti']), [fname ' - beta-Ti'], plot_max_tth);
    nm = numel(unique(out.T.spot));
    summary(end+1, :) = {fname, size(out.expt, 1), nm, nm / size(out.expt, 1), out.result.phi0, ...
        mean(out.T.d_tth), std(out.T.d_beta), out.cell(1), esd_a(out.cell_log)}; %#ok<SAGROW>
    close all;
end
disp(summary);
writetable(summary, fullfile(here, 'results', 'test_beta_Ti_summary.csv'));

function e = esd_a(cell_log)
% esd of a from the last refinement; NaN if the cell was not refined (< 6 matched spots)
if isempty(cell_log), e = NaN; else, e = cell_log(end, 5); end
end

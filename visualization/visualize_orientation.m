%% VISUALIZE_ORIENTATION  What the orientation (q, phi0) means in the experiment.
%
% Reads orientation.mat written by run_orientation_search / orientation_pipeline
% and shows in 3D, step by step, how the found orientation is built:
%
%   1. START       the unit cell in the reference position: a || x (beam),
%                  b || y, c || z (rotation axis of the sample); the image
%                  plate (cylinder around z), the beam and the beam centre
%   2. QUATERNION  q turns the cell so that one crystal direction n = [uvw]
%                  lies along the rotation axis z. q is a rotation by an angle
%                  about an axis in the horizontal (xy) plane - it only tilts
%                  the cell, it does not turn it about z (q(4) = 0)
%   3. PHI0        the sample then turns about z by phi0: the orientation of
%                  the crystal in the MIDDLE of the exposure
%   4. EXPOSURE    during the exposure the sample rotates from phi0 +
%                  rotation_range(1) to phi0 + rotation_range(2); every matched
%                  reflection diffracts at its own angle phi in this range and
%                  sends a beam to its spot on the plate (drawn)
%
% So the orientation = which crystal direction is along the rotation axis
% (2 angles, in q) + how far the crystal is turned about that axis (phi0).
% Crystal frame -> laboratory in the middle of the exposure:
%     L = Rz(phi0) * quat2dcm(q)'
% Laboratory: x = beam direction, z = rotation axis (up), y = x cross z side;
% the image plate is a cylinder around z with the sample in its axis.

clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here) || ~isfile(fullfile(here, 'visualize_orientation.m')), here = fileparts(which('visualize_orientation.m')); end
here = fileparts(here);                  % the tool folder (this script is in visualization/)
addpath(fullfile(here, 'src'), fullfile(here, 'external', 'RigakuRapidUtils'));

%% Settings ===============================================================
result_dir = fullfile(here, 'results', 'Ti24Nb_1_beta_Ti');   % output of run_orientation_search / test_beta_Ti
grain      = 1;             % which grain (several grains: S.n_grains > 1)
n_rays     = 12;            % number of diffracted beams drawn in panel 4 (strongest = lowest 2theta)

%% Orientation and detector =============================================
r = load(fullfile(result_dir, 'orientation.mat'), 'grains_saved', 'phase', 'S', 'image_file', 'beam_center');
G = r.grains_saved(grain);
d = RAPIDread(char(r.image_file));
d.zero = r.beam_center;
geom = RAPIDutils.getDiffractionGeometry(d);
[~, frame_name] = fileparts(r.image_file);

fig = plot_orientation_3d(G.q, G.phi0, r.phase, geom, G.T, r.S, sprintf('%s, grain %d', frame_name, grain), n_rays);
exportgraphics(fig, fullfile(result_dir, sprintf('orientation_3d_grain%d.png', grain)), 'Resolution', 150);
fprintf('Saved to %s\n', fullfile(result_dir, sprintf('orientation_3d_grain%d.png', grain)));

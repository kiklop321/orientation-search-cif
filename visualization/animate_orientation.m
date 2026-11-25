%% ANIMATE_ORIENTATION  How the vector n, the quaternion q and phi0 turn the unit cell.
%
% The search describes an orientation by
%   n     a direction FIXED IN THE CRYSTAL (e.g. [2 1 5]); it turns together
%         with the cell
%   q     the quaternion of the shortest rotation that brings n to the
%         rotation axis z: a tilt of the whole cell by the angle between n and
%         z, about the horizontal axis perpendicular to both (q(4) = 0)
%   phi0  then a turn of the whole cell about z; n stays on z, the rest of the
%         cell (a, b, c) goes round z
%
% This script shows it on the result of a search: a figure with six stages
% (start, half tilt, n on z, half phi0, final, and the final cell seen from
% above) and an animated GIF of the same movement, both saved in result_dir.
%
% Laboratory: x = beam direction, z = rotation axis of the sample.

clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here) || ~isfile(fullfile(here, 'animate_orientation.m')), here = fileparts(which('animate_orientation.m')); end
here = fileparts(here);                  % the tool folder (this script is in visualization/)
addpath(fullfile(here, 'src'));

%% Settings ===============================================================
result_dir = fullfile(here, 'results', 'Ti24Nb_1_beta_Ti');   % output of run_orientation_search / test_beta_Ti
grain      = 1;
make_gif   = true;          % also play the movement and save it as .gif and .mp4 (about 10 s)

%% Orientation ============================================================
r = load(fullfile(result_dir, 'orientation.mat'), 'grains_saved', 'phase', 'image_file');
G = r.grains_saved(grain);
[~, frame_name] = fileparts(r.image_file);
files = {};
if make_gif
    files = {'gif', fullfile(result_dir, sprintf('orientation_animation_grain%d.gif', grain)), ...
             'mp4', fullfile(result_dir, sprintf('orientation_animation_grain%d.mp4', grain))};
end
fig = plot_orientation_steps(G.q, G.phi0, r.phase, sprintf('%s, grain %d', frame_name, grain), 'animate', make_gif, files{:});
exportgraphics(fig, fullfile(result_dir, sprintf('orientation_steps_grain%d.png', grain)), 'Resolution', 150);
fprintf('Saved to %s\n', result_dir);

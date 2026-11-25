%% SHOW_VARIANTS_3D  Look at saved variants in 3D (no MTEX needed).
%
% Opens the 3D pictures of the variants saved by make_variants_mtex.m:
%   - overview: the c axes of all variants in the parent cell, and variant 1
%   - every variant with its parallel planes and directions
%   - a viewer: choose the variant in the menu, turn it with the mouse
% Parallel PLANES: parent grey, child coloured; parallel DIRECTIONS: parent
% black line, child coloured dashed line; their indices are in the titles.

clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here) || ~isfile(fullfile(here, 'show_variants_3d.m')), here = fileparts(which('show_variants_3d.m')); end
addpath(fullfile(here, 'src'));

%% Settings ===============================================================
variants_file = fullfile(here, 'variants', 'omega_in_beta.mat');    % or alpha2_in_beta.mat, alpha_hcp_in_beta.mat

%% Show ===================================================================
v = load(variants_file);
cifs = {v.parent_cif, v.child_cif};
for k = 1:2                                      % moved: look in cif/ next to the variants folder
    if ~isfile(cifs{k})
        [~, nm, ext] = fileparts(cifs{k});
        cifs{k} = fullfile(fileparts(fileparts(variants_file)), 'cif', [nm ext]);
    end
end
pc = read_cif(cifs{1});  cc = read_cif(cifs{2});
if isfield(v, 'parent_label'), pl = v.parent_label; else, pl = sprintf('%s (%s)', pc.name, pc.space_group); end
if isfield(v, 'child_label'),  cl = v.child_label;  else, cl = sprintf('%s (%s)', cc.name, cc.space_group); end
OR = [];
if isfield(v, 'OR_saved') && isstruct(v.OR_saved) && ~isfield(v.OR_saved, 'mtex_euler'), OR = v.OR_saved; end
fprintf('%s: %d variants of %s in %s\n', variants_file, size(v.V, 3), cl, pl);
[f1, f2, f3] = plot_variants_3d(v.V, pc, cc, OR, pl, cl);
figure(f3);                                      % the viewer in front

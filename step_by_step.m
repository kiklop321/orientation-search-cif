%% STEP_BY_STEP  The orientation search one step at a time, with a check after each step.
%
% The same calculation as run_orientation_search.m (orientation_pipeline.m),
% but split into sections that are run one after another with Ctrl+Enter
% ("Run Section"). Every section draws its intermediate result and prints
% what to look at, so that a wrong setting is seen in the step where it
% matters:
%
%   0. settings
%   1. read the frame      -> do the spots lie on the Debye rings of the phase?
%   2. find the spots      -> are all visible spots found, and no noise?
%      2b. one spot up close, 2c. number of spots vs. threshold
%   3. 2theta of the spots -> do they sit at allowed 2theta? (no orientation yet)
%   4. orientation search  -> is the best orientation clearly the best?
%   5. matches, residuals  -> are the residuals small and without trends?
%   6. lattice parameter   -> refined cell, residuals after
%   7. final picture
%   8. what the orientation means -> the cell turned by q and phi0, in the
%                                     beam and the detector (figures + animation)
%
% Run section 0 first; after changing a setting, run section 0 again and the
% sections from the one where the setting is used.

%% 0. Settings ============================================================
% (meaning of every setting: run_orientation_search.m, section 0)

clear; close all;
% folder of this script; with "Run Section" MATLAB runs a temporary copy, so
% mfilename does not give it - then it is found on the MATLAB path / current folder
here = fileparts(mfilename('fullpath'));
if isempty(here) || ~isfile(fullfile(here, 'step_by_step.m')), here = fileparts(which('step_by_step.m')); end
if isempty(here), error('Go to the folder of the tool (orientation-search-cif) first (cd), then run this section again.'); end
addpath(fullfile(here, 'src'), fullfile(here, 'external', 'RigakuRapidUtils'), fullfile(here, 'external', 'sphere'));

% --- what is analysed
image_file  = fullfile(here, 'data', 'beta_single_grain', 'Ti24Nb_1.img');  % the frame (.img)
cif_file    = fullfile(here, 'cif', 'beta_Ti_COD9012924.cif');                     % the phase: beta-Ti
phase       = phase_from_cif(cif_file, 'hkl_max', 6, 'cell', [3.30 3.30 3.30 90 90 90]);
%   'hkl_max', 6  reflections with every index from -6 to 6 are predicted
%   'cell', [...] starting lattice parameters [a b c alpha beta gamma] (A, deg),
%                 replace those of the CIF (Ti-24Nb: a ~ 3.30 A); refined in section 6

% --- instrument (sections 1-7)
S.wavelength     = 0.7093;          % X-ray wavelength (A); Mo K-alpha
beam_center      = [1274, 3247];    % detector pixel [row, column] where the direct beam
                                    % hits the plate (calibrated; the .img header value
                                    % [1287, 3246] is only nominal). The sample-detector
                                    % distance is read from the .img header

% --- rotation of the sample (sections 4-6)
S.rotation_range = [-30 30];        % how far the sample rotated during the exposure,
                                    % relative to the middle of the rotation (deg);
                                    % [-30 30] = a 60 deg oscillation
S.phi0_range     = [0 360];         % where the middle of the rotation, phi0, is searched
                                    % (deg); keep [0 360] (phi0 is not the goniometer angle)
S.phi0_step      = 0.25;            % step of phi0 in the search (deg)

% --- orientation search and matching (sections 2c, 3-6)
S.tth_range      = [0 180];         % only predicted spots with 2theta in this range are
                                    % used (deg); [0 180] = all
S.mesh_degree    = [7 10];          % how fine the trial orientations are: subdivisions of an
                                    % icosahedron, coarse stage and fine stage
                                    % (7 = about 7000 orientations for cubic, 10 = fine)
S.tolerance      = [0.25 0.3];      % largest |d2theta| and |dbeta| (deg) for a predicted and a
                                    % measured spot to count as the same spot, coarse and
                                    % fine stage; the 2nd value is also the "on a Debye
                                    % ring" limit of sections 2c and 3
S.unmatched_tol  = 2;               % (deg) a measured spot farther than this from every
                                    % prediction counts as unexplained (only reported)
S.search_region  = phase.parent_rotations;  % symmetry of the crystal: only one of the
                                    % symmetry-equivalent orientations is tried (24x fewer
                                    % for cubic); leave as it is

% --- spot search (section 2), on the frame compressed to log(log(I+1))
P.min_height           = 1.7;       % THRESHOLD: a spot must reach at least this value at its
                                    % maximum; 1.7 = 237 counts (printed in section 2).
                                    % Lower = more (weaker) spots, also more noise
P.min_prominence_map   = 0.09;      % how much a peak must rise above its neighbourhood in the
                                    % (2theta, beta) map (same units); higher = fewer spots
P.min_prominence_pixel = 0.01;      % the same for the search directly in detector pixels
P.min_distance_px      = 6;         % two peaks closer than this (px) are one spot (the
                                    % brighter one is kept)
P.refine_radius_px     = 3;         % a peak found in the map is moved to the brightest pixel
                                    % within this distance (px) on the detector
P.tth_beta_map         = [0.6, pi/2, 1.2, 2*pi];   % the (2theta, beta) map of the first search:
                                    % [2theta centre, beta centre, 2theta width, beta width]
                                    % in rad = 2theta 0-69 deg, all beta

% helpers
counts     = @(h) exp(exp(h)) - 1;                % log(log(I+1)) -> counts I
ring_index = @(t, rings) knnsearch(rings(:), t(:));          % nearest allowed 2theta
tth_of_g   = @(g) 2 * asind(vecnorm(g, 2, 2) * S.wavelength / (4 * pi));
fig_pos    = [60 60 1150 560];                % size of the figures (px), text added below

% --- section 8
play_animation = false;             % true = also play the movement of the cell (q, then phi0) in
                                    % a separate window; slower (about 10 s)

%% 1. Read the frame and check the geometry ===============================
% The Debye rings of all allowed reflections of the phase are drawn on the
% frame. Before any search, every real spot of this phase must lie ON a ring.
%   spots beside all rings, all to the same side   -> lattice parameter ('cell')
%   on the rings on one side of the frame, off on the other -> beam centre
%   many spots far from any ring                   -> other phase / other grain / noise

d = RAPIDread(char(image_file));
d.zero = beam_center;
geom = RAPIDutils.getDiffractionGeometry(d);
img_ll = log(log(double(d.image) + 1));
fprintf('Frame %d x %d px, distance %.2f mm (header), beam centre [%d %d]\n', ...
    size(d.image), geom.distance, beam_center);

[g, hkl, variant] = phase_reflections(phase, S.wavelength);
tth_all = tth_of_g(g);
[rings, ia] = uniquetol(tth_all, 1e-6);           % allowed 2theta, one per family
ring_hkl = sort(abs(hkl(ia, :)), 2, 'descend');
fprintf('%s: a b c = %.4f %.4f %.4f A, %d reflections, %d different 2theta\n', ...
    phase.name, phase.cell(1:3), size(g, 1), numel(rings));
fprintf('  first rings: %s\n', strjoin(compose('%d%d%d %.2f', [ring_hkl(1:min(8, end), :), rings(1:min(8, end))]), ',  '));

figure(1); clf; set(gcf, 'Color', 'w', 'Name', '1 frame and Debye rings', 'Position', fig_pos);
imagesc(img_ll); colormap(parula); axis image; hold on
h = RAPIDutils.drawImgDebyeRing(geom, rings(rings < 120) * pi / 180, gca);
set(h, 'Color', [1 1 1 0.5], 'LineWidth', 0.5);
plot(beam_center(2), beam_center(1), 'r+', 'MarkerSize', 15, 'LineWidth', 2);
for k = 1:min(8, numel(rings))                    % labels at the top of the first rings
    [x, z] = RAPIDutils.diff2ipPixels(geom, rings(k) * pi / 180, pi / 2);
    if x > 1 && x < size(img_ll, 2) && z > 1 && z < size(img_ll, 1)
        text(x, z, sprintf('%d%d%d', ring_hkl(k, :)), 'Color', 'w', 'FontSize', 8, 'VerticalAlignment', 'bottom');
    end
end
title('Debye rings of the phase (white), beam centre (red)');
xlabel('column (px)'); ylabel('row (px)');
describe_figure(gcf, { ...
    'WHAT: the measured frame, shown as log(log(I+1)) so that weak and strong spots are both visible (brighter = more counts). White circles: the 2theta of every reflection the phase allows (Debye rings), labelled hkl at the top. Red +: the direct-beam position (beam_center). The dark bar from the centre is the beam-stop holder. A single crystal gives spots only at a few places on each ring.', ...
    'HOW TO READ: every real spot of this phase must lie ON a white circle. Spots just inside or just outside the circles, all on the same side -> the lattice parameter is wrong (set ''cell'' in section 0). Spots on the circles on one side of the beam but off on the other side -> the beam centre is wrong. Spots far from all circles -> another phase, another grain or noise.'});

%% 2. Find the spots ======================================================
% Two searches are merged: (o) in the frame resampled to (2theta, beta) and
% (x) directly in detector pixels; (+) the spots kept.
% Zoom in (magnifier in the figure toolbar) and check:
%   visible spots without a +          -> lower P.min_height or the prominences
%   + on background, streaks, beam stop -> raise them
% The threshold in counts is printed; compare it with the counts of a weak
% spot in section 2b.

tic;
[expt, peaks_px, img_ll, info] = find_bragg_peaks(d, geom, P);
fprintf('Spots: %d from the map, %d from the pixels -> %d kept (%.1f s)\n', ...
    size(info.from_map, 1), size(info.from_pixels, 1), size(expt, 1), toc);
fprintf('Threshold P.min_height = %.2f  =  %.0f counts; found spots: %.0f to %.0f counts\n', ...
    P.min_height, counts(P.min_height), counts(min(peaks_px(:, 3))), counts(max(peaks_px(:, 3))));

figure(2); clf; set(gcf, 'Color', 'w', 'Name', '2 spot search', 'Position', fig_pos);
subplot(1, 3, 1:2);
imagesc(img_ll); colormap(parula); axis image; hold on
plot(info.from_map(:, 2), info.from_map(:, 1), 'co', 'MarkerSize', 9, 'DisplayName', 'search in (2\theta, \beta) map');
plot(info.from_pixels(:, 2), info.from_pixels(:, 1), 'mx', 'MarkerSize', 9, 'DisplayName', 'search in pixels');
plot(peaks_px(:, 2), peaks_px(:, 1), 'r+', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'kept');
legend('Location', 'southoutside', 'Orientation', 'horizontal');
title(sprintf('%d spots', size(expt, 1)));
subplot(1, 3, 3);
v = img_ll(isfinite(img_ll));
v = v(randperm(numel(v), min(2e5, numel(v))));
histogram(v, 100, 'Normalization', 'probability', 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'DisplayName', 'all pixels'); hold on
histogram(peaks_px(:, 3), 20, 'Normalization', 'probability', 'FaceColor', 'r', 'EdgeColor', 'none', 'DisplayName', 'spot maxima');
xline(P.min_height, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
set(gca, 'YScale', 'log'); xlabel('log(log(I+1))'); legend('Location', 'northwest');
title(sprintf('threshold (dashed) = %.0f counts', counts(P.min_height)));
describe_figure(gcf, { ...
    'LEFT: the frame with the result of the spot search. Cyan o = peak found in the frame resampled to (2theta, beta); magenta x = peak found directly in the detector pixels; red + = the spots kept (of two closer than P.min_distance_px only the brighter one). RIGHT: histogram of all pixel values (grey, mostly background) and of the maxima of the found spots (red), both as log(log(I+1)); dashed line = the threshold P.min_height (its value in counts in the title).', ...
    'HOW TO READ: zoom in (magnifier in the toolbar). Every spot you see by eye should have a red +; a + on plain background, on streaks or on the beam stop is a false spot. The red histogram should lie clearly to the right of the grey one. Spots visible but without + -> lower P.min_height (or the prominences); many + on noise, threshold cutting into the grey tail -> raise it. Sections 2b and 2c help to choose.'});

%% 2b. (optional) One spot up close ======================================
% Put the [row, column] of a spot you think is missed (read it from figure 2
% with the data tip); empty = the weakest spot found. Shows the counts
% around it, the threshold and what was found there.

spot_rc = [];                                   % [row column], e.g. [1500 3800]
half = 30;                                      % half-size of the window (px)

% Which spot:
%   spot_rc = []          -> NOT random: always the WEAKEST spot that was found
%                            (lowest peaks_px(:,3), i.e. lowest log(log(I+1)) at
%                            its maximum). It sits closest to the threshold, so it
%                            shows what "just found" looks like and what weaker
%                            structure next to it is left out. Same frame and
%                            settings -> always the same spot (Ti24Nb_1: row 1229,
%                            column 1837, 250 counts vs. threshold 237).
%   spot_rc = [row col]   -> the window is centred there. It does not have to be a
%                            found spot: any place on the frame works, e.g. a spot
%                            you see in figure 2 without a red + (read its position
%                            with Data Tips: X = column, Y = row).
if isempty(spot_rc), [~, k] = min(peaks_px(:, 3)); spot_rc = peaks_px(k, 1:2); end
r = max(spot_rc(1) - half, 1):min(spot_rc(1) + half, size(d.image, 1));
c = max(spot_rc(2) - half, 1):min(spot_rc(2) + half, size(d.image, 2));
win = double(d.image(r, c));
[~, im] = max(win(:));  [ir, ic] = ind2sub(size(win), im);
fprintf('Window around [%d %d]: max %.0f counts at [%d %d], median %.0f, threshold %.0f counts\n', ...
    spot_rc, max(win(:)), r(ir), c(ic), median(win(:)), counts(P.min_height));

figure(3); clf; set(gcf, 'Color', 'w', 'Name', '2b one spot', 'Position', fig_pos);
subplot(1, 3, 1);
imagesc(c, r, log(log(win + 1))); axis image; colorbar; hold on
contour(c, r, log(log(win + 1)), [P.min_height P.min_height], 'w', 'LineWidth', 1.2);
in = peaks_px(:, 1) >= r(1) & peaks_px(:, 1) <= r(end) & peaks_px(:, 2) >= c(1) & peaks_px(:, 2) <= c(end);
plot(peaks_px(in, 2), peaks_px(in, 1), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
title('white = threshold, + = found');
subplot(1, 3, 2);
plot(c, win(ir, :), 'b.-'); hold on; yline(counts(P.min_height), 'k--', 'threshold');
xlabel('column (px)'); ylabel('counts'); title(sprintf('row %d', r(ir)));
subplot(1, 3, 3);
plot(r, win(:, ic), 'b.-'); hold on; yline(counts(P.min_height), 'k--', 'threshold');
xlabel('row (px)'); ylabel('counts'); title(sprintf('column %d', c(ic)));
describe_figure(gcf, { ...
    sprintf(['WHICH SPOT: [row %d, column %d]. With spot_rc = [] it is not random but always the WEAKEST spot that was found (lowest intensity ' ...
    'at its maximum), the one closest to the threshold - so you see what "just found" looks like and what weaker structure next to it is ' ...
    'left out. To look at another place, e.g. a spot without a red + in figure 2, set spot_rc = [row column] (Data Tips in figure 2: ' ...
    'X = column, Y = row); it does not have to be a found spot.'], round(spot_rc)), ...
    'WHAT: a window of +-half px around that spot. LEFT: log(log(I+1)) of the window; white contour = the threshold P.min_height; red + = spots found in the window. MIDDLE and RIGHT: the counts along the row and along the column through the brightest pixel; dashed line = the threshold in counts.', ...
    'HOW TO READ: a spot is found only if its maximum is above the dashed line and it is a local maximum along both its row and its column (with at least the set prominence). A spot that you see but whose peak stays below the dashed line needs a lower P.min_height. Put the [row column] of a missed spot (data tip in figure 2) into spot_rc and run this section again.'});

%% 2c. (optional) Number of spots vs. threshold ===========================
% Repeats the search for several P.min_height (a few seconds each) and counts
% how many spots lie ON A DEBYE RING of the phase, i.e. with 2theta within
% +-S.tolerance(2) of an allowed 2theta (white circles of figure 1). A good
% threshold keeps most of the on-ring spots while the off-ring ones (not this
% phase: noise, artefacts, other phase) stay few.

heights = 1.5:0.1:2.1;
t = linspace(min(expt(:, 1)), max(expt(:, 1)), 20000)';      % fraction of the 2theta range
chance_frac = mean(abs(t - rings(ring_index(t, rings))) <= S.tolerance(2));   % covered by the rings
n_all = zeros(size(heights));  n_on = n_all;
for k = 1:numel(heights)
    P2 = P;  P2.min_height = heights(k);
    e2 = find_bragg_peaks(d, geom, P2);
    n_all(k) = size(e2, 1);
    n_on(k) = sum(abs(e2(:, 1) - rings(ring_index(e2(:, 1), rings))) <= S.tolerance(2));
    fprintf('  min_height %.1f (%5.0f counts): %4d spots, %4d on a Debye ring\n', heights(k), counts(heights(k)), n_all(k), n_on(k));
end
figure(4); clf; set(gcf, 'Color', 'w', 'Name', '2c threshold', 'Position', fig_pos);
plot(heights, n_all, 'ko-', heights, n_on, 'bo-', heights, n_all - n_on, 'ro-', 'LineWidth', 1.5);
xline(P.min_height, 'k--', 'current');
legend('all spots', sprintf('on a Debye ring of the phase (within \\pm%.2f%s in 2\\theta)', S.tolerance(2), char(176)), ...
    'on no Debye ring (not this phase: noise, artefact, other phase)');
xlabel('P.min\_height  (log(log(I+1)))'); ylabel('number of spots'); grid on
describe_figure(gcf, { ...
    'WHAT: how many spots the search finds for different thresholds P.min_height (bottom axis, in log(log(I+1)); the counts are printed in the command window). Black = all spots. Blue = spots ON A DEBYE RING of the phase (the white circles of figure 1), i.e. their 2theta is within +-S.tolerance(2) of a 2theta the phase allows - consistent with the phase, most likely real reflections. Red = spots on no Debye ring - not this phase: noise, artefacts or another phase. Blue + red = black. Dashed line = the current setting.', ...
    sprintf(['HOW STRICT "ON A RING" IS: +-%.2f deg in 2theta, about +-%.0f px on the detector; real spots lie much closer (typically 0.03 deg). ' ...
    'All the bands together cover about %.0f %% of the 2theta range, so a spot unrelated to the phase still lands on a ring by chance in about ' ...
    '%.0f %% of cases: blue is "consistent with the phase", not proof for a single spot (that comes with the orientation, section 4).'], ...
    S.tolerance(2), S.tolerance(2) * pi / 180 * geom.distance / geom.pixSize(1), 100 * chance_frac, 100 * chance_frac), ...
    'HOW TO READ: going to a lower threshold (to the left) should add blue spots. Where it adds mainly red ones, it only adds noise. A good threshold is where the red curve is still low while the blue one is as high as possible.'});

%% 3. Are the spots at allowed 2theta? ====================================
% 2theta depends only on the cell, not on the orientation, so this checks
% phase, cell and geometry BEFORE the search (like a powder pattern).
%   many spots on rings (far more than by chance)  -> phase OK, go on
%   offsets growing with 2theta, same sign         -> cell: give 'cell' in section 0
%   offsets going up and down with beta (sine)     -> beam centre
%   few spots on rings                             -> wrong phase, or noise (section 2)

off = expt(:, 1) - rings(ring_index(expt(:, 1), rings));
on = abs(off) <= S.tolerance(2);
t = linspace(min(expt(:, 1)), max(expt(:, 1)), 20000)';
chance = mean(abs(t - rings(ring_index(t, rings))) <= S.tolerance(2));
fprintf('%d of %d spots (%.0f %%) within %.2f deg of an allowed 2theta; by chance %.0f %%\n', ...
    sum(on), numel(on), 100 * mean(on), S.tolerance(2), 100 * chance);
fprintf('Offset of the on-ring spots: median %+.3f deg (systematic -> cell / distance)\n', median(off(on)));

figure(5); clf; set(gcf, 'Color', 'w', 'Name', '3 2theta of the spots', 'Position', fig_pos);
subplot(2, 1, 1);
yl = [0 max(histcounts(expt(:, 1), 0:0.2:max(expt(:, 1)) + 1)) + 1];
for k = find(rings <= max(expt(:, 1)) + 1)'
    plot([rings(k) rings(k)], yl, '-', 'Color', [1 0.6 0.6]); hold on
    text(rings(k), yl(2), sprintf('%d%d%d', ring_hkl(k, :)), 'Rotation', 90, 'FontSize', 7, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Color', 'r');
end
histogram(expt(:, 1), 0:0.2:max(expt(:, 1)) + 1, 'FaceColor', [0.1 0.1 0.8], 'EdgeColor', [0.1 0.1 0.8]);
ylim(yl);
xlabel('2\theta (deg)'); ylabel('spots'); title('2\theta of the spots (blue) and allowed reflections (red)');
subplot(2, 2, 3);
plot(expt(~on, 1), off(~on), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 10); hold on
plot(expt(on, 1), off(on), 'b.', 'MarkerSize', 10);
yline(0, 'k-'); ylim([-1 1]); xlabel('2\theta (deg)'); ylabel('offset from ring (deg)');
title('trend with 2\theta -> cell');
subplot(2, 2, 4);
plot(expt(on, 2), off(on), 'b.', 'MarkerSize', 10); yline(0, 'k-');
ylim(S.tolerance(2) * [-1.5 1.5]); xlabel('\beta (deg)'); ylabel('offset from ring (deg)');
title('sine in \beta -> beam centre');
describe_figure(gcf, { ...
    'TOP: histogram of the 2theta of all found spots (blue bars) and the 2theta of all allowed reflections of the phase (red lines, hkl). It is like a powder pattern made from the spots and does not need the orientation. BOTTOM: the distance of each spot from the nearest allowed 2theta, against 2theta (left) and against beta, the angle around the beam (right); blue = within the tolerance, grey = further away.', ...
    'HOW TO READ: the blue bars should stand on red lines (the fraction on the lines and the fraction expected by chance are printed). Distances growing with 2theta, all of one sign -> the lattice parameter is off. Distances going up and down with beta (like a sine) -> the beam centre is off. Distances scattered around zero -> phase and geometry are fine; go on with the search.'});

%% 4. Orientation search ==================================================
% Two stages (coarse grid in the fundamental region, then fine around the
% best). Printed: the five best orientations of stage 1 - the best should
% be clearly above the rest; if several are about equal, the solution is
% not unique (too few spots, or several grains: S.n_grains in run_orientation_search).

tic;
result = search_orientation(expt, g, S);
fprintf('Best: q = [%.5f %.5f %.5f %.5f], phi0 = %.2f deg, %d spots matched (%.0f s)\n', ...
    result.q, result.phi0, result.score, toc);
sc = result.stage(1).scores;
fprintf('Stage 1, best five scores (matched spots): %s\n', mat2str(sc(1:min(5, end), 2)'));

[T, calc] = match_spots(result.q, result.phi0, g, hkl, variant, expt, S);
[x, z] = RAPIDutils.diff2ipPixels(geom, calc(:, 1) * pi / 180, calc(:, 2) * pi / 180);
on_plate = x >= 1 & x <= size(img_ll, 2) & z >= 1 & z <= size(img_ll, 1);

% the frame resampled to (2theta, beta) as background; NaN = not on the plate
tth_map = 0:0.1:max(expt(:, 1)) + 5;
beta_map = -100:0.15:300;
[Xt, Xb] = meshgrid(tth_map, beta_map);
[xm, zm] = RAPIDutils.diff2ipPixels(geom, Xt * pi / 180, Xb * pi / 180);
xm = round(xm);  zm = round(zm);
ok = xm >= 1 & xm <= size(img_ll, 2) & zm >= 1 & zm <= size(img_ll, 1);
frame_map = nan(size(Xt));
frame_map(ok) = img_ll(sub2ind(size(img_ll), zm(ok), xm(ok)));

figure(6); clf; set(gcf, 'Color', 'w', 'Name', '4 orientation', 'Position', fig_pos);
imagesc(tth_map, beta_map, frame_map, 'AlphaData', ~isnan(frame_map));
set(gca, 'YDir', 'normal', 'Color', 'w');                      % off the plate: white
colormap(gca, parula);                                          % as in figures 1 and 2
clim([min(frame_map(ok)) max(frame_map(ok))]);
hold on
plot(calc(on_plate, 1), calc(on_plate, 2), 'o', 'Color', [0 0.6 0], 'MarkerSize', 6, 'DisplayName', 'predicted (on the plate)');
plot(expt(:, 1), expt(:, 2), 'r+', 'MarkerSize', 6, 'LineWidth', 1.2, 'DisplayName', 'measured');
plot(T.tth_calc, T.beta_calc, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'predicted and matched');
xlabel('2\theta (deg)'); ylabel('\beta (deg)'); legend('Location', 'eastoutside'); grid on
xlim([0 max(expt(:, 1)) + 5]);
title(sprintf('best orientation: %d of %d spots matched', numel(unique(T.spot)), size(expt, 1)));
describe_figure(gcf, { ...
    'WHAT: the result of the orientation search in the (2theta, beta) plane; beta is the angle around the beam (the azimuth on the detector). Background: the frame itself, recalculated to (2theta, beta), in the same colours as figure 1 (log(log(I+1)): blue = background, yellow = intense); the Debye rings are vertical lines here. White = directions that do not hit the image plate (the beam hits the cylindrical plate off its middle: on one side (beta ~ 0) it reaches 2theta ~ 146 deg, on the other (beta ~ 180) only ~ 61 deg, up and down about 45 deg), so no spot can be there. Red + = measured spots. Green o = spots predicted by the best orientation (only those that hit the plate). Black o = predicted spots with a measured spot within the tolerance (matched).', ...
    'HOW TO READ: most red + should carry a black circle. Green circles without a + are reflections predicted but not seen - usually weak ones at high 2theta, which is normal. Red + without any circle are spots this orientation does not explain (another grain, noise). If only few spots match, check section 3 first; if the best scores printed in the command window are all about equal, the solution is not unique.'});

%% 5. Matches and residuals ===============================================
% Every matched pair: predicted minus measured 2theta and beta.
%   d2theta trend with 2theta           -> cell (section 6 fixes it)
%   d2theta going up and down with beta -> beam centre
%   dbeta large at low 2theta only      -> normal (beta is ill-defined near the beam)
%   residuals close to the tolerance    -> tolerance too tight / wrong orientation

matched = unique(T.spot);
fprintf('%d predicted reflections matched to %d of %d spots; unexplained %d\n', ...
    height(T), numel(matched), size(expt, 1), size(expt, 1) - numel(matched));
fprintf('d2theta %+.3f +- %.3f deg, dbeta %+.3f +- %.3f deg\n', mean(T.d_tth), std(T.d_tth), mean(T.d_beta), std(T.d_beta));
disp(head(T(:, {'h', 'k', 'l', 'tth_calc', 'd_tth', 'beta_calc', 'd_beta', 'phi_minus_phi0'}), 12));
rest = setdiff(1:size(expt, 1), matched);
fprintf('Unexplained spots: %d on an allowed 2theta (other grain?), %d off (noise / other phase)\n', ...
    sum(on(rest)), sum(~on(rest)));

figure(7); clf; set(gcf, 'Color', 'w', 'Name', '5 residuals', 'Position', fig_pos);
subplot(2, 2, 1); plot(T.tth_calc, T.d_tth, 'b.', 'MarkerSize', 12); yline(0); grid on
xlabel('2\theta (deg)'); ylabel('\Delta2\theta (deg)'); title('trend -> cell');
subplot(2, 2, 2); plot(T.beta_calc, T.d_tth, 'b.', 'MarkerSize', 12); yline(0); grid on
xlabel('\beta (deg)'); ylabel('\Delta2\theta (deg)'); title('sine -> beam centre');
subplot(2, 2, 3); plot(T.tth_calc, T.d_beta, 'b.', 'MarkerSize', 12); yline(0); grid on
xlabel('2\theta (deg)'); ylabel('\Delta\beta (deg)');
subplot(2, 2, 4); histogram(max(abs(T.d_tth), abs(T.d_beta)), 0:0.02:S.tolerance(2) + 0.02); grid on
xline(S.tolerance(2), 'k--', 'tolerance'); xlabel('max(|\Delta2\theta|, |\Delta\beta|) (deg)'); ylabel('pairs');
describe_figure(gcf, { ...
    'WHAT: the residuals (predicted minus measured) of every matched spot. Top left: d2theta against 2theta. Top right: d2theta against beta. Bottom left: dbeta against 2theta. Bottom right: histogram of the larger of |d2theta| and |dbeta|, with the matching tolerance (dashed).', ...
    'HOW TO READ: the points should scatter around zero without a pattern. d2theta growing or falling with 2theta -> the lattice parameter (section 6 refines it). d2theta following a sine in beta -> the beam centre (S.refine_beam_center in run_orientation_search). dbeta is naturally larger at low 2theta, where beta is poorly defined. Pairs piling up close to the tolerance -> tolerance too tight, or the orientation is not right.'});

%% 6. Refine the lattice parameter ========================================
% The cell is fitted to the measured 2theta of the matched spots, the
% reflections recomputed and the orientation refined (fine stage only);
% repeated until the cell does not change. Compare figure 8 with figure 7.

T0 = T;
if height(T) < 6
    error('Only %d matched spots - too few to refine the lattice parameter (see sections 2-5).', height(T));
end
for it = 1:10
    [cell_new, ci] = refine_cell(phase.cell, [T.h T.k T.l], T.tth_meas, S.wavelength);
    change = max(abs(cell_new(1:3) - phase.cell(1:3)));
    phase.cell = cell_new;
    [g, hkl, variant] = phase_reflections(phase, S.wavelength);
    S2 = S;  S2.start_q = result.q;
    result = search_orientation(expt, g, S2);
    [T, calc] = match_spots(result.q, result.phi0, g, hkl, variant, expt, S);
    fprintf('round %d: a b c = %.4f %.4f %.4f A (esd %.4f), rms d2theta %.3f -> %.3f deg, %d spots matched\n', ...
        it, cell_new(1:3), ci.esd(1), ci.rms_before, ci.rms_after, numel(unique(T.spot)));
    if change < 1e-5 * max(cell_new(1:3)), break; end
end
rings = uniquetol(tth_of_g(g), 1e-6);              % rings of the refined cell for the next sections

figure(8); clf; set(gcf, 'Color', 'w', 'Name', '6 after refinement', 'Position', fig_pos);
plot(T0.tth_calc, T0.d_tth, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 12); hold on
plot(T.tth_calc, T.d_tth, 'b.', 'MarkerSize', 12); yline(0); grid on
legend(sprintf('before (%d spots)', numel(unique(T0.spot))), sprintf('after (%d spots)', numel(unique(T.spot))));
xlabel('2\theta (deg)'); ylabel('\Delta2\theta (deg)');
title(sprintf('a b c = %.4f %.4f %.4f A', phase.cell(1:3)));
describe_figure(gcf, { ...
    'WHAT: d2theta (predicted minus measured 2theta) of the matched spots against 2theta, before (grey) and after (blue) the lattice parameter is refined to these spots; the refined cell is in the title, the number of matched spots in the legend.', ...
    'HOW TO READ: after the refinement the blue points should scatter around zero without a trend, and at least as many spots should be matched as before. A trend that stays -> look at the beam centre (sine in beta in figure 7).'});

%% 7. Final picture =======================================================
% The same figure as orientation_overlay.png of run_orientation_search.

grains = struct('q', result.q, 'phi0', result.phi0, 'T', T, 'calc', calc, 'n_spots', numel(unique(T.spot)));
fig = plot_overlay(img_ll, geom, peaks_px, grains, beam_center, ...
    sprintf('%s: %d of %d spots, a = %.4f A', phase.name, grains.n_spots, size(expt, 1), phase.cell(1)));
describe_figure(fig, { ...
    'WHAT: the final result, the same picture as orientation_overlay.png of run_orientation_search. Red + = measured spots; black o = reflections of the found orientation with a measured spot within the tolerance; green o = reflections predicted but without a spot; cyan + = beam centre.', ...
    'HOW TO READ: a good solution has black circles on most red crosses all over the frame; green circles are mostly weak reflections at high 2theta.'});
% exportgraphics(fig, fullfile(here, 'results', 'step_by_step_overlay.png'), 'Resolution', 200);

%% 8. What the orientation means: the unit cell, the beam and the detector =
% The orientation found is described by a crystal direction n (a unit
% vector, 2 parameters) that the quaternion q tilts onto the rotation axis z,
% and a turn phi0 about z (1 parameter). Two figures show it:
%   - the six stages: start, tilt by q (n goes to z), turn by phi0, final;
%     and an animation of the same movement in a separate window
%   - the cell inside the cylindrical image plate with the beam, and the
%     diffracted beams of the matched reflections ending at their spots
% The axis and angle of q, the crystal directions along x, y, z and the
% Euler angles (Bunge) are printed in the command window.
% To save the animation, give file names, e.g.
%   plot_orientation_steps(result.q, result.phi0, phase, name, 'animate', true, ...
%       'gif', fullfile(here, 'results', 'orientation.gif'), 'mp4', fullfile(here, 'results', 'orientation.mp4'));

[~, name] = fileparts(image_file);
plot_orientation_steps(result.q, result.phi0, phase, name, 'animate', play_animation);
plot_orientation_3d(result.q, result.phi0, phase, geom, T, S, name, 12);

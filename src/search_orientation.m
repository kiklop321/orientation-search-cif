function result = search_orientation(expt, g, S)
%SEARCH_ORIENTATION  Orientation of the parent crystal that reproduces the spots.
%
%   result = search_orientation(expt, g, S)
%
%   expt : N x 2 measured spots [2theta, beta] (deg), from FIND_BRAGG_PEAKS
%   g    : M x 3 reciprocal vectors of all variants in the parent frame
%          (from PHASE_REFLECTIONS)
%   S    : search settings (see run_orientation_search.m)
%            S.wavelength       (A)
%            S.rotation_range   [from to] rotation of the sample during the
%                               exposure, relative to its centre phi0 (deg)
%            S.phi0_range       [from to] where the centre phi0 is searched (deg)
%            S.phi0_step        step of phi0 (deg)
%            S.tth_range        [min max] 2theta used for matching (deg)
%            S.mesh_degree      [stage1 stage2] subdivisions of the icosahedron
%            S.tolerance        [stage1 stage2] max |d2theta| and |dbeta| (deg)
%            S.unmatched_tol    tolerance for counting unexplained spots (deg)
%            S.search_region    'octant' or 3x3xN proper rotations of the crystal
%                               (phase.parent_rotations): directions tested
%                               only in the fundamental region (ORIENTATION_GRID)
%
%   Method
%   An orientation of the parent crystal = a direction n on the unit
%   sphere (two angles) + the rotation phi0 about the sample rotation axis.
%   Stage 1 tests the directions of a subdivided icosahedron that lie in
%   the fundamental region of the crystal symmetry (S.search_region; e.g.
%   6 998 directions for a cubic crystal and mesh degree 7), stage 2 the
%   4^(d2-d1+2) directions of a finer mesh closest to the stage-1 winner.
%   For every direction the reciprocal vectors are rotated, their
%   diffraction positions (2theta, beta) and rotation angles phi computed
%   (EWALD_ANGLES) and compared with the measured spots. For every phi0 on
%   the grid, the number of DIFFERENT measured spots matched by predictions
%   with phi - phi0 inside S.rotation_range is counted; the best count is
%   the score. If several phi0 give the best count, phi0 is the centre of
%   the largest contiguous (circular) range of them. Orientations whose
%   total number of matches is lower than the best score so far cannot win
%   and are skipped.
%
%   result.q          best orientation (quaternion [w x y z])
%   result.phi0       centre of the rotation range of the best orientation (deg)
%   result.score      number of matched measured spots
%   result.stage(s)   .q (all tested orientations, best first) and .scores:
%                     [n_matches_any_phi, score, phi0_first, phi0_last,
%                      n_best_phi0, (stage 2: trimmed mean squared beta
%                      residual,) unexplained low-angle spots]
%                     phi0_first..phi0_last is the contiguous range of the
%                     best phi0 (it may wrap through 0/360 deg)

phi0_grid = S.phi0_range(1):S.phi0_step:S.phi0_range(2);
phi0_grid = phi0_grid(phi0_grid < S.phi0_range(1) + 360);       % 0 and 360 are the same

if isfield(S, 'start_q') && ~isempty(S.start_q)
    % only stage 2 around a known orientation (e.g. after refining the cell)
    result.stage(1).q = S.start_q;
    result.stage(1).scores = [];
else
    q1 = orientation_grid(S.mesh_degree(1), S.search_region);
    [result.stage(1).scores, result.stage(1).q] = run_stage(q1, g, expt, phi0_grid, S.tolerance(1), false, S);
end

n2 = 4^(S.mesh_degree(2) + 2 - S.mesh_degree(1));
q2 = closest_quaternions(result.stage(1).q(1,:), S.mesh_degree(2), S.search_region);
q2 = q2(1:n2, :);
[result.stage(2).scores, result.stage(2).q] = run_stage(q2, g, expt, phi0_grid, S.tolerance(2), true, S);

result.q     = result.stage(2).q(1,:);
first = result.stage(2).scores(1,3);  last = result.stage(2).scores(1,4);
result.phi0  = mod(first + mod(last - first, 360) / 2, 360);   % centre, also across 0/360
result.score = result.stage(2).scores(1,2);
end

%% ------------------------------------------------------------------------
function [scores, Q] = run_stage(Q, g, expt, phi0_grid, tol, with_loss, S)
fprintf('  %d orientations\n', size(Q, 1));
scores = zeros(size(Q, 1), 6 + with_loss);
best_so_far = 0;
n_check = round(size(expt, 1) / 3);
width = S.rotation_range(2) - S.rotation_range(1);

for i = 1:size(Q, 1)
    [tth, b1, b2, p1, p2] = ewald_angles(quatrotate(quatconj(Q(i,:)), g), S.wavelength);
    calc = [tth, b1, p1; tth, b2, p2];
    calc = calc(calc(:,1) >= S.tth_range(1) & calc(:,1) <= S.tth_range(2), :);
    calc = wrap_beta(calc);

    [hit, idx] = ismembertol(calc(:,1:2), expt, tol, 'ByRows', true, 'DataScale', [1 1]);
    explained = ismembertol(expt(1:n_check,:), calc(:,1:2), S.unmatched_tol, 'ByRows', true, 'DataScale', [1 1]);
    if sum(hit) < best_so_far
        continue
    end

    % matched predictions inside the rotation range, for every phi0
    in_range = mod(calc(hit,3) - phi0_grid - S.rotation_range(1), 360) <= width;
    matched = idx(hit) .* in_range;
    s = sort(matched, 1);
    n_distinct = sum(diff([zeros(1, numel(phi0_grid)); s], 1, 1) ~= 0 & s ~= 0, 1);

    score = max(n_distinct);
    best_so_far = max(best_so_far, score);
    [first, last, nbest] = best_phi0_range(phi0_grid, n_distinct == score, S);
    row = [sum(hit), score, first, last, nbest];
    if with_loss
        % beta residual at the phi0 that is returned (centre of the best run)
        centre = mod(first + mod(last - first, 360) / 2, 360);
        [~, I] = min(abs(mod(phi0_grid - centre + 180, 360) - 180));
        cm = calc(hit, :);
        em = expt(idx(hit), :);
        r2 = sort((mod(cm(in_range(:,I), 2) - em(in_range(:,I), 2) + 180, 360) - 180).^2);
        row = [row, mean(r2(1:round(3 * numel(r2) / 4)))]; %#ok<AGROW>
    end
    scores(i,:) = [row, sum(~explained)];
    if mod(i, 5000) == 0, fprintf('    %d / %d\n', i, size(Q, 1)); end
end

if with_loss
    % stage 2: among orientations with the same score the smallest beta residual wins
    [~, order] = sortrows(scores, [-2 6]);
else
    [~, order] = sort(scores(:,2), 'descend');
end
scores = scores(order, :);
Q = Q(order, :);
end

function calc = wrap_beta(calc)
% beta runs from -90 to 270 deg (EWALD_ANGLES, RAPIDutils): -90 and 270 are the
% same direction (straight below the beam). Predictions within 5 deg of this
% seam are added once more with beta -+ 360, so that a spot on the other
% side of the seam can still be matched.
calc = [calc; calc(calc(:,2) > 265, :) - [0 360 0]; calc(calc(:,2) < -85, :) + [0 360 0]];
end

function [first, last, n] = best_phi0_range(grid, is_best, S)
% Largest contiguous run of best phi0 values; runs at both ends of the grid
% are joined when the grid covers the full circle.
idx = find(is_best);
n = numel(idx);
breaks = [0, find(diff(idx) > 1), numel(idx)];
runs = arrayfun(@(k) idx(breaks(k)+1:breaks(k+1)), 1:numel(breaks)-1, 'UniformOutput', false);
full_circle = abs(S.phi0_range(2) - S.phi0_range(1) - 360) < 1e-9;
if full_circle && numel(runs) > 1 && runs{1}(1) == 1 && runs{end}(end) == numel(grid)
    runs{1} = [runs{end}, runs{1}];          % join the run that wraps through 360/0
    runs(end) = [];
end
[~, longest] = max(cellfun(@numel, runs));
r = runs{longest};
first = grid(r(1));
last = grid(r(end));
end

function out = orientation_pipeline(image_file, beam_center, phase, S, P, out_dir, name, plot_max_tth)
%ORIENTATION_PIPELINE  Spots, orientation search and output for one frame and one phase.
%
%   out = orientation_pipeline(image_file, beam_center, phase, S, P, out_dir, name, plot_max_tth)
%
%   Steps 1-6 of RUN_ORIENTATION_SEARCH (see there for the settings S, P):
%   read the frame, find the spots, build the reflections of the phase,
%   search the orientation, refine the lattice parameters, optionally
%   search further grains in the spots that are still unexplained, plot and
%   save to out_dir.
%
%   out.result, out.T      orientation and matched spots of grain 1
%   out.grains             struct array, one element per grain found:
%                          q, phi0, T (matched spots, spot = row of expt),
%                          n_spots, misorientation_to_1 (deg)
%   out.expt               all measured spots [2theta, beta]
%   out.cell, out.cell_log refined cell and the refinement history
%   out.beam_center, out.bc_log  beam centre used [row, column] and the
%                          history of its refinement [it, row, col, esd_row,
%                          esd_col, rms_before, rms_after]

%% 1. Read the frame ======================================================
if ~isfile(image_file)
    error('orientation:noFrame', ['Frame not found:\n  %s\n' ...
        'The example frames are not part of the repository: put them into data/ as listed in ' ...
        'data/README.md, or set image_file (section 0) to your own .img frame.'], image_file);
end
if ~isfolder(out_dir), mkdir(out_dir); end
d = RAPIDread(image_file);
d.zero = beam_center;
geom = RAPIDutils.getDiffractionGeometry(d);

%% 2. Find the spots ======================================================
[expt, peaks_px, img_ll] = find_bragg_peaks(d, geom, P);
fprintf('%s: %d spots found\n', name, size(expt, 1));

%% 3. Reflections of the phase ============================================
[g, hkl, variant] = phase_reflections(phase, S.wavelength);
fprintf('%s: %d reflections in %d variants\n', phase.name, size(g, 1) / size(phase.variants, 3), size(phase.variants, 3));

%% 4. Orientation search (grain 1) ========================================
tic;
result = search_orientation(expt, g, S);
fprintf('Best orientation q = [%.5f %.5f %.5f %.5f], phi0 = %.2f deg, %d spots matched (%.0f s)\n', ...
    result.q, result.phi0, result.score, toc);
[T, calc] = match_spots(result.q, result.phi0, g, hkl, variant, expt, S);
fprintf('%d predicted reflections matched to %d of %d measured spots, ', height(T), numel(unique(T.spot)), size(expt, 1));
fprintf('mean d2theta = %.3f deg, mean dbeta = %.3f deg\n', mean(T.d_tth), mean(T.d_beta));

%% 5. Refine the lattice parameters, the beam centre and the orientation =
% 2theta depends only on the cell, so the cell is fitted to the measured
% 2theta of the matched spots (REFINE_CELL). Optionally the direct-beam
% position is fitted to the same spots first (REFINE_BEAM_CENTER, cell and
% orientation kept; in the first round after a coarse scan of
% +-S.beam_center_scan px for the most matched spots) and the measured
% spots are recomputed with it. Then the
% reflections are recomputed and the orientation refined around the
% previous one (stage 2 only). Repeated until the cell and the beam centre
% no longer change, at most S.refine_cell times (10 if only the beam centre
% is refined).
cell_cif = phase.cell;
beam_center_start = beam_center;
cell_log = zeros(0, 9);
bc_log = zeros(0, 7);
if ~isfield(S, 'refine_cell'), S.refine_cell = 0; end
if ~isfield(S, 'refine_beam_center'), S.refine_beam_center = false; end
if ~isfield(S, 'beam_center_scan'), S.beam_center_scan = 20; end
n_it = S.refine_cell;
if S.refine_beam_center && n_it == 0, n_it = 10; end
for it = 1:n_it
    if height(T) < 6
        warning('orientation_pipeline:fewSpots', ['Only %d matched spots - the lattice ' ...
            'parameters and the beam centre are not refined.'], height(T));
        break
    end
    bc_old = beam_center;
    if S.refine_beam_center && it == 1 && S.beam_center_scan > 0
        % coarse scan first: the centre with the most matched spots (the fit
        % alone can stay in a local minimum if the start is several px off)
        [beam_center, n_best] = scan_beam_center(beam_center, d, peaks_px, result, g, hkl, variant, S);
        fprintf('beam centre scan (+-%g px): start [%.0f %.0f], %d spots matched\n', ...
            S.beam_center_scan, beam_center, n_best);
        % spots and their partners at the scanned centre, so that the fit
        % starts from these matches and not from those of the old centre
        d.zero = beam_center;
        geom = RAPIDutils.getDiffractionGeometry(d);
        [tth, beta] = RAPIDutils.ipPixels2diff(geom, peaks_px(:,2), peaks_px(:,1));
        expt = [tth, beta] * 180 / pi;
        [~, o] = sort(expt(:,1));                    % keep expt sorted by 2theta,
        expt = expt(o, :);  peaks_px = peaks_px(o, :);   % peaks_px in the same order
        [T, calc] = match_spots(result.q, result.phi0, g, hkl, variant, expt, S);
    end
    if S.refine_beam_center
        [beam_center, bi] = refine_beam_center(beam_center, geom, peaks_px, T, g, hkl, variant, ...
            result.q, S.wavelength);
        d.zero = beam_center;
        geom = RAPIDutils.getDiffractionGeometry(d);
        [tth, beta] = RAPIDutils.ipPixels2diff(geom, peaks_px(:,2), peaks_px(:,1));
        expt = [tth, beta] * 180 / pi;
        [~, o] = sort(expt(:,1));                    % keep expt sorted by 2theta,
        expt = expt(o, :);  peaks_px = peaks_px(o, :);   % peaks_px in the same order
        [T, calc] = match_spots(result.q, result.phi0, g, hkl, variant, expt, S);
        bc_log(end+1, :) = [it, beam_center, bi.esd, bi.rms_before, bi.rms_after]; %#ok<AGROW>
        fprintf('beam centre %d: [row column] = [%.2f %.2f] px (esd %.2f %.2f), rms residual %.3f -> %.3f deg\n', ...
            it, beam_center, bi.esd, bi.rms_before, bi.rms_after);
    end
    if S.refine_cell > 0
        [cell_new, ci] = refine_cell(phase.cell, [T.h T.k T.l], T.tth_meas, S.wavelength);
        phase.cell = cell_new;
        cell_log(end+1, :) = [it, cell_new(1:3), ci.esd, ci.rms_before, ci.rms_after]; %#ok<AGROW>
        [g, hkl, variant] = phase_reflections(phase, S.wavelength);
    end
    S2 = S;  S2.start_q = result.q;
    result = search_orientation(expt, g, S2);
    [T, calc] = match_spots(result.q, result.phi0, g, hkl, variant, expt, S);
    if S.refine_cell > 0
        fprintf('cell refinement %d: a b c = %.4f %.4f %.4f A (esd %.4f %.4f %.4f), rms d2theta %.3f -> %.3f deg; %d spots matched, mean d2theta = %.3f deg\n', ...
            it, cell_new(1:3), ci.esd, ci.rms_before, ci.rms_after, numel(unique(T.spot)), mean(T.d_tth));
    end
    cell_done = S.refine_cell == 0 || (size(cell_log, 1) > 1 && ...
        max(abs(cell_log(end, 2:4) - cell_log(end-1, 2:4))) < 1e-5 * max(phase.cell(1:3)));
    bc_done = ~S.refine_beam_center || (it > 1 && max(abs(beam_center - bc_old)) < 0.01);
    if it > 1 && cell_done && bc_done
        break                                   % converged
    end
    % the cell can alternate between two values when two orientations explain
    % the same spots (the difference is within the esd): stop then as well
    if S.refine_cell > 0 && size(cell_log, 1) > 2 && bc_done && ...
            max(abs(cell_log(end, 2:4) - cell_log(end-2, 2:4))) < 1e-5 * max(phase.cell(1:3))
        fprintf('cell refinement: alternates between two values (difference %.4f A, within the esd) - stopped\n', ...
            max(abs(cell_log(end, 2:4) - cell_log(end-1, 2:4))));
        break
    end
end
T.grain = ones(height(T), 1);
grains = struct('q', result.q, 'phi0', result.phi0, 'T', T, 'calc', calc, ...
                'n_spots', numel(unique(T.spot)), 'misorientation_to_1', 0);

%% 6. Further grains (optional) ===========================================
% The spots explained by the grains found so far are removed and the
% orientation search is repeated on the rest (with the refined cell). A new
% grain is accepted if it explains at least S.min_spots_per_grain spots;
% otherwise, or after S.n_grains grains, the search stops.
if ~isfield(S, 'n_grains'), S.n_grains = 1; end
if ~isfield(S, 'min_spots_per_grain'), S.min_spots_per_grain = 15; end
used = unique(T.spot);
for k = 2:S.n_grains
    rest = setdiff(1:size(expt, 1), used);
    if numel(rest) < S.min_spots_per_grain, break; end
    Sk = S;  Sk.start_q = [];
    fprintf('Grain %d: search in the %d unexplained spots\n', k, numel(rest));
    rk = search_orientation(expt(rest, :), g, Sk);
    [Tk, ck] = match_spots(rk.q, rk.phi0, g, hkl, variant, expt(rest, :), S);
    nk = numel(unique(Tk.spot));
    if nk < S.min_spots_per_grain
        fprintf('Grain %d: best orientation explains only %d spots (< %d) - no further grain\n', ...
            k, nk, S.min_spots_per_grain);
        break
    end
    Tk.spot = reshape(rest(Tk.spot), [], 1);           % back to rows of expt
    Tk.grain = k * ones(height(Tk), 1);
    m = misorientation(result.q, result.phi0, rk.q, rk.phi0, phase.parent_rotations);
    grains(k) = struct('q', rk.q, 'phi0', rk.phi0, 'T', Tk, 'calc', ck, ...
                       'n_spots', nk, 'misorientation_to_1', m); %#ok<AGROW>
    used = [used; unique(Tk.spot)]; %#ok<AGROW>
    fprintf('Grain %d: q = [%.5f %.5f %.5f %.5f], phi0 = %.2f deg, %d spots, misorientation to grain 1 = %.1f deg\n', ...
        k, rk.q, rk.phi0, nk, m);
end

%% 7. Plot and save =======================================================
for k = 1:numel(grains)
    grains(k).calc = grains(k).calc(grains(k).calc(:,1) < plot_max_tth, :);
end
if numel(grains) == 1
    ttl = sprintf('%s: best orientation', name);
else
    ttl = sprintf('%s: %d grains', name, numel(grains));
end
fig = plot_overlay(img_ll, geom, peaks_px, grains, beam_center, ttl);
exportgraphics(fig, fullfile(out_dir, 'orientation_overlay.png'), 'Resolution', 200);

T_all = vertcat(grains.T);
T_all = movevars(T_all, 'grain', 'Before', 1);
writetable(array2table(expt, 'VariableNames', {'tth_deg', 'beta_deg'}), fullfile(out_dir, 'spots_measured.csv'));
writetable(T_all, fullfile(out_dir, 'spots_matched.csv'));
grains_saved = rmfield(grains, 'calc');
save(fullfile(out_dir, 'orientation.mat'), 'result', 'grains_saved', 'expt', 'peaks_px', 'phase', ...
    'S', 'P', 'beam_center', 'bc_log', 'cell_log', 'image_file');

fid = fopen(fullfile(out_dir, 'summary.txt'), 'w');
fprintf(fid, 'Frame:            %s\n', image_file);
fprintf(fid, 'Phase:            %s\n', phase.name);
fprintf(fid, 'Rotation range:   %g to %g deg around phi0\n', S.rotation_range);
fprintf(fid, 'Spots found:      %d\n', size(expt, 1));
fprintf(fid, 'Cell of the CIF:  %.4f %.4f %.4f %g %g %g\n', cell_cif);
fprintf(fid, 'Beam centre:      [%.2f %.2f] px (row, column)', beam_center);
if S.refine_beam_center && ~isempty(bc_log)
    fprintf(fid, ' +- [%.2f %.2f], refined from [%.2f %.2f]\n', bc_log(end, 4:5), beam_center_start);
else
    fprintf(fid, '\n');
end
for k = 1:size(cell_log, 1)
    fprintf(fid, 'Refined cell %d:   a b c = %.4f %.4f %.4f A (esd %.4f %.4f %.4f), rms d2theta %.3f -> %.3f deg\n', ...
        cell_log(k, :));
end
for k = 1:numel(grains)
    Tk = grains(k).T;
    fprintf(fid, '\nGrain %d\n', k);
    fprintf(fid, '  Orientation q:  [%.6f %.6f %.6f %.6f]\n', grains(k).q);
    fprintf(fid, '  phi0:           %.2f deg\n', grains(k).phi0);
    fprintf(fid, '  Matched spots:  %d predicted reflections, %d measured spots\n', height(Tk), grains(k).n_spots);
    fprintf(fid, '  Residuals:      d2theta %.3f +- %.3f deg, dbeta %.3f +- %.3f deg\n', ...
        mean(Tk.d_tth), std(Tk.d_tth), mean(Tk.d_beta), std(Tk.d_beta));
    if k > 1
        fprintf(fid, '  Misorientation to grain 1: %.1f deg\n', grains(k).misorientation_to_1);
    end
end
fprintf(fid, '\nUnexplained spots: %d of %d\n', size(expt, 1) - numel(unique(T_all.spot)), size(expt, 1));
fclose(fid);
fprintf('Saved to %s\n', out_dir);

out = struct('result', result, 'expt', expt, 'T', T, 'out_dir', out_dir, ...
    'cell', phase.cell, 'cell_log', cell_log, 'beam_center', beam_center, 'bc_log', bc_log, ...
    'grains', rmfield(grains, 'calc'));
end

%% ------------------------------------------------------------------------
function m = misorientation(q1, phi1, q2, phi2, R)
% smallest rotation angle between two orientations, allowing for the
% rotations R of the crystal
Rz = @(a) [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
L1 = Rz(phi1) * quat2dcm(q1)';
L2 = Rz(phi2) * quat2dcm(q2)';
m = inf;
for k = 1:size(R, 3)
    D = L1' * L2 * R(:,:,k);
    m = min(m, acosd(max(-1, min(1, (trace(D) - 1) / 2))));
end
end


%% ------------------------------------------------------------------------
function [best, n_best] = scan_beam_center(bc0, d, peaks_px, result, g, hkl, variant, S)
% beam centre on a grid of +-S.beam_center_scan px (step 1 px, refined from
% a 4 px grid) around bc0 with the most measured spots matched by the
% orientation result; ties are resolved by the smallest rms residual
best = bc0;  n_best = -1;  rms_best = inf;
for step = [4 1]
    if step == 4, w = S.beam_center_scan; else, w = 4; end
    centre = best;
    for dr = -w:step:w
        for dc = -w:step:w
            bc = centre + [dr dc];
            d.zero = bc;
            geom = RAPIDutils.getDiffractionGeometry(d);
            [tth, beta] = RAPIDutils.ipPixels2diff(geom, peaks_px(:,2), peaks_px(:,1));
            T = match_spots(result.q, result.phi0, g, hkl, variant, [tth, beta] * 180/pi, S);
            n = numel(unique(T.spot));  rms = sqrt(mean(T.d_tth.^2));
            if n > n_best || (n == n_best && rms < rms_best)
                best = bc;  n_best = n;  rms_best = rms;
            end
        end
    end
end
end

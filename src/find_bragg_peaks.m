function [expt, peaks_px, img_ll, info] = find_bragg_peaks(d, geom, P)
%FIND_BRAGG_PEAKS  Positions of the Bragg spots on an image-plate frame.
%
%   [expt, peaks_px, img_ll, info] = find_bragg_peaks(d, geom, P)
%
%   d, geom : frame from RAPIDread (d.zero = beam centre) and its geometry
%   P       : peak-search settings (see run_orientation_search.m)
%             P.min_height, P.min_prominence_map, P.min_prominence_pixel,
%             P.min_distance_px, P.refine_radius_px, P.tth_beta_map
%
%   The frame is compressed to log(log(I + 1)). A pixel is a peak when it
%   is a local maximum (MATLAB findpeaks) both along its row and along its
%   column. The search is done twice and the results merged:
%     1) in the frame resampled to a regular (2theta, beta) grid; every
%        hit is mapped back to the detector and moved to the brightest
%        pixel within P.refine_radius_px;
%     2) directly in detector pixels.
%   Of two peaks closer than P.min_distance_px only the brighter is kept.
%
%   expt     : N x 2 peak positions [2theta, beta] (deg), sorted by 2theta
%   peaks_px : N x 3 [row, column, log(log(I+1))] of the same spots (row k
%              of peaks_px is spot k of expt)
%   img_ll   : the log(log(I+1)) image (for plotting)
%   info     : intermediate results (for checking the search, see
%              step_by_step.m): from_map, from_pixels ([row, column, value]
%              found by search 1 and 2), map, tth_axis, beta_axis (the
%              (2theta, beta) map and its axes, rad)

img_ll = log(log(double(d.image) + 1));

% 1) search in the (2theta, beta) map
d_ll = d;
d_ll.image = img_ll;
m = P.tth_beta_map;
[map, tth_axis, beta_axis] = RAPIDutils.getTThBtIntensityMap(d_ll, m(1), m(2), m(3), m(4));
hits = peaks_2d(map, P.min_height, P.min_prominence_map, P.min_distance_px);
[ipx, ipy] = RAPIDutils.diff2ipPixels(geom, tth_axis(hits(:,2)), beta_axis(hits(:,1)));
from_map = brightest_pixel(img_ll, ipx', ipy', P.refine_radius_px);

% 2) search in detector pixels
from_pixels = peaks_2d(img_ll, P.min_height, P.min_prominence_pixel, P.min_distance_px);

% merge
peaks_px = remove_close([from_map; from_pixels], P.min_distance_px);
peaks_px = unique(peaks_px, 'rows');
[~, order] = sort(peaks_px(:,3), 'descend');
peaks_px = peaks_px(order, :);

[tth, beta] = RAPIDutils.ipPixels2diff(geom, peaks_px(:,2), peaks_px(:,1));
expt = [tth, beta] * 180 / pi;
[~, order] = sort(expt(:,1));
expt = expt(order, :);
peaks_px = peaks_px(order, :);
info = struct('from_map', from_map, 'from_pixels', from_pixels, 'map', map, ...
              'tth_axis', tth_axis, 'beta_axis', beta_axis);
end

%% ------------------------------------------------------------------------
function peaks = peaks_2d(img, min_height, min_prominence, min_dist)
% local maxima along rows AND along columns: [row, column, value]
[pk_r, loc_r] = findpeaks(double(img(:)), 'MinPeakHeight', min_height, ...
    'MinPeakProminence', min_prominence);
img_t = img';
[pk_c, loc_c] = findpeaks(double(img_t(:)), 'MinPeakHeight', min_height, ...
    'MinPeakProminence', min_prominence);
[r1, c1] = ind2sub(size(img), loc_r);
[c2, r2] = ind2sub(size(img_t), loc_c);
peaks = intersect([r1, c1, pk_r], [r2, c2, pk_c], 'rows');
peaks = remove_close(peaks, min_dist);
end

function peaks = remove_close(peaks, min_dist)
% of every pair of peaks closer than min_dist remove the weaker one
dist = triu(pdist2(peaks(:,1:2), peaks(:,1:2)));
[r, c] = find((dist < min_dist) & (dist > eps));
remove = c .* (peaks(r,3) > peaks(c,3)) + r .* (peaks(r,3) <= peaks(c,3));
peaks(remove, :) = [];
end

function out = brightest_pixel(img, x, y, radius)
% brightest pixel within +-radius around each (x, y): [row, column, value]
out = zeros(numel(x), 3);
for k = 1:numel(x)
    x0 = round(max(x(k) - radius, 1));  x1 = round(min(x(k) + radius, size(img, 2)));
    y0 = round(max(y(k) - radius, 1));  y1 = round(min(y(k) + radius, size(img, 1)));
    win = img(y0:y1, x0:x1);
    [val, idx] = max(win(:));
    [r, c] = ind2sub(size(win), idx);
    out(k, :) = [y0 + r - 1, x0 + c - 1, val];
end
end

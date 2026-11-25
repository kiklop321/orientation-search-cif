function fig = plot_overlay(img_ll, geom, peaks_px, grains, beam_center, title_str)
%PLOT_OVERLAY  Frame with measured spots and predicted reflections (cf. Fig. 13).
%
%   fig = plot_overlay(img_ll, geom, peaks_px, grains, beam_center, title_str)
%
%   grains : struct array with fields calc (predicted spots), T (matched
%            spots) and n_spots for every grain (see ORIENTATION_PIPELINE)
%
%   Red +   measured spots
%   One grain (style of the article):
%     black o   predicted reflection with a measured spot within the tolerance
%     green o   predicted reflection without a measured spot
%   Several grains: one colour per grain (black, magenta, cyan, yellow, ...)
%     thick o   predicted reflection of that grain with a measured spot
%     thin  <>  predicted reflection of that grain without a measured spot
%   "Predicted" = recorded during the rotation, 2theta below the plot limit
%   and on the image plate.

colours = [0 0 0; 1 0 1; 0 1 1; 1 1 0; 1 0.5 0; 1 1 1];
on = @(x, z) x >= 1 & x <= size(img_ll, 2) & z >= 1 & z <= size(img_ll, 1);
single = numel(grains) == 1;

fig = figure('Color', 'w', 'Position', [100 100 1400 760]);
imagesc(img_ll); colormap(parula); axis image; hold on
plot(peaks_px(:,2), peaks_px(:,1), 'r+', 'DisplayName', 'measured spots');

for k = 1:numel(grains)
    c = grains(k).calc;
    T = grains(k).T;
    matched = ismember(c(:, [1 2]), [T.tth_calc, T.beta_calc], 'rows');
    [ux, uz] = RAPIDutils.diff2ipPixels(geom, c(~matched, 1) * pi/180, c(~matched, 2) * pi/180);
    [mx, mz] = RAPIDutils.diff2ipPixels(geom, T.tth_calc * pi/180, T.beta_calc * pi/180);
    ku = on(ux, uz);  km = on(mx, mz);
    if single
        plot(ux(ku), uz(ku), 'go', 'LineWidth', 1, 'DisplayName', 'predicted, no spot within tolerance');
        plot(mx(km), mz(km), 'ko', 'LineWidth', 1.5, 'DisplayName', 'predicted, matched');
    else
        col = colours(min(k, size(colours, 1)), :);
        plot(ux(ku), uz(ku), 'd', 'Color', col, 'LineWidth', 0.75, 'MarkerSize', 5, ...
            'DisplayName', sprintf('grain %d: predicted, no spot', k));
        plot(mx(km), mz(km), 'o', 'Color', col, 'LineWidth', 1.75, 'MarkerSize', 8, ...
            'DisplayName', sprintf('grain %d: matched (%d spots)', k, grains(k).n_spots));
    end
end

plot(beam_center(2), beam_center(1), 'c+', 'MarkerSize', 12, 'DisplayName', 'beam centre');
xlim([1 size(img_ll, 2)]); ylim([1 size(img_ll, 1)]);
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', min(3 + 2 * numel(grains), 4));
xlabel('x (px)'); ylabel('z (px)');
title(title_str, 'Interpreter', 'none');
end

function describe_figure(fig, paragraphs)
%DESCRIBE_FIGURE  Explanatory text below the plots of a figure.
%
%   describe_figure(fig, paragraphs)
%
%   paragraphs : cell array of strings, one per paragraph
%
%   The figure is made taller, all axes are moved up into the original area
%   and the text is put into the space below them. Calling it again on the
%   same figure (e.g. after re-running a section) does not grow the figure
%   further.

if ~isappdata(fig, 'base_height')
    setappdata(fig, 'base_height', fig.Position(4));
end
H = getappdata(fig, 'base_height');
delete(findall(fig, 'Tag', 'description'));

% height of the text: about 7 px per character, 15 px per line
width_px = fig.Position(3);
chars_per_line = max(floor((width_px - 30) / 6.5), 40);
n_lines = sum(cellfun(@(p) ceil(numel(p) / chars_per_line), paragraphs)) + numel(paragraphs) - 1;
h = n_lines * 15 + 20;

fig.Position(2) = max(fig.Position(2) - (H + h - fig.Position(4)), 0);
fig.Position(4) = H + h;
% SUBPLOT listens to its axes and would move the last one back to its grid
% position at the next redraw; remove its listeners, the axes are placed here
if isappdata(fig, 'SubplotListenersManager')
    delete(getappdata(fig, 'SubplotListenersManager'));
    rmappdata(fig, 'SubplotListenersManager');
end
squeeze = @(p) [p(1), h / (H + h) + p(2) * H / (H + h), p(3), p(4) * H / (H + h)];
for ax = findobj(fig, 'Type', 'axes')'
    if ax.Parent ~= fig, continue; end           % axes in a tiledlayout: moved with it below
    ax.Units = 'normalized';
    ax.Position = squeeze(ax.Position);
end
for tl = findobj(fig, 'Type', 'tiledlayout')'
    if tl.Parent ~= fig, continue; end
    tl.Units = 'normalized';
    tl.OuterPosition = squeeze(tl.OuterPosition);
end
annotation(fig, 'textbox', [0.01 0.005 0.98 (h - 5) / (H + h)], 'String', strjoin(paragraphs, [newline newline]), ...
    'Interpreter', 'none', 'FontSize', 9, 'EdgeColor', 'none', 'VerticalAlignment', 'top', ...
    'Tag', 'description');
end

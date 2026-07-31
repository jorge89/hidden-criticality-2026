function violinBoxByRecording(ax, V, labels, recColors, recMarkers, varargin)
% violinBoxByRecording(ax, V, labels, recColors, recMarkers)
% V: [nRec x nGroups] values (rows=recordings, cols=signals)
% labels: 1 x nGroups cellstr
% recColors: [nRec x 3] RGB per recording
% recMarkers: 1 x nRec cellstr markers per recording (e.g., {'o','s','^',...})
%
% Optional name-value:
%   'YScale' : 'linear' (default) or 'log'
%   'Title'  : ''
%   'YLabel' : ''
%   'ShowLegend' : false

p = inputParser;
addParameter(p, 'YScale', 'linear');
addParameter(p, 'Title', '');
addParameter(p, 'YLabel', '');
addParameter(p, 'ShowLegend', false);
parse(p, varargin{:});

hold(ax, 'on');
box(ax, 'off');
grid(ax, 'on');

[nRec, nG] = size(V);
xpos = 1:nG;

% Violin appearance
violinColor = [0.87 0.87 0.87];
violinEdge  = 'none';
violinAlpha = 1.0;
maxWidth = 0.35;      % half-width of violin

% Point appearance
jit = 0.10;
ms = 38;              % point size

% Summary appearance: median + IQR
medSize = 90;         % median dot size
iqrLW  = 5;           % IQR bar line width

% Helper: is marker "fillable"?
fillable = @(mk) ismember(mk, {'o','s','d','^','v','>','<','p','h'});

for g = 1:nG
    yAll = V(:,g);
    yAll = yAll(:);

    % --- violin KDE built from finite values ---
    y = yAll(isfinite(yAll));
    if numel(y) < 3
        continue;
    end

    if strcmpi(p.Results.YScale, 'log')
        y = y(y > 0);
        if numel(y) < 3, continue; end
        z = log10(y);
        [f, zz] = ksdensity(z, 'Function', 'pdf', 'NumPoints', 200);
        yy = 10.^zz;
    else
        [f, yy] = ksdensity(y, 'Function', 'pdf', 'NumPoints', 200);
    end

    f = f / max(f);
    w = maxWidth * f;

    Xv = [xpos(g)-w, fliplr(xpos(g)+w)];
    Yv = [yy,        fliplr(yy)];
    patch(ax, Xv, Yv, violinColor, 'EdgeColor', violinEdge, 'FaceAlpha', violinAlpha);

    % --- per-recording points: color + marker shape ---
    for r = 1:nRec
        yr = yAll(r);
        if ~isfinite(yr), continue; end
        if strcmpi(p.Results.YScale, 'log') && yr <= 0, continue; end

        xr = xpos(g) + jit * randn(1);
        mk = recMarkers{r};

        if fillable(mk)
            scatter(ax, xr, yr, ms, ...
                'Marker', mk, ...
                'MarkerFaceColor', recColors(r,:), ...
                'MarkerEdgeColor', recColors(r,:), ...
                'MarkerFaceAlpha', 0.35, ...
                'MarkerEdgeAlpha', 0.8);
        else
            % non-fillable markers like 'x','+','*'
            scatter(ax, xr, yr, ms, ...
                'Marker', mk, ...
                'MarkerEdgeColor', recColors(r,:), ...
                'LineWidth', 1.2);
        end
    end

    % --- median + IQR (Q1..Q3) ---
    yUse = yAll(isfinite(yAll));
    if strcmpi(p.Results.YScale, 'log')
        yUse = yUse(yUse > 0);
    end
    if numel(yUse) >= 3
        q1 = prctile(yUse, 25);
        med = prctile(yUse, 50);
        q3 = prctile(yUse, 75);

        plot(ax, [xpos(g) xpos(g)], [q1 q3], '-', 'LineWidth', iqrLW, 'Color', [0 0 0]);
        scatter(ax, xpos(g), med, medSize, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
    end
end

% Axes styling
ax.XLim = [0.5 nG+0.5];
ax.YLim = [0 4];
ax.XTick = xpos;
ax.XTickLabel = labels;
ax.FontSize = 14;
ax.LineWidth = 1.2;
set(ax, 'YScale', p.Results.YScale);

if ~isempty(p.Results.Title)
    title(ax, p.Results.Title, 'FontWeight','bold');
end
if ~isempty(p.Results.YLabel)
    ylabel(ax, p.Results.YLabel);
end

% Optional legend (one entry per recording: can be huge)
if p.Results.ShowLegend
    h = gobjects(nRec,1);
    for r=1:nRec
        mk = recMarkers{r};
        h(r) = plot(ax, nan, nan, ...
            'LineStyle','none', 'Marker', mk, ...
            'MarkerFaceColor', recColors(r,:), ...
            'MarkerEdgeColor', recColors(r,:), ...
            'MarkerSize', 7);
    end
    legend(ax, h, arrayfun(@(r)sprintf('Rec %d',r), 1:nRec, 'UniformOutput',false), ...
        'Location','eastoutside');
end
end

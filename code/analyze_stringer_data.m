%% analyze_stringer_data.m
%
% Analyzes the mouse visual cortex recordings (Stringer et al., 2019;
% doi:10.25378/janelia.6163622.v4) to produce the panels underlying
% Fig. 4 (main text).
%
% Requires the following .mat files to be present in the working
% directory (downloaded separately; not included in this repository —
% see data/README.md):
%   spont_*.mat   (one file per recording, each containing variable "Fsp")
%
% Requires on the MATLAB path:
%   checkCriticalitySignaturesMulti.m
%   violinBoxByRecording.m
%
% NOTE: Sections 1 and 3 below appear to overlap with parts of Section 2
% and Section 4 respectively (single-recording PCA test vs. the full
% batch loop; generic boxplot vs. the per-recording colored-symbol plot
% matching the published Fig. 4g-j style). Both are kept here since it
% was not confirmed whether they were used for anything beyond an
% intermediate diagnostic step. Confirm before removing.

clear; close all; clc;

%% ---- Config ----
dataDir = pwd;
varName = 'Fsp';           % activity matrix variable name inside each .mat
k = 50;                    % number of PCs to compute
nPCplot = 5;                % number of PCs used for criticality analysis
thresholdQuantile = 0.6;    % avalanche threshold quantile
SF = 2.5;                   % sampling frequency (Hz)
fitPCmin = 2;               % PC range used for explained-variance power-law fit
fitPCmax = 50;
saveResults = true;

%% ---- Section 1 (diagnostic): PCA and explained variance for one recording ----
% Loads a single Fsp matrix already in the workspace and inspects PCA
% explained variance and the first few PC timecourses. This appears to
% duplicate the explained-variance step performed automatically for
% every recording in Section 2 below; kept here as it may document the
% rationale for choosing k and nPCplot.

if exist('Fsp', 'var')
    [nNeurons, nSamples] = size(Fsp); %#ok<*ASGLU>
    fprintf('FR size = %d neurons x %d samples\n', nNeurons, nSamples);
    Fsp = double(Fsp);

    mu    = nanmean(Fsp, 2);
    sigma = nanstd(Fsp, 0, 2);
    sigma(sigma == 0) = 1;
    Z = (Fsp - mu) ./ sigma;
    popSum = sum(Z, 1);

    X = single(Z).';
    opts = struct('tol', 1e-4, 'maxit', 500, 'disp', 1);
    [U, S, V] = svds(double(X), k, 'largest', opts);

    score  = U * S;
    latent = diag(S).^2 / (size(X, 1) - 1);
    explained = 100 * latent / sum(var(double(X), 0, 1));

    figure;
    x = 1:length(explained);
    y = explained(:)';
    valid = y > 0 & x > 0;
    p = polyfit(log10(x(valid)), log10(y(valid)), 1);
    alpha = p(1); %#ok<NASGU>
    xf = logspace(log10(min(x)), log10(max(x)), 200);
    yf = 10.^polyval(p, log10(xf));
    loglog(x, y, 'o', 'MarkerSize', 5, 'LineWidth', 1.2); hold on;
    loglog(xf, yf, '--r', 'LineWidth', 2);
    xlabel('PC #'); ylabel('Explained variance (%)');
    title('PCA explained variance'); grid on; box off;
    annotation('textbox', [0.60 0.70 0.25 0.12], ...
        'String', sprintf('Explained \\sim PC^{%.2f}', alpha), ...
        'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 11);

    figure;
    plot(score(:, 1:5));
    xlabel('Sample'); ylabel('Score'); title('First PCs (timecourses)');
    legend(arrayfun(@(i) sprintf('PC%d', i), 1:min(5, k), 'UniformOutput', false));
else
    disp('Section 1 skipped: no Fsp matrix in workspace.');
end

%% ---- Section 2: batch analysis across all recordings ----
% For each spont_*.mat recording: z-score, run PCA, build population-sum
% and PC1-PC5 signals, and compute avalanche/DFA/PSD/ACF criticality
% measures for each. This generates the per-recording example panels
% (Fig. 4b-f style) and stores all metrics for the across-recording
% summary in Section 4.

files = dir(fullfile(dataDir, 'spont_*.mat'));
if isempty(files)
    error('No spont_*.mat files found in %s', dataDir);
end

results = struct();

for r = 1:numel(files)
    disp(r)
    fname = fullfile(files(r).folder, files(r).name);
    S = load(fname);

    if ~isfield(S, varName)
        warning('Skipping %s (missing variable "%s").', files(r).name, varName);
        continue;
    end

    Fsp = double(S.(varName));
    [nNeurons, nSamples] = size(Fsp);

    mu    = nanmean(Fsp, 2);
    sigma = nanstd(Fsp, 0, 2);
    sigma(sigma == 0) = 1;
    Z = (Fsp - mu) ./ sigma;

    popSum = sum(Z, 1);
    X = single(Z).';

    opts = struct('tol', 1e-4, 'maxit', 500, 'disp', 0);
    [U, Sig, V] = svds(double(X), k, 'largest', opts); %#ok<ASGLU>

    score  = U * Sig;
    latent = diag(Sig).^2 / (size(X, 1) - 1);
    explained = 100 * latent / sum(var(double(X), 0, 1));

    % Build z-scored population-sum and PC1-PC5 signals
    pcMat = score(:, 1:min(nPCplot, k));
    popZ  = (popSum(:) - mean(popSum)) / std(popSum);
    pcZ   = pcMat;
    for j = 1:size(pcZ, 2)
        pcZ(:, j) = (pcZ(:, j) - mean(pcZ(:, j))) / std(pcZ(:, j));
    end

    signals = cell(1, 1 + size(pcZ, 2));
    labels  = cell(1, 1 + size(pcZ, 2));
    signals{1} = popZ(:)'; labels{1} = 'popSum';
    for j = 1:size(pcZ, 2)
        signals{1+j} = pcZ(:, j)';
        labels{1+j}  = sprintf('PC%d', j);
    end

    % ---- Build per-recording panel figure (Fig. 4b-f style) ----
    fig = figure('Color', 'w', 'Name', files(r).name, 'NumberTitle', 'off');
    tl = tiledlayout(fig, 2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

    axExpl  = nexttile(tl, 1, [1 2]);
    axTrace = nexttile(tl, 3, [1 4]);
    axSize  = nexttile(tl, 7);
    axDFA   = nexttile(tl, 8);
    axPSD   = nexttile(tl, 9);
    axACF   = nexttile(tl, 10);

    axExpl.PlotBoxAspectRatio  = [1 1 1];
    axSize.PlotBoxAspectRatio  = [1 1 1];
    axDFA.PlotBoxAspectRatio   = [1 1 1];
    axPSD.PlotBoxAspectRatio   = [1 1 1];
    axACF.PlotBoxAspectRatio   = [1 1 1];
    axTrace.PlotBoxAspectRatio = [3 1 1];

    % Explained variance + power-law fit (diagnostic; not itself a
    % published panel — see note at top of file)
    x = (1:numel(explained))';
    y = explained(:);
    pcMaxUse = min([fitPCmax, numel(y), k]);
    fitMask = (x >= fitPCmin) & (x <= pcMaxUse) & (y > 0);
    p = polyfit(log10(x(fitMask)), log10(y(fitMask)), 1);
    alpha = p(1);
    xf = logspace(log10(min(x(fitMask))), log10(max(x(fitMask))), 200);
    yf = 10.^polyval(p, log10(xf));

    loglog(axExpl, x, y, 'o', 'MarkerSize', 4, 'LineWidth', 1); hold(axExpl, 'on');
    loglog(axExpl, xf, yf, '--', 'LineWidth', 2);
    grid(axExpl, 'on'); box(axExpl, 'off');
    xlabel(axExpl, 'PC #'); ylabel(axExpl, 'Explained var (%)');
    title(axExpl, 'Explained variance');
    text(axExpl, 0.05, 0.08, sprintf('\\alpha = %.2f  (PC %d-%d)', alpha, fitPCmin, pcMaxUse), ...
        'Units', 'normalized', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3, 'FontSize', 10);

    title(tl, sprintf('%s | %d neurons x %d samples', files(r).name, nNeurons, nSamples), ...
        'Interpreter', 'none', 'FontWeight', 'bold');

    % Stacked population-sum + PC1-5 traces (Fig. 4b)
    t = (1:nSamples)';
    gap = 3;
    offsets = gap * (0:(size(pcZ, 2) - 1));
    plot(axTrace, t, popZ, 'k', 'LineWidth', 1.3); hold(axTrace, 'on');
    for j = 1:size(pcZ, 2)
        plot(axTrace, t, pcZ(:, j) + offsets(j) + gap, 'LineWidth', 1.0);
    end
    xlabel(axTrace, 'Sample'); ylabel(axTrace, 'z-score (stacked)');
    title(axTrace, 'popSum + PC1..PC5');
    ytick = [0, offsets + gap];
    ylab = [{'popSum'}, arrayfun(@(j) sprintf('PC%d', j), 1:size(pcZ, 2), 'UniformOutput', false)];
    set(axTrace, 'YTick', ytick, 'YTickLabel', ylab);
    grid(axTrace, 'on'); box(axTrace, 'off');

    % Avalanche / DFA / PSD / ACF panels (Fig. 4c-f)
    outCrit = checkCriticalitySignaturesMulti(signals, labels, thresholdQuantile, SF, ...
        axSize, axDFA, axPSD, axACF);

    if saveResults
        recID = matlab.lang.makeValidName(files(r).name);
        results.(recID).file = files(r).name;
        results.(recID).alphaExplained = alpha;
        results.(recID).explained = explained;
        results.(recID).crit = outCrit;
    end
end

if saveResults
    save(fullfile(dataDir, 'criticality_batch_results.mat'), 'results');
end

%% ---- Section 3 (diagnostic): generic boxplot summary across recordings ----
% Standard MATLAB boxplot of each criticality metric across all
% recordings and signals (popSum, PC1-PC5). Superseded in appearance by
% Section 4 below, which matches the colored-symbol + median/IQR style
% used in the published Fig. 4g-j; kept here in case it was used as an
% intermediate check. Confirm before removing.

recNames = fieldnames(results);
nRec = numel(recNames);
if nRec == 0
    error('results is empty.');
end

firstCrit = results.(recNames{1}).crit;
labels = firstCrit.labels;
nSig = numel(labels);

metrics = {'tau', 'tauT', 'rho', 'rhoT', 'rsnz', 'corrTime', 'alfa', 'beta'};
metricTitles = {'\tau (size exponent)', '\tau_T (duration exponent)', ...
    '\rho (size PL range)', '\rho_T (duration PL range)', ...
    '1/(\sigma\nu z) (rsnz)', 'corrTime (lag@ACF<0.1)', ...
    '\alpha (PSD slope)', '\beta (DFA slope)'};

nM = numel(metrics);
nCols = 4;
nRows = ceil(nM / nCols);

fig = figure('Color', 'w');
tl = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Criticality metrics across %d recordings', nRec), 'FontWeight', 'bold');

for m = 1:nM
    metric = metrics{m};
    Vmat = nan(nRec, nSig);
    for r = 1:nRec
        crit = results.(recNames{r}).crit;
        if ~isequal(crit.labels, labels)
            [tf, idx] = ismember(labels, crit.labels);
            if ~all(tf)
                warning('Recording %s missing some labels. Filling with NaN.', recNames{r});
            end
        else
            idx = 1:nSig; tf = true(1, nSig); %#ok<NASGU>
        end
        if isfield(crit, metric)
            val = crit.(metric);
            if ~isequal(crit.labels, labels)
                tmp = nan(1, nSig);
                tmp(tf) = val(idx(tf));
                val = tmp;
            end
            Vmat(r, :) = val(:).';
        else
            warning('Metric "%s" not found in recording %s', metric, recNames{r});
        end
    end

    ax = nexttile(tl, m);
    y = Vmat(:);
    g = repelem(1:nSig, nRec)';
    keep = isfinite(y);
    y = y(keep); g = g(keep);
    boxplot(ax, y, g, 'Labels', labels, 'Symbol', 'k.');
    grid(ax, 'on'); box(ax, 'off');
    title(ax, metricTitles{m}, 'Interpreter', 'tex');
    if ismember(metric, {'corrTime'})
        set(ax, 'YScale', 'log');
    end
end

%% ---- Section 4: per-recording colored-symbol summary (Fig. 4g-j) ----
% Matches the published figure style: one colored symbol per recording,
% median +/- IQR shown as a black oval, for each criticality metric.

recColors = lines(nRec);
markerPool = {'o', 's', '^', 'v', 'd', '>', '<', 'p', 'h', 'x', '+'};
recMarkers = cell(1, nRec);
for r = 1:nRec
    recMarkers{r} = markerPool{mod(r-1, numel(markerPool)) + 1};
end

metrics = {'tau', 'tauT', 'rho', 'rhoT', 'rsnz', 'corrTime', 'alfa', 'beta'};
metricTitles = {'\tau', '\tau_T', '\rho', '\rho_T', 'rsnz', 'corrTime', '\alpha (PSD)', '\beta (DFA)'};

fig = figure('Color', 'w');
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(sprintf('Criticality metrics across %d recordings (median + IQR)', nRec), 'FontWeight', 'bold');

for m = 1:numel(metrics)
    metric = metrics{m};
    Vmat = nan(nRec, nSig);
    for r = 1:nRec
        crit = results.(recNames{r}).crit;
        if isfield(crit, metric)
            Vmat(r, :) = crit.(metric)(:).';
        end
    end

    ax = nexttile(tl, m);
    if strcmp(metric, 'corrTime')
        violinBoxByRecording(ax, Vmat, labels, recColors, recMarkers, ...
            'YScale', 'log', 'Title', metricTitles{m}, 'YLabel', 'Value');
    else
        violinBoxByRecording(ax, Vmat, labels, recColors, recMarkers, ...
            'YScale', 'linear', 'Title', metricTitles{m}, 'YLabel', 'Value');
    end
end

% Export the Section 4 summary figure (Fig. 4g-j) as a vector file.
% NOTE: previously used a hardcoded figure(12) handle to select this
% figure for export, which is fragile if figure numbering changes;
% replaced with the actual handle created above.
exportgraphics(fig, 'criticality_boxpanel.svg', 'ContentType', 'vector');

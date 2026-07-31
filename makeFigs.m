%% makeFigs.m
%
% Generates the plots underlying Fig. 2 (main text), Fig. 3 (main text),
% and Supplementary Fig. S2/S3, from the quantities produced by
% compute_criticality_quantities.m.
%
% Run compute_criticality_quantities.m first (once per case) so that the
% following variables exist in the workspace:
%   eigvals, IdxX, Idx, dataSum, dataModes, dataModesEnv,
%   rho, rhoT, alfa, beta, corrTime, d2Env,
%   rhoDS, rhoTDS, tauDS, tauTDS, rsnzDS, corrTimeDS, alfaDS, betaDS,
%   varAlfa, SumV, simDrop   (see Section 0 below)
%
% Requires plotCrit.m on the MATLAB path.

%% ---- Section 0: quantities used only for plotting (not already saved) ----
% These support the Supplementary Fig. S2 panels and are computed here
% rather than in compute_criticality_quantities.m, since they are purely
% plotting-support quantities (component sums, variances, and the
% leave-one-mode-out impact on the population sum).

recOns = zeros(size(dataModes));
varAlfa = zeros(size(IdxX));
SumV    = zeros(size(IdxX));
ii = 1;
for i = IdxX
    if isreal(V(:, i))
        recOns(ii, :) = dataModes(ii, :) * sum(V(:, i));
        varAlfa(ii)   = var(dataModes(ii, :));
    else
        recOns(ii, :) = 2 * real(dataModes(ii, :) * sum(V(:, i)));
        varAlfa(ii)   = var(dataModes(ii, :));
    end
    SumV(ii) = sum(V(:, i));
    ii = ii + 1;
end

num_modes = size(recOns, 1);
simDrop = zeros(1, num_modes);
for i = 1:num_modes
    recOns_removed = recOns;
    recOns_removed(i, :) = 0;
    partialSum = sum(recOns_removed, 1);
    r = corr(dataSum(:), partialSum(:));
    simDrop(i) = 1 - abs(r);
end

%% ---- Figure: Supplementary Fig. S2 (mode-level predictors of hidden criticality) ----
% Panel a: alpha(t) variance and loading component sum vs Re(lambda)
% Panel b: impact on population average vs |loading component sum|
figure;
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(1); % panel a
yyaxis left
semilogy(real(eigvals(IdxX)), varAlfa, '*');
ylabel('\alpha(t) Variance');
yyaxis right
plot(real(eigvals(IdxX)), real(SumV));
ylabel('v_i Component Sum');
xlabel('Re(\lambda_i)');

nexttile(2); % panel b
colorValues = real(eigvals(IdxX));
scatter(abs(real(SumV)), simDrop, 100, colorValues, 'filled');
colormap('jet');
c = colorbar;
c.Label.String = '\lambda';
caxis([-1 1]);
set(gca, 'XScale', 'log');
xlabel('|v_i Component Sum|');
ylabel('Impact on Pop Avg');

sgtitle('Supplementary Fig. S2 panels (run once per case, Case 2 -> a,b and Case 3 -> c,d)');

%% ---- Figure: Fig. 2 (main text) + Supplementary Fig. S3 (same case) ----
% Fig. 2 panels: avalanche size distribution, PSD, fluctuation function (DFA), autocorrelation
% Supp. Fig. S3 panels: avalanche duration distribution, avalanche size vs duration
figure;
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

modeSpecs = { ...
    dataModesEnv(1, :),      'green',   1;    % bulk / farthest-from-critical mode
    dataSum,                 'k',       10;   % population average
    dataModesEnv(end-1, :),  'magenta', 100;  % second-nearest-to-critical mode
    dataModesEnv(end, :),    'b',       1000  % near-critical mode
    };

nexttile(1); hold on; % Fig. 2: avalanche size distribution
for m = 1:size(modeSpecs, 1)
    plotCrit(modeSpecs{m,1}, 0.5, 22, 'sizeDist', modeSpecs{m,2}, modeSpecs{m,3});
end
hold off;

nexttile(2); hold on; % Supp. Fig. S3: avalanche duration distribution
for m = 1:size(modeSpecs, 1)
    plotCrit(modeSpecs{m,1}, 0.5, 22, 'durationDist', modeSpecs{m,2}, modeSpecs{m,3});
end
hold off;

nexttile(3); hold on; % Supp. Fig. S3: avalanche size vs duration (crackling noise scaling)
for m = 1:size(modeSpecs, 1)
    plotCrit(modeSpecs{m,1}, 0.5, 22, 'sizeVsDuration', modeSpecs{m,2}, modeSpecs{m,3});
end
hold off;

nexttile(4); hold on; % Fig. 2: fluctuation function (DFA)
for m = 1:size(modeSpecs, 1)
    plotCrit(modeSpecs{m,1}, 0.5, 22, 'dfa', modeSpecs{m,2});
end
hold off;

nexttile(5); hold on; % Fig. 2: power spectral density
for m = 1:size(modeSpecs, 1)
    plotCrit(modeSpecs{m,1}, 0.5, 22, 'psd', modeSpecs{m,2});
end
hold off;

nexttile(6); hold on; % Fig. 2: autocorrelation function
for m = 1:size(modeSpecs, 1)
    plotCrit(modeSpecs{m,1}, 0.5, 22, 'acf', modeSpecs{m,2});
end
hold off;

sgtitle('Fig. 2 (panels 1,4,5,6) and Supp. Fig. S3 (panels 2,3) for this case');

%% ---- Figure: Fig. 3 (main text) ----
% Population-average-vs-single-mode comparison: d2 plotted against each
% empirical criticality measure, colored by Re(lambda); population
% average shown as a distinct magenta reference point.
EiG = eigvals(IdxX);
dcc    = rsnz - rsnzPredicted;
dccDS  = rsnzDS - ((tauTDS - 1) / (tauDS - 1));

figure;
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(1);
scatter(dcc, d2Env, 30, real(EiG), 'filled'); hold on;
scatter(dccDS, d2DS, 40, 'magenta', 'filled');
colorbar; colormap(turbo); caxis([-1 1]);
xlabel('DCC'); ylabel('d2');

nexttile(2);
scatter(rho, d2Env, 30, real(EiG), 'filled'); hold on;
scatter(rhoDS, d2DS, 40, 'magenta', 'filled');
colorbar; colormap(turbo); caxis([-1 1]);
xlabel('Size PL range (dec)'); ylabel('d2');

nexttile(3);
scatter(rhoT, d2Env, 30, real(EiG), 'filled'); hold on;
scatter(rhoTDS, d2DS, 40, 'magenta', 'filled');
colorbar; colormap(turbo); caxis([-1 1]);
xlabel('Duration PL range (dec)'); ylabel('d2');

nexttile(4);
scatter(alfa, d2Env, 30, real(EiG), 'filled'); hold on;
scatter(alfaDS, d2DS, 40, 'magenta', 'filled');
colorbar; colormap(turbo); caxis([-1 1]);
xlabel('PSD exponent'); ylabel('d2');

nexttile(5);
scatter(beta, d2Env, 30, real(EiG), 'filled'); hold on;
scatter(betaDS, d2DS, 40, 'magenta', 'filled');
colorbar; colormap(turbo); caxis([-1 1]);
xlabel('DFA exponent'); ylabel('d2');

nexttile(6);
scatter(corrTime, d2Env, 30, real(EiG), 'filled'); hold on;
scatter(corrTimeDS, d2DS, 40, 'magenta', 'filled');
colorbar; colormap(turbo); caxis([-1 1]);
xlabel('ACF cutoff'); ylabel('d2');

sgtitle('Fig. 3 panels (b or c, depending on case) for this case');

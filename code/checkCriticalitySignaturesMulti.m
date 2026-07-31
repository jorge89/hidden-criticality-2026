function out = checkCriticalitySignaturesMulti(signals, labels, thresholdQuantile, SF, axSize, axDFA, axPSD, axACF)
% out = checkCriticalitySignaturesMulti(signals, labels, thresholdQuantile, SF, axSize, axDFA, axPSD, axACF)
%
% signals: cell array {1xT row vectors} OR matrix [nSignals x T]
% labels:  cell array of strings (nSignals)
% thresholdQuantile: e.g. 0.8 (threshold is quantile(signal, thresholdQuantile))
% SF: sampling frequency for PSD (Hz). If your time axis is samples, set SF=1.
% ax*: axes handles for plotting.

% ---- normalize input to cell array of row vectors ----
if ~iscell(signals)
    % assume [nSignals x T]
    sigMat = signals;
    nSig = size(sigMat,1);
    signals = cell(1,nSig);
    for i=1:nSig
        signals{i} = sigMat(i,:);
    end
else
    nSig = numel(signals);
end

if nargin < 2 || isempty(labels)
    labels = arrayfun(@(i)sprintf('sig%d',i), 1:nSig, 'UniformOutput', false);
end

% ----eplace any column vectors with row vectors
for i=1:nSig
    s = signals{i};
    if size(s,1) > size(s,2)
        s = s';
    end
    signals{i} = s;
end

% ---- output struct ----
out = struct();
out.labels = labels;
out.thresholdQuantile = thresholdQuantile;

% per-signal arrays
out.tau      = nan(1,nSig);
out.tauT     = nan(1,nSig);
out.rho      = nan(1,nSig);
out.rhoT     = nan(1,nSig);
out.rsnz     = nan(1,nSig);
out.corrTime = nan(1,nSig);
out.alfa     = nan(1,nSig);
out.beta     = nan(1,nSig);

out.mm   = nan(1,nSig); out.MM   = nan(1,nSig);
out.mmT  = nan(1,nSig); out.MMT  = nan(1,nSig);

out.sizes    = cell(1,nSig);
out.durations= cell(1,nSig);

out.dfa.lengths = cell(1,nSig);
out.dfa.F       = cell(1,nSig);

out.psd.freq = cell(1,nSig);
out.psd.pxx  = cell(1,nSig);

out.acf.lags = cell(1,nSig);
out.acf.acf  = cell(1,nSig);

% ---- compute for each signal ----
gof = 0.9;
acfCutoff = 0.1;

for i=1:nSig
    sig = signals{i};

    % threshold is quantile of THIS signal (common in avalanche definitions)
    %thr = quantile(sig, thresholdQuantile);

    % avalanches
    [sizes, durations] = getAvalanches(sig, thresholdQuantile, 1);
    out.sizes{i} = sizes;
    out.durations{i} = durations;

    % PL fits
    [out.tau(i), out.rho(i), mm, MM, ~, ~, cont]     = plfit2023(sizes, gof, 0);
    [out.tauT(i),out.rhoT(i),mmT,MMT,~, ~, contT]    = plfit2023(durations, gof, 0);
    out.mm(i)=mm; out.MM(i)=MM; out.mmT(i)=mmT; out.MMT(i)=MMT;

    % rsnz via size vs duration within fitted duration range
    mask = (mmT < durations) & (durations <= MMT);
    x = durations(mask);
    y = sizes(mask);
    if numel(x) >= 5 && all(x>0) && all(y>0)
        coeff = polyfit(log10(x), log10(y), 1);
        out.rsnz(i) = coeff(1);
    end

    % ACF + corr time
    try
        [acf, lags] = autocorr(sig, 'NumLags', round(length(sig)/10));
    catch
        % fallback if autocorr not available
        maxL = round(length(sig)/10);
        [acf, lags] = xcorr(sig-mean(sig), maxL, 'coeff');
        mid = maxL+1;
        acf = acf(mid:end);
        lags = 0:maxL;
    end
    out.acf.acf{i} = acf;
    out.acf.lags{i} = lags;

    idx = find(acf < acfCutoff, 1, 'first');
    if ~isempty(idx)
        out.corrTime(i) = lags(idx);
    else
        out.corrTime(i) = 1000;
    end

    % PSD (Welch)
    [pxx, freq] = pwelch(sig-mean(sig), round(length(sig)/10), [], [], SF);
    out.psd.pxx{i} = pxx;
    out.psd.freq{i} = freq;

    % PSD slope alpha (fit top 70% of log-freq, like your earlier trimming)
    freqNZ = freq(freq>0);
    if ~isempty(freqNZ)
        log_freq = logspace(log10(min(freqNZ)), log10(max(freq)), length(freq));
        cutoff_log_freq = quantile(log_freq, 0.3);
        valid = freq >= cutoff_log_freq & freq > 0 & pxx > 0;
        if nnz(valid) >= 5
            coeffPS = polyfit(log10(freq(valid)), log10(pxx(valid)), 1);
            out.alfa(i) = coeffPS(1);
        end
    end

    % ----- DFA beta (robust to DFA output length mismatch) -----
    lengths = round(logspace(log10(4), log10(length(sig)/10), ...
                     max(5, round(10*log10(length(sig)/10/4)))));  % ensure >=5 points
    
    [F, ~] = DFA(sig, lengths);
    
    % Force column vectors
    lengths = lengths(:);
    F = F(:);
    
    % If DFA returns fewer/more points than requested, align safely
    m = min(numel(lengths), numel(F));
    lengths = lengths(1:m);
    F = F(1:m);
    
    out.dfa.lengths{i} = lengths;
    out.dfa.F{i} = F;
    
    % Build valid mask safely (same length for all terms)
    maskBasic = (lengths > 0) & (F > 0) & isfinite(lengths) & isfinite(F);
    if nnz(maskBasic) >= 5
        log_lengths = log10(lengths(maskBasic));
        cutoff_log_length = quantile(log_lengths, 0.7);
    
        valid = maskBasic & (log10(lengths) <= cutoff_log_length);
    
        if nnz(valid) >= 5
            coeffDFA = polyfit(log10(lengths(valid)), log10(F(valid)), 1);
            out.beta(i) = coeffDFA(1);
        end
    end

end

% =========================
% PLOTTING (4 square panels)
% =========================
% colors (MATLAB default lines)
cols = lines(nSig);

% ---- (1) Size distributions with vertical shifts ----
axes(axSize); cla(axSize); hold(axSize,'on');
shiftDecades = 0.6;
cols = lines(nSig);

for i=1:nSig
    s = out.sizes{i};
    s = s(s>0 & isfinite(s));
    if numel(s) < 20, continue; end

    nbins = 25;
    edges = logspace(log10(min(s)), log10(max(s)), nbins+1);
    counts = histcounts(s, edges, 'Normalization','pdf');
    xc = sqrt(edges(1:end-1).*edges(2:end));
    y = counts(:)';

    ok = (xc>0) & (y>0) & isfinite(xc) & isfinite(y);
    xc = xc(ok);
    y  = y(ok) * 10^(shiftDecades*(i-1));

    plot(axSize, xc, y, '-', 'Color', cols(i,:), 'LineWidth', 1.3);
end

set(axSize,'XScale','log','YScale','log');
xlabel(axSize,'Avalanche size S');
ylabel(axSize,'PDF (shifted)');
title(axSize,'Size distributions (shifted)');
grid(axSize,'on'); box(axSize,'off');
legend(axSize, labels, 'Location','southwest');


% ---- (2) DFA ----
axes(axDFA); cla(axDFA); hold(axDFA,'on');
cols = lines(nSig);

for i=1:nSig
    L = out.dfa.lengths{i};  L = L(:);
    F = out.dfa.F{i};        F = F(:);

    m = min(numel(L), numel(F));
    L = L(1:m); F = F(1:m);

    ok = (L>0) & (F>0) & isfinite(L) & isfinite(F);
    if nnz(ok) < 5, continue; end

    plot(axDFA, L(ok), F(ok), '-', 'Color', cols(i,:), 'LineWidth', 1.3);
end

set(axDFA,'XScale','log','YScale','log');
xlabel(axDFA,'Window length w');
ylabel(axDFA,'F(w)');
title(axDFA,'DFA');
grid(axDFA,'on'); box(axDFA,'off');

% ---- (3) PSD ----
axes(axPSD); cla(axPSD); hold(axPSD,'on');
cols = lines(nSig);

for i=1:nSig
    f = out.psd.freq{i}(:);
    p = out.psd.pxx{i}(:);

    m = min(numel(f), numel(p));
    f = f(1:m); p = p(1:m);

    ok = (f>0) & (p>0) & isfinite(f) & isfinite(p);
    if nnz(ok) < 5, continue; end

    plot(axPSD, f(ok), p(ok), '-', 'Color', cols(i,:), 'LineWidth', 1.3);
end

set(axPSD,'XScale','log','YScale','log');
xlabel(axPSD,'Frequency (Hz)');
ylabel(axPSD,'PSD');
title(axPSD,'PSD (Welch)');
grid(axPSD,'on'); box(axPSD,'off');

% ---- (4) Autocorrelation ----
axes(axACF); cla(axACF); hold(axACF,'on');
for i=1:nSig
    l = out.acf.lags{i};
    a = out.acf.acf{i};
    if numel(l) < 5, continue; end
    plot(axACF, l, a, '-', 'Color', cols(i,:), 'LineWidth', 1.3);
end
yline(axACF, acfCutoff, '--k', 'LineWidth', 1.0);
set(axACF,'XScale','log');
xlabel(axACF,'Lag');
ylabel(axACF,'ACF');
title(axACF,'Autocorrelation');
grid(axACF,'on'); box(axACF,'off');

end

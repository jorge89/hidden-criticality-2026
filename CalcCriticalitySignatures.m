function [rho, rhoT, tau,predS, tauT,predT, rsnz,rsnzPredicted,alfa,beta, corrTime] = CalcCriticalitySignatures(signal, threshold,SF)
%checkCriticalitySignatures takes a time series (signal) and the avalanche
%threshold (threshold, specified as a quantile of signal) and returns 
%the avalanche size distribution power-law range (rho), the duration
%distribution power-law range (rhoT), the size (tau)
%and duration (tauT) exponents, the exponent 1/(sigma*nu*z) for size vs.
%duration (rsnz), and the autocorrelation time (corrTime), defined as the
%lag at which the autocorrelation first dips below 0.1.

if size(signal, 1) > size(signal, 2)
    signal = signal';
end

gof = 0.95;%goodness of fit criterion for avalanche distribution power-law fitting
[sizes, durations] = getAvalanches(signal, threshold, 1);%get avalanche sizes and durations
[tau, rho, mm, MM, nsamp, bestfracgood, continuous] = plfit2023(sizes, gof, 0);%fit size distribution to a power-law
[tauT, rhoT, mmT, MMT, nsampT, bestfracgoodT, continuousT] = plfit2023(durations, gof, 0);%fit duration distribution to a power-law

%calculate 1/(sigma*nu*z)
x = durations(mmT<durations & durations<=MMT);
y = sizes(mmT<durations & durations<=MMT);
coeff = polyfit(log10(x), log10(y), 1);
rsnz = coeff(1);
rsnzPredicted = (tauT-1)/(tau-1);
predT = 1+(tauT-1)/rsnz;
predS = 1+(tau-1)*rsnz;
%calculate autocor
acfCutoff = 0.1;
[acf, lags] = autocorr(signal, 'NumLags', round(length(signal)/10));%autocorrelation function
if length(lags(acf<acfCutoff))>0
    corrTime = min(lags(acf<acfCutoff));%autocorrelation time
else
    corrTime = 1000;
    disp("Warning: autocorrelation never falls below cutoff.")
end
% power spectrum
[pxx, freq] = pwelch(signal-mean(signal), round(length(signal)/10), [], [], SF);%power spectral density with Welch's method
% Exclude lowest 20% of the frequency points
log_freq = logspace(log10(min(freq(freq > 0))), log10(max(freq)), length(freq));
cutoff_log_freq = quantile(log_freq, 0.3);
valid_idx_freq = freq >= cutoff_log_freq;

% Create reduced vectors
xx = freq(valid_idx_freq);
yy = pxx(valid_idx_freq);
% Fit in log-log space
coeffPS = polyfit(log10(xx), log10(yy), 1);
alfa = coeffPS(1);
%detrended fluctuation analysis (DFA) with 50% overlap between adjacent windows
% Determine max X (scale) and cutoff threshold (70% of that)
lengths = round(logspace(log10(4), log10(length(signal)/10), round(10*log10(length(signal)/10/4))));
[F, numBins] = DFA(signal, lengths);
log_lengths = log10(lengths);

% Determine cutoff in log space (exclude top 30%)
cutoff_log_length = quantile(log_lengths, 0.7);  % or: max(log_lengths) * 0.8 if linear threshold in log

% Find indices of values below the cutoff
valid_idx = log_lengths <= cutoff_log_length;  

% Apply trimming
xxx = lengths(valid_idx);
yyy = F(valid_idx);

% Fit on trimmed log-log data
coeffDFA = polyfit(log10(xxx), log10(yyy), 1);
beta = coeffDFA(1);
end
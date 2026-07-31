function [] = plotCrit(signal,threshold,SF,plotSelection,co,shift)
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


%fit duration distribution to a power-law




  





if ischar(plotSelection)
    plotSelection = {plotSelection}; % Ensure it's a cell array
end

%figure;

if any(strcmp(plotSelection, 'sizeDist'))
    [sizes, durations] = getAvalanches(signal, threshold, 1);%get avalanche sizes and durations
    [tau, rho, mm, MM, nsamp, bestfracgood, continuous] = plfit2023(sizes, gof, 0);%fit size distribution to a power-law
    [tauT, rhoT, mmT, MMT, nsampT, bestfracgoodT, continuousT] = plfit2023(durations, gof, 0);
    x = durations(mmT<durations & durations<=MMT);
    y = sizes(mmT<durations & durations<=MMT);
    coeff = polyfit(log10(x), log10(y), 1);
    rsnz = coeff(1);
    hold on;
    loglogDist(sizes, co, 'o', 3, false, false, mm, MM, shift, true);
    pred = 1 + (tauT - 1) / rsnz;
    if continuous
        normfac = 1 / (1.00001 - tau) * (MM^(1.00001 - tau) - mm^(1.00001 - tau));
        normfac2 = 1 / (1.00001 - pred) * (MM^(1.00001 - pred) -mm^(1.00001 - pred));
    else
        xx = mm:MM;
        normfac = sum(xx.^(-tau));
        normfac2 = sum(xx.^(-pred));
    end
    hold on;
    x = logspace(log10(mm), log10(MM), 100);
    plot(x, shift*x.^(-tau) / normfac, '--r', 'LineWidth', 1.5);
    plot(x, shift*x.^(-pred) / normfac2, '--g', 'LineWidth', 1.5);
    xlabel('Avalanche Size S')
    ylabel('Probability Density P(S)')
end

if any(strcmp(plotSelection, 'durationDist'))
    [sizes, durations] = getAvalanches(signal, threshold, 1);%get avalanche sizes and durations
    [tau, rho, mm, MM, nsamp, bestfracgood, continuous] = plfit2023(sizes, gof, 0);%fit size distribution to a power-law
    [tauT, rhoT, mmT, MMT, nsampT, bestfracgoodT, continuousT] = plfit2023(durations, gof, 0);
    x = durations(mmT<durations & durations<=MMT);
    y = sizes(mmT<durations & durations<=MMT);
    coeff = polyfit(log10(x), log10(y), 1);
    rsnz = coeff(1);
    hold on;
    loglogDist(durations, co, 'o', 3, false, false, mmT, MMT, shift, true);
    pred = 1 + (tau - 1) * rsnz;    
    if continuousT
        normfac = 1 / (1.00001 - tauT) * (MMT^(1.00001 - tauT) - mmT^(1.00001 - tauT));
        normfac2 = 1 / (1.00001 - pred) * (MMT^(1.00001 - pred) - mmT^(1.00001 - pred));
    else
        xx = mmT:MMT;
        normfac = sum(xx.^(-tauT));
        normfac2 = sum(xx.^(-pred));
    end

    x = logspace(log10(mmT), log10(MMT), 100);
    plot(x, shift*x.^(-tauT) / normfac, '--r', 'LineWidth', 1.5);
    plot(x, shift*x.^(-pred) / normfac2, 'g', 'LineWidth', 1.5);
    xlabel('Avalanche Duration T')
    ylabel('Probability Density P(T)')
end

if any(strcmp(plotSelection, 'sizeVsDuration'))
    [sizes, durations] = getAvalanches(signal, threshold, 1);%get avalanche sizes and durations
    [tau, rho, mm, MM, nsamp, bestfracgood, continuous] = plfit2023(sizes, gof, 0);%fit size distribution to a power-law
    [tauT, rhoT, mmT, MMT, nsampT, bestfracgoodT, continuousT] = plfit2023(durations, gof, 0);
    x = durations(mmT<durations & durations<=MMT);
    y = sizes(mmT<durations & durations<=MMT);
    coeff = polyfit(log10(x), log10(y), 1);
    rsnz = coeff(1);

    %    subplot(2, 3, 3);
    hold on;
    scatter(durations, shift*sizes,2,co)
        %co,'LineWidth', 0.001, 'LineStyle', 'none','Marker','o', 'MarkerSize', 2);
    x = logspace(log10(mmT), log10(MMT), 100);
    plot(x, shift*10^coeff(2)*x.^rsnz, 'r', 'LineWidth', 1.5);
    rsnzPredicted = (tauT-1)/(tau-1);
    plot(x, shift*10^coeff(2)*x.^rsnzPredicted, 'g', 'LineWidth', 1.5);
    xlabel('Avalanche Duration T');
    ylabel('Avalanche Size S')
    set(gca, 'xscale', 'log', 'yscale', 'log')
end

if any(strcmp(plotSelection, 'dfa'))
    %detrended fluctuation analysis (DFA) with 50% overlap between adjacent windows
    lengths = round(logspace(log10(4), log10(length(signal)/10), round(10*log10(length(signal)/10/4))));
    [F, numBins] = DFA(signal, lengths);
%    subplot(2, 3, 4);
    hold on;
    plot(lengths, F, co, 'LineWidth', 1.5);
    set(gca, 'xscale', 'log', 'yscale', 'log')
    xlabel('Window Length w (time step)')
    ylabel('Fluctuation Function F(w)')
end

if any(strcmp(plotSelection, 'psd'))
    [pxx, freq] = pwelch(signal-mean(signal), round(length(signal)/10), [], [], SF);%power spectral density with Welch's method
%    subplot(2, 3, 5);
    hold on;
    plot(freq, pxx, co, 'LineWidth', 1.5);
    set(gca, 'xscale', 'log', 'yscale', 'log')
    xlabel('Frequency (time step)^{-1}')
    ylabel('Power Spectral Density')
end

if any(strcmp(plotSelection, 'acf'))
    acfCutoff = 0.1;
    [acf, lags] = autocorr(signal, 'NumLags', round(length(signal)/10));%autocorrelation function
    
    if length(lags(acf<acfCutoff))>0
        corrTime = min(lags(acf<acfCutoff));%autocorrelation time
    else
        corrTime = 1000;
        disp("Warning: autocorrelation never falls below cutoff.")
    end
%    subplot(2, 3, 6);
    plot(lags, acf, co, 'LineWidth', 1.5);
    yline(acfCutoff, '--b', 'LineWidth', 1.5);
    xlabel('Lag (time steps)')
    ylabel('Autocorrelation')
    set(gca, 'xscale', 'log')
    box off;
end



end
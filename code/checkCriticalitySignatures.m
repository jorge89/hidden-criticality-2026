function [rho, rhoT, tau, tauT, rsnz, corrTime, alfa, beta] = checkCriticalitySignatures(signal, threshold,SF,ploT)
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
display(tau)
display(tauT)
%calculate 1/(sigma*nu*z)
x = durations(mmT<durations & durations<=MMT);
y = sizes(mmT<durations & durations<=MMT);
coeff = polyfit(log10(x), log10(y), 1);
rsnz = coeff(1);

acfCutoff = 0.1;
[acf, lags] = autocorr(signal, 'NumLags', round(length(signal)/10));%autocorrelation function
if length(lags(acf<acfCutoff))>0
    corrTime = min(lags(acf<acfCutoff));%autocorrelation time
else
    corrTime = 1000;
    disp("Warning: autocorrelation never falls below cutoff.")
end
if ploT==1  
    [pxx, freq] = pwelch(signal-mean(signal), round(length(signal)/10), [], [], SF);%power spectral density with Welch's method
    
    
    log_freq = logspace(log10(min(freq(freq > 0))), log10(max(freq)), length(freq));
    cutoff_log_freq = quantile(log_freq, 0.3);
    valid_idx_freq = freq >= cutoff_log_freq;

    % Create reduced vectors
    xxxx = freq(valid_idx_freq);
    yyyy = pxx(valid_idx_freq);
    % Fit in log-log space
    coeffPS = polyfit(log10(xxxx), log10(yyyy), 1);
    alfa = coeffPS(1)  
    %detrended fluctuation analysis (DFA) with 50% overlap between adjacent windows
    lengths = round(logspace(log10(4), log10(length(signal)/10), round(10*log10(length(signal)/10/4))));
    [F, numBins] = DFA(signal, lengths);
    % Determine max X (scale) and cutoff threshold (80% of that)
    log_lengths = log10(lengths);
    
    % Determine cutoff in log space (exclude top 20%)
    cutoff_log_length = quantile(log_lengths, 0.7);  % or: max(log_lengths) * 0.8 if linear threshold in log
    
    % Find indices of values below the cutoff
    valid_idx = log_lengths <= cutoff_log_length;    
    % Apply trimming
    xxx = lengths(valid_idx);
    yyy = F(valid_idx);

% Fit on trimmed log-log data
    coeffDFA = polyfit(log10(xxx), log10(yyy), 1);
    beta = coeffDFA(1)
    
    figure;
    subplot(2, 3, 1);
    hold on;
    loglogDist(sizes, 'b', 'o', 6, false, false, mm, MM, 1, true);
    xlabel('Avalanche Size S')
    ylabel('Probability Density P(S)')
    %xlim([1 40000]);        % Set x-axis limits from 0 to 10
    %ylim([0.00000005 1]);   
    subplot(2, 3, 2);
    hold on;
    loglogDist(durations, 'b', 'o', 6, false, false, mmT, MMT, 1, true);
    xlabel('Avalanche Duration T')
    ylabel('Probability Density P(T)');
    %xlim([1 1000]);        % Set x-axis limits from 0 to 10
    %ylim([0.0000005 2]);  
    subplot(2, 3, 3);
    hold on;
    plot(durations, sizes, '.b', 'MarkerSize', 6);
    x = logspace(log10(min(durations)), log10(max(durations)), 100);
    plot(x, 10^coeff(2)*x.^rsnz, '--r', 'LineWidth', 1.5);
    rsnzPredicted = (tauT-1)/(tau-1); %1/(sigma*nu*z) exponent predicted by the crackling noise scaling relation
    plot(x, 10^coeff(2)*x.^(rsnzPredicted), 'g', 'LineWidth', 1.5);
    xlabel('Avalanche Duration T');
    ylabel('Avalanche Size S')
    %xlim([1 1000]);
    %ylim([1 40000]);
    %legend('', '$\frac{1}{\sigma\nu z}$', '$\frac{\tau_t-1}{\tau-1}$', 'Interpreter', 'latex', 'FontSize', 12)
    set(gca, 'xscale', 'log', 'yscale', 'log')
    
    subplot(2, 3, 4);
    hold on;
    plot(lengths, F, 'b', 'LineWidth', 1.5);
    plot(xxx, 10^coeffDFA(2)*xxx.^(coeffDFA(1)), 'g', 'LineWidth', 1.5);
    set(gca, 'xscale', 'log', 'yscale', 'log')
    xlabel('Window Length w (time step)')
    ylabel('Fluctuation Function F(w)')
     % Set y-axis limits from -1.5 to 1.5
    subplot(2, 3, 5);
    hold on;
    plot(freq, pxx, '*', 'LineWidth', 1.5);%10^coeffPS(2)*xxxx.^(coeffPS(1))
    plot(xxxx, 10.^(polyval(coeffPS, log10(xxxx))), 'g', 'LineWidth', 1.5);
    set(gca, 'xscale', 'log', 'yscale', 'log')
    xlabel('Frequency (time step)^{-1}')
    ylabel('Power Spectral Density')
    
    subplot(2, 3, 6);
    plot(lags, acf, 'b', 'LineWidth', 1.5);
    yline(acfCutoff, '--b', 'LineWidth', 1.5);
    xlabel('Lag (time steps)')
    ylabel('Autocorrelation')
    set(gca, 'xscale', 'log')
    box off;
    
    %plot binned pdfs of best-fit power-laws
    for whichOne = [1, 2]
        if whichOne == 1
            temp = sizes;
            expo = tau;
            pred = 1+(tauT-1)/rsnz;
            cont = continuous;
            subplot(2, 3, 1);
            hold on;
        else
            temp = durations;
            expo = tauT;
            cont = continuousT;
            pred = 1+(tau-1)*rsnz;
            subplot(2, 3, 2);
            hold on;
        end
    
        if cont
            normfac = 1/(1-expo)*(max(temp)^(1-expo) - min(temp)^(1-expo));
            normfac2 = 1/(1-pred)*(max(temp)^(1-pred) - min(temp)^(1-pred));
        else
            xx = min(temp):max(temp);
            normfac = sum(xx.^(-expo));
            normfac2 = sum(xx.^(-pred));
        end
        
        x = logspace(log10(min(temp)), log10(max(temp)), 100);
        h = zeros(1, 2);
        h(1) = plot(x, x.^(-expo)/normfac, '--r', 'LineWidth', 1.5);
        h(2) = plot(x, x.^(-pred)/normfac2, 'g', 'LineWidth', 1.5);
    
        if whichOne == 2
            legend(h, 'best-fit', 'predicted');
        end
    end
    
    
    set(gcf, 'Color', 'w')
end

end
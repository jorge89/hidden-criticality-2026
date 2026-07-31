function [X_binned, bin_edges] = bin_timeseries(X, bin_size)
    % BIN_TIMESERIES bins the time series data in X over time windows
    % 
    % Inputs:
    %   X         - N x T matrix of time series data (neurons x time)
    %   dt        - Time step of simulation (e.g., 0.01)
    %   bin_width - Width of bin in time units (e.g., 10)
    %
    % Outputs:
    %   X_binned  - N x num_bins matrix of binned data (mean over bin)
    %   bin_edges - Time vector indicating the start of each bin

    [N, T] = size(X);
    bin_width = 23;%round(bin_width / dt); % number of time steps per bin
    num_bins = floor(T / bin_size);

    X_binned = zeros(N, num_bins);
    bin_edges = (0:num_bins-1) * bin_width;

    for b = 1:num_bins
        idx = (1:bin_size) + (b-1)*bin_size;
        X_binned(:, b) = mean(X(:, idx), 2);
    end
end

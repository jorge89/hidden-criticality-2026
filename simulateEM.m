function XX = simulateEM(n,T,sigma,A,iniCond,dt,tau)
% Step 6: Simulate autoregressive dynamics
    X = zeros(n, T);
    X(:,1) = iniCond; % Initial condition
    I = eye(n);
    update_matrix = (A - I);
    noise_t = randn(n,T);
    % Simulation loop
%     dt = 0.002;              % time step
%     tau = 0.05;  
    bin_width = 23;% time constant
    % Euler–Maruyama integration
    for t = 1:T-1
        drift = (dt / tau) * update_matrix * X(:,t);
        diffusion = (sigma) * sqrt(dt/ tau) * noise_t(:,t);
        X(:,t+1) = X(:,t) + drift + diffusion;
    end
    %XX = X;%no corse grain applied
    [XX,~] = bin_timeseries(X, bin_width);
end
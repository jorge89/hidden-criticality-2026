function J = build_IE_matrix(N, alpha, p, w, g, seed)
%BUILD_IE_MATRIX_WITH_SEED Constructs synaptic connectivity matrix J with seed control
%
%   J = BUILD_IE_MATRIX_WITH_SEED(N, alpha, p, w, g, seed)
%
%   Inputs:
%       N     - Total number of neurons
%       alpha - Fraction of inhibitory neurons (e.g., 0.2)
%       p     - Connection probability
%       w     - Maximum excitatory synaptic weight
%       g     - Inhibitory to excitatory weight ratio
%       seed  - Seed for the random number generator
%
%   Output:
%       J     - NxN connectivity matrix (row = target, col = source)

    if nargin == 6
        rng(seed);  % Set random seed
    end

    % Number of inhibitory and excitatory neurons
    N_I = round(alpha * N);
    N_E = N - N_I;

    % Initialize matrix
    J = zeros(N, N);

    % Assign excitatory synapses
    for i = 1:N_E
        targets = rand(N, 1) < p;
        J(targets, i) = rand(sum(targets), 1) * w;
    end

    % Assign inhibitory synapses
    for i = N_E+1:N
        targets = rand(N, 1) < p;
        J(targets, i) = -rand(sum(targets), 1) * g * w;
    end
end

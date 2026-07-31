function W = make_winner_take_all_connectivity(N_e, N_i, p_dense, p_sparse, normalize)
%MAKE_WINNER_TAKE_ALL_CONNECTIVITY Generate connectivity matrix from Figure 8A of Jones et al. (2023)
%
%   W = MAKE_WINNER_TAKE_ALL_CONNECTIVITY(N_e, N_i, p_dense, p_sparse, normalize)
%
%   INPUTS:
%     N_e       - number of excitatory neurons per group (e+ and e-)
%     N_i       - number of inhibitory neurons per group (i+ and i-)
%     p_dense   - connection probability for dense connections (e.g., 0.5)
%     p_sparse  - connection probability for sparse connections (e.g., 0.05)
%     normalize - true/false: normalize W so that max eigenvalue = 1
%
%   OUTPUT:
%     W         - full (N_total x N_total) connectivity matrix

    % Total number of neurons
    N_total = 2 * (N_e + N_i);
    W = zeros(N_total);

    % Neuron group indices
    idx.eplus  = 1:N_e;
    idx.eminus = N_e + (1:N_e);
    idx.iplus  = 2*N_e + (1:N_i);
    idx.iminus = 2*N_e + N_i + (1:N_i);

    % Weights
    w_E = 1;   % excitatory
    w_I = -1;  % inhibitory

    % --- Excitatory connections ---
    % Within-group dense
    W(idx.eplus,  idx.eplus)  = rand(N_e) < p_dense;
    W(idx.eminus, idx.eminus) = rand(N_e) < p_dense;

    % Across-group sparse
    W(idx.eplus,  idx.eminus) = (rand(N_e) < p_sparse);
    W(idx.eminus, idx.eplus)  = (rand(N_e) < p_sparse);

    % --- Crossing inhibition (e -> i) ---
    W(idx.iminus, idx.eplus)  = rand(N_i, N_e) < p_dense;
    W(idx.iplus,  idx.eminus) = rand(N_i, N_e) < p_dense;

    % --- Feedback inhibition (i -> e) ---
    W(idx.eplus,  idx.iplus)  = rand(N_e, N_i) < p_dense;
    W(idx.eminus, idx.iminus) = rand(N_e, N_i) < p_dense;

    % --- Inhibitory self-group (optional) ---
    W(idx.iplus,  idx.iplus)  = rand(N_i) < p_dense;
    W(idx.iminus, idx.iminus) = rand(N_i) < p_dense;

    % Apply weights: E = +1, I = -1
    W(:, [idx.iplus, idx.iminus]) = W(:, [idx.iplus, idx.iminus]) * w_I;
    W(:, [idx.eplus, idx.eminus]) = W(:, [idx.eplus, idx.eminus]) * w_E;

    % Optional normalization
    if normalize
        lambda_max = max(real(eig(W)));
        W = W * (0.999 / lambda_max);
        %W = W / lambda_max;
    end
end

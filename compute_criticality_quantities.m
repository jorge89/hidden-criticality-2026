%% compute_criticality_quantities.m
%
% Computes the quantities plotted by makeFigs.m for one of the three
% cases described in the paper's Methods:
%   Case 1 — non-hidden, non-oscillatory criticality (build_IE_matrix)
%   Case 2 — hidden criticality with anticorrelations (make_winner_take_all_connectivity)
%   Case 3 — hidden criticality with mixed/oscillatory modes (saved dense random matrix)
%
% Requires (place in the same folder or on the MATLAB path):
%   build_IE_matrix.m
%   make_winner_take_all_connectivity.m
%   simulateEM.m
%   CalcCriticalitySignatures.m
%   checkCriticalitySignatures.m
%   getD2v2.m
%
% Output variables used downstream by makeFigs.m:
%   eigvals, IdxX, Idx, dataSum, dataModes, dataModesEnv,
%   rho, rhoT, alfa, beta, corrTime, d2Env,
%   rhoDS, rhoTDS, tauDS, tauTDS, rsnzDS, corrTimeDS, alfaDS, betaDS

%% ---- 1. Choose case and set parameters ----
caseID = 3;         % 1, 2, or 3 — selects which interaction matrix to build
n      = 500;       % number of neurons
T      = 200000;    % number of simulated time steps per ensemble run
sigma  = 1;         % noise level
dt     = 0.002;     % simulation time step
tau    = 0.05;      % simulation time constant
n_runs = 20;         % number of ensemble repetitions

%% ---- 2. Construct the interaction matrix A ----
switch caseID
    case 1
        % Two excitatory/inhibitory populations with reciprocal coupling
        A = build_IE_matrix(n, 0.2, 0.2, 0.016, 1.02, 1);
    case 2
        % Two excitatory groups competing via crossing inhibition
        A = make_winner_take_all_connectivity(200, 50, 0.5, 0.05, 1);
    case 3
        % Dense random matrix (Gaussian-distributed weights), precomputed
        % and stored in saved_matrices.mat — load it here.
        load('saved_matrices.mat', 'saved_matrices');
        A = saved_matrices{2};
    otherwise
        error('caseID must be 1, 2, or 3.');
end

% Normalize so the real part of the largest eigenvalue sits at the
% target value used throughout the paper (near-critical).
% NOTE: paper Methods states 0.999 for all three cases; this script
% historically used 0.998 — confirm which value was actually used to
% generate the published results before finalizing.
target_re_lambda = 0.999;
eigsA = eig(A);
max_eig = max(real(eigsA));
if max_eig > 0
    A = A * (target_re_lambda / max_eig);
end

%% ---- 3. Eigendecomposition and biorthogonal normalization ----
[Vuns, Duns, Wuns] = eig(A);
eigvals = diag(Duns);
[~, ind] = sort(real(eigvals));
eigvals  = eigvals(ind);
V = Vuns(:, ind);
W = Wuns(:, ind);

% Normalize so that W' * V = I (biorthogonal system)
for i = 1:n
    norm_factor = W(:, i)' * V(:, i);
    V(:, i) = V(:, i) / norm_factor;
end

%% ---- 4. Build index list over unique modes (skip complex conjugates) ----
x_int = 1:n;
ii = 1;
IdxX = [];
while ii <= length(x_int)
    i = x_int(ii);
    IdxX = [IdxX i]; %#ok<AGROW>
    if ~isreal(eigvals(i)) && ii < length(x_int)
        ii = ii + 2;  % skip the complex-conjugate partner
    else
        ii = ii + 1;
    end
end

%% ---- 5. Ensemble simulation ----
iniCond   = randn(n, 1);
dataModes = [];
dataSum   = [];

for jj = 1:n_runs
    disp(jj)
    MA = simulateEM(n, T, sigma, A, iniCond, dt, tau);

    % Project activity onto each unique mode
    coeffs = W' * MA;
    relevant_proj = [];
    Idx = [];
    ii = 1;
    while ii <= length(x_int)
        i = x_int(ii);
        actSum = coeffs(i, :);
        relevant_proj = [relevant_proj; actSum]; %#ok<AGROW>
        if ~isreal(eigvals(i)) && ii < length(x_int)
            Idx = [Idx; 1]; %#ok<AGROW>
            ii = ii + 2;
        else
            Idx = [Idx; 0]; %#ok<AGROW>
            ii = ii + 1;
        end
    end

    dataSum   = [dataSum sum(MA, 1)]; %#ok<AGROW>
    dataModes = [dataModes relevant_proj]; %#ok<AGROW>
    iniCond   = MA(:, end);
end

%% ---- 6. Envelope for oscillatory (complex) modes ----
dataModesEnv = zeros(size(dataModes));
for i = 1:length(Idx)
    if Idx(i) == 1
        dataModesEnv(i, :) = abs(dataModes(i, :));   % amplitude envelope
    else
        dataModesEnv(i, :) = dataModes(i, :);         % real mode, unchanged
    end
end

%% ---- 7. Per-mode criticality signatures ----
results(size(dataModesEnv, 1)) = struct();
for i = 1:size(dataModesEnv, 1)
    disp(i)
    [rho_i, rhoT_i, tau_i, predS_i, tauT_i, predT_i, rsnz_i, rsnzPredicted_i, alfa_i, beta_i, corrTime_i] = ...
        CalcCriticalitySignatures(dataModesEnv(i, :), 0.5, 22);

    results(i).rho           = rho_i;
    results(i).rhoT          = rhoT_i;
    results(i).tau           = tau_i;
    results(i).predS         = predS_i;
    results(i).tauT          = tauT_i;
    results(i).predT         = predT_i;
    results(i).rsnz          = rsnz_i;
    results(i).rsnzPredicted = rsnzPredicted_i;
    results(i).alfa          = alfa_i;
    results(i).beta          = beta_i;
    results(i).corrTime      = corrTime_i;
end

rho           = arrayfun(@(x) x.rho, results);
rhoT          = arrayfun(@(x) x.rhoT, results);
tau           = arrayfun(@(x) x.tau, results); %#ok<NASGU>
predS         = arrayfun(@(x) x.predS, results); %#ok<NASGU>
tauT          = arrayfun(@(x) x.tauT, results); %#ok<NASGU>
predT         = arrayfun(@(x) x.predT, results); %#ok<NASGU>
rsnz          = arrayfun(@(x) x.rsnz, results); %#ok<NASGU>
rsnzPredicted = arrayfun(@(x) x.rsnzPredicted, results); %#ok<NASGU>
alfa          = arrayfun(@(x) x.alfa, results);
beta          = arrayfun(@(x) x.beta, results);
corrTime      = arrayfun(@(x) x.corrTime, results);

%% ---- 8. d2 per mode ----
d2Env = [];
for i = 1:size(dataModes, 1)
    disp(i)
    d2Env = [d2Env getD2v2(dataModesEnv(i, :), 20, 0.0435)]; %#ok<AGROW>
end

%% ---- 9. Reference values from the population sum ----
[rhoDS, rhoTDS, tauDS, tauTDS, rsnzDS, corrTimeDS, alfaDS, betaDS] = ...
    checkCriticalitySignatures(dataSum, 0.5, 22, 1);

% d2 for the population sum (same call used per-mode in Section 8),
% used as the magenta reference point in Fig. 3.
d2DS = getD2v2(dataSum, 20, 0.0435);

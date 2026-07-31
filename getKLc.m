function [KLc, closest_model] = getKLc(phi, sigma)

options = struct; options.verbosity = 0;
problem = struct;
warning('off', 'manopt:getHessian:approx');
warning('off', 'manopt:getGradient:approx')

problem.M = euclideanfactory(length(phi)-1);
problem.cost = @(phic) costKLc(phic, phi, sigma);
if length(phi)>25
    cme = [1 zeros(1, length(phi)-1)]'; 
else
    [~, cme] = getFixedPointDistance(length(phi), 2, phi);
end
[x, xcost, info, ~] = rlbfgs(problem, cme(1:(end-1)), options);
KLc = xcost;
closest_model = struct;
x = x(:);
closest_model.phic = [x ; 1-sum(x)];
closest_model.sigmac = getSigmaC(phi, sigma, closest_model.phic);

end
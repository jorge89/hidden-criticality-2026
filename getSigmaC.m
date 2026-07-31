function [sigmaC] = getSigmaC(phi, sigma, phic)
gc = @(omega) get_g_no_sigma(phic, omega);
g = @(omega) get_g_no_sigma(phi, omega);
fun = @(omega)sigma^2/(2*pi)*gc(omega)./g(omega);
sigmaC = sqrt(integral(fun, -pi, pi));
end
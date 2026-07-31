function [cost] = costKLc(phic, phi, sigma)
phic = phic(:);
phic = [phic ; 1-sum(phic)];
sigmac = getSigmaC(phi, sigma, phic);
fun = @(omega) 1/(4*pi)*(log(get_g_no_sigma(phi, omega).*sigmac.^2./get_g_no_sigma(phic, omega)./sigma.^2) + get_g_no_sigma(phic, omega).*sigma.^2./get_g_no_sigma(phi, omega)./sigmac.^2-1);
cost = integral(fun, -pi, pi);
end
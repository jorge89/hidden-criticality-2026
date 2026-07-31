function [g] = get_g_no_sigma(phi, omega)
phi_tilde = 0;
for t = 1:length(phi)
    phi_tilde = phi_tilde + phi(t)*exp(-1i*omega*t);
end
g = abs(1-phi_tilde).^2;
end
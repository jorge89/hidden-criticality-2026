% example: calculate d2 for AR(1) models
phi_list = linspace(0, 1, 100);
sigma = 1;
d2 = zeros(size(phi_list));
for i = 1:length(phi_list)
    % first argument: history kernel, second argument: noise standard deviation
    % remark: d2 is independent of the noise standard deviation
    [d2(i), ~] = getKLc(phi_list(i), sigma);
end

figure;
plot(phi_list, d2, '.k', 'MarkerSize', 12);
xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$d_2$', 'Interpreter', 'latex', 'FontSize', 12);
hmax_vector = linspace(10 * 10^-6, 18 * 10^-6, 101);
TPP = zeros(1, numel(hmax_vector));
UL = zeros(1, numel(hmax_vector));

% Parameters
R = 330 * 10^-6;
t = 0.1 * 10^-6;
E = 170 * 10^9;
v = 0.28;

for i = 1 : numel(hmax_vector)
    % Taking g = hmax_vector / 2 + (1 * 10^-6)
    g = hmax_vector(i) / 2 + (1 * 10^-6);

    % Calculating rs
    rs = hmax_vector(i) / 4 + R^2 / hmax_vector(i);
    
    % Calculating D
    D = E * (hmax_vector(i) / 2)^3 / (12 * (1 - v^2));

    % Calculating TPP
    TPP(i) = 64 * D * (g - t - rs + sqrt(rs^2 - R^2)) / R^4;

    % Calculating UL
    UL(i) = 40000 * D * (g - t - sqrt(rs^2 - 0.64 * R^2) + sqrt(rs^2 - R^2)) / (81 * R^4);
end

% Plotting the graph
plot(hmax_vector, TPP, 'b', 'LineWidth', 1.5);
hold on
plot(hmax_vector, UL, 'r', 'LineWidth', 1.5);
xlim([10 * 10^-6, 18 * 10^-6]);
xlabel('Max. Thickness (in m)');
ylabel('Working Range (in Pa)');
title('Working Range Vs Max. Thickness');
legend('TPP', 'Upper limit');
grid on;

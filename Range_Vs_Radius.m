R_vector = linspace(275 * 10^-6, 450 * 10^-6, 101);
TPP = zeros(1, numel(R_vector));
UL = zeros(1, numel(R_vector));

% Parameters
hmax = 15 * 10^-6;
g = 8.5 * 10^-6;
t = 0.1 * 10^-6;
E = 170 * 10^9;
v = 0.28;

for i = 1 : numel(R_vector)
    R = R_vector(i);

    % Calculating rs
    rs = hmax / 4 + R^2 / hmax;
    
    % Calculating D
    D = E * (hmax / 2)^3 / (12 * (1 - v^2));

    % Calculating TPP
    TPP(i) = 64 * D * (g - t - rs + sqrt(rs^2 - R^2)) / R^4;

    % Calculating UL
    UL(i) = 40000 * D * (g - t - sqrt(rs^2 - 0.64 * R^2) + sqrt(rs^2 - R^2)) / (81 * R^4);
end

% Plotting the graph
plot(R_vector, TPP, 'b', 'LineWidth', 1.5);
hold on
plot(R_vector, UL, 'r', 'LineWidth', 1.5);
xlim([275 * 10^-6, 450 * 10^-6]);
xlabel('Radius (in m)');
ylabel('Working Range (in Pa)');
title('Working Range Vs Radius');
legend('TPP', 'Upper limit');
grid on;

format long

P_vector = linspace(0.3 * 10^5, 11 * 10^5, 101);
C = zeros(1, numel(P_vector));

% Parameters
R = 330 * 10^-6;
hmax = 14 * 10^-6;
g = 8 * 10^-6;
t = 0.1 * 10^-6;
E = 170 * 10^9;
v = 0.28;
eps_o = 8.854 * 10^-12;
eps_i = 3.7;

% Calculating rs
rs = hmax / 4 + R^2 / hmax;

% Calculating D
D = E * (hmax / 2)^3 / (12 * (1 - v^2));

% Calculating r_n1
syms s p
assume(s >= 0 & s < R);
equation = p * (R^2 - s^2)^2 / (64 * D) - g + t + sqrt(rs^2 - s^2) - sqrt(rs^2 - R^2);
for i = 1 : numel(P_vector)
    eq = subs(equation, p, P_vector(i));
    solutions = solve(eq, s, 'Real', true);
    r_n1 = double(solutions);
    disp(r_n1);

    % Calculating w_rn1
    w_rn1 = g - t - sqrt(rs^2 - r_n1^2) + sqrt(rs^2 - R^2);
    
    % Calculating m
    term = 2 * t + eps_i * g - eps_i * sqrt(4 * rs^2 - (r_n1 + R)^2) / 4 + eps_i * sqrt(rs^2 - R^2) - eps_i * t;
    m = -eps_i * w_rn1 / term;
    
    % Calculating k1 and k2
    sqrt_m = sqrt(m);
    k1 = atan(sqrt_m) / (2 * sqrt_m);
    t1 = sqrt(m) - m;
    t2 = sqrt(m) + m;
    k2 = (atan(sqrt(m / t1)) / sqrt(t1) + atanh(sqrt(m / t2)) / sqrt(t2)) / 2;
    
    % Calculating Capacitance, C
    r_n2 = R - r_n1;
    C(i) = pi * eps_o * eps_i * r_n1^2 / (2 * t) + 2 * pi * eps_o * eps_i * r_n2 * (r_n2 * k1 + r_n1 * k2) / term;
end

% Plotting the graph
plot(P_vector, C, 'r', 'LineWidth', 1.5);
xlim([0, 11 * 10^5]);
xlabel('Pressure (in Pa)');
ylabel('Capacitance (in Farad)');
title('Capacitance Vs Pressure');
grid on;

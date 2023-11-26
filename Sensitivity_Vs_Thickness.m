format long

hmax = linspace(7 * 10^-6, 18 * 10^-6, 101);
S = zeros(1, numel(hmax));

% Constants
R = 330 * 10^-6;
g = 8.5 * 10^-6;
t = 0.1 * 10^-6;
E = 170 * 10^9;
v = 0.28;
eps_i = 3.7;
eps_o = 8.854 * 10^-12;
P = 1 * 10^5;

for i = 1 : numel(hmax)
    % Taking g = hmax / 2 + (1 * 10^-6)
    g = hmax(i) / 2 + 1 * 10^-6;
    
    % Calculating rs
    rs = hmax(i) / 4 + R^2 / hmax(i);
    
    % Calculating D
    D = E * (hmax(i) / 2)^3 / (12 * (1 - v^2));
    
    % Calculating r_n1
    syms s
    assume(s >= 0 & s < R);
    equation = P * (R^2 - s^2)^2 / (64 * D) - g + t + sqrt(rs^2 - s^2) - sqrt(rs^2 - R^2);
    solution = solve(equation, s, 'Real', true);
    r_n1 = double(solution);
    disp(r_n1)

    % Calculating dr_n1/dp
    dr_n1BydP = (R^2 - r_n1^2)^2 / (r_n1 * (4 * P * (R^2 - r_n1^2) + 64 * D / sqrt(rs^2 - r_n1^2)));

    % Calculating w_rn1
    w_rn1 = g - t - sqrt(rs^2 - r_n1^2) + sqrt(rs^2 - R^2);
    
    term = 2 * t + eps_i * g - eps_i * sqrt(4 * rs^2 - (r_n1 + R)^2) / 4 + eps_i * sqrt(rs^2 - R^2) - eps_i * t;    % Frequently used term

    % Calculating m
    m = -eps_i * w_rn1 / term;

    % Calculating dm/dr_n1
    t1 = eps_i^2 * w_rn1 * (r_n1 + R) / (4 * sqrt(4 * rs^2 - (r_n1 + R)^2));
    t2 = eps_i * r_n1 * term / sqrt(rs^2 - r_n1^2);
    dmBydr_n1 = (t1 - t2) / term^2;
    
    sqrt_m = sqrt(m);

    % Calculating k1 and k2
    k1 = atan(sqrt_m) / (2 * sqrt_m);
    t1 = sqrt(m) - m;
    t2 = sqrt(m) + m;
    k2 = (atan(sqrt(m / t1)) / sqrt(t1) + atanh(sqrt(m / t2)) / sqrt(t2)) / 2;
    
    % Calculating dk1/dm and dk2/dm
    dk1Bydm = (sqrt_m - atan(sqrt_m) * (1 + m)) / (4 * sqrt_m * m * (1 + m));

    dk2Bydm = (sqrt(t1) - (1 - 2 * sqrt(m)) * atan(sqrt(m / t1)) / (4 * t1) + sqrt(t2) - (1 + 2 * sqrt(m)) * atanh(sqrt(m / t2)) / (4 * t2)) / (2 * sqrt(m));
    
    const = 2 * pi * eps_o * eps_i;
    
    % Calculating dC/dr_n1
    t1 = const * ((R - 2 * r_n1) * k2 - 2 * (R - r_n1) * k1) * term;
    t2 = const * (R - r_n1) * ((R - r_n1) * k1 - r_n1 * k2) * eps_i * (r_n1 + R) / (4 * sqrt(4 * rs^2 - (r_n1 + R)^2));
    dCBydr_n1 = const * r_n1 / (2 * t) + (t1 - t2) / term^2;

    r_n2 = R - r_n1;

    % Calculating dC/dk1
    dCBydk1 = const * r_n2^2 / term;

    % Calculating dC/dk2
    dCBydk2 = const * r_n1 * r_n2 / term;

    % Calculating S
    S(i) = (dCBydr_n1 + dCBydk1 * dk1Bydm * dmBydr_n1 + dCBydk2 * dk2Bydm * dmBydr_n1) * dr_n1BydP;
end

% Plotting the graph
plot(hmax, S, 'b', 'LineWidth', 1.5);
xlim([7 * 10^-6, 18 * 10^-6]);
xlabel('Thickness (in m)');
ylabel('Capacitive Sensitivity (in Fm^2 /N)');
title('Capacitive Sensitivity Vs Thickness');
grid on;

clc; % Clear command windowSSS

% Define global variables for the neuron model
global alpha delta mu phi A omega

%% Parameters for Period 1
alpha = 0.7;
delta = 0.8;
mu = 0.1;
phi = 0.15;
A = 1;
omega = 2;

%% Uncomment these lines for different periods
% Parameters for Period 2
% alpha = 1.1;
% delta = 0.8;
% mu = 0.3;
% phi = 0.15;
% A = 1.3;
% omega = 0.942;

%% Parameters for Period 3
% alpha = 0;
% delta = 0.8;
% mu = 0.1;
% phi = 0.2;
% A = 1;
% omega = 0.942;

% Initial conditions for the system
initial_conditions = [displacement(1,1); displacement(2,1)];
time_step = 0.01; % Time step for the ODE solver
time_span = 0:time_step:10*2*pi/omega; % Time span for the simulation

% Options for the ODE solver
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

% Solve the differential equations using ode45 solver
[time, solution] = ode45(@neuron, time_span, initial_conditions, options);

% Extract x1 and x2 from the solution matrix
x1 = solution(:, 1);
x2 = solution(:, 2);

% Calculate the derivatives of x1 and x2
dx1 = x1 * (1 - phi) - (1/3) * x1.^3 - x2 + A * cos(omega * time);
dx2 = mu * (x1 + alpha - delta * x2);

% Plot the response of x1 over time
figure;
plot(time(1:end), solution(1:end, 1), 'k.', 'MarkerSize', 10);
legend_x1 = legend('$${\rm RK}-x_1$$');
set(legend_x1, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Plot the response of x2 over time
figure;
plot(time(1:end), solution(1:end, 2), 'k.', 'MarkerSize', 10);
legend_x2 = legend('$${\rm RK}-x_2$$');
set(legend_x2, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Plot phase diagrams for x1 and x2
figure;
phase_plot_x1 = plot(x1(fix(7*2*pi/omega/time_step):16:end), dx1(fix(7*2*pi/omega/time_step):16:end), 'k.', 'MarkerSize', 10);
hold on;
phase_plot_x2 = plot(x2(fix(7*2*pi/omega/time_step):25:end), dx2(fix(7*2*pi/omega/time_step):25:end), 'k.', 'MarkerSize', 10);
legend_phase = legend([phase_plot_x1, phase_plot_x2], {'$$ {\rm RK}-x_1 $$','$${\rm RK}-x_2$$'});
set(legend_phase, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);


clear; clc;
close all;
tic;

% Global variables for the system parameters
global final_time step_size coefficient_r1 coefficient_r2 num_degrees_freedom num_harmonics time_data
global matrix_C matrix_K omega frequency_updated harmonic_coefficients alpha delta mu phi A

% Initialize parameters
num_degrees_freedom = 2;
num_harmonics = 30; % Adjusting num_harmonics affects solution accuracy, higher num_harmonics results in higher accuracy.

%% Figure 2-3 Related parameters
% alpha = 0.7;delta = 0.8;mu = 0.1;phi = 0.15;A = 1;
% omega = 2;
% frequency_updated = omega; % Period 1

%% Figure 4-7 Related parameters
% alpha=1.1;delta=0.8;mu=0.3;phi=0.15;A=1.3;
% omega=0.942;
% frequency_updated=omega/2; % Period 2

%% Figure 8-10 Related parameters
alpha=0;delta=0.8;mu=0.1;phi=0.2;A=1;
omega=0.942;
frequency_updated=omega/3; % Period 3

%% Initialize coefficients
coefficient_r1 = 1/3;
coefficient_r2 = 0;

% Define system matrices
matrix_C = [1, 0; 0, 1];
matrix_K = [phi-1, 1; -mu, delta*mu];

% Define time data
final_time = 2 * pi / frequency_updated;
step_size = 2 * pi / (frequency_updated * 1000);
time_data = (0:step_size:final_time);

% Initialize harmonic coefficients matrix
harmonic_coefficients = zeros(num_harmonics + 1, 2 * num_degrees_freedom);
harmonic_coefficients(1,:) = [-1, 0, -0.5, 0]; % Period 1
harmonic_coefficients(2,:) = [0, 0.5, 0, 0];

% parameter_set(1,:)=[-0.46,0,0.29,0];% Period 2
% parameter_set(2,:)=[1.76,1.16,-0.21,0.6];

% parameter_set(1,:)=[-0.46,0,0.29,0];% Period 3
% parameter_set(2,:)=[-4.5,1.16,-0.21,0.6];

initial_harmonic_coefficients = harmonic_coefficients;
iterations = length(time_data);

% Define a function to check if parameters are within a specified range
coefficients_check = @(parameter_set)(abs(parameter_set) < 10);

% Set trust region algorithm parameters
gamma_trust = 1.414;
rho_trust = 0.5;

% Record of parameter updates
harmonic_coefficients_record = harmonic_coefficients;
trust_region_record = [];

% Maximum iterations for response sensitivity
max_iterations_rs = 500;
max_iterations_tr = 20;

%% Iterate for response sensitivity
for iteration_index = 1:max_iterations_rs
    
    % Relative error tolerance for convergence
    error_tolerance = 1e-12;

    % Calculate residual for the current parameter matrix
    current_residual = calculate_residual(harmonic_coefficients);
    
    % Extract response sensitivity matrix
    sensitivity_parameter_updates = reshape(current_residual(:, num_degrees_freedom+1:2*num_degrees_freedom), num_degrees_freedom * length(time_data), 1);
    for i = 1:2*num_harmonics*num_degrees_freedom + 1
        sensitivity_parameter_updates = [sensitivity_parameter_updates, reshape(current_residual(:, num_degrees_freedom*(i+1)+1:num_degrees_freedom*(i+2)), num_degrees_freedom * length(time_data), 1)];
    end
     % Calculate neg_residual (negative of the first column of residual)
    neg_residual = -reshape(current_residual(:,1:num_degrees_freedom), num_degrees_freedom * length(time_data), 1);
    
    % Compute Tikhonov regularization using SVD
    [U, s, V] = csvd(sensitivity_parameter_updates);
    lambda_optimal = l_curve(U, s, neg_residual);
    temp_harmonic_coefficients = harmonic_coefficients;
    
    %% Trust-region algorithm
    for trust_index = 1:max_iterations_tr
       % Compute real_coeff_updates (optimal parameter increment)
        real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_optimal);

        % Reshape real_coeff_updates to match the harmonic coefficients matrix structure
        coeff_updates = reshape(real_coeff_updates, 2, num_degrees_freedom * num_harmonics + 1)';
        
        % Extract sensitivity_matrix (formerly sensitivity_parameter_da)
        temp_sensitivity_parameter_da = coeff_updates(1:num_harmonics, 1:2);
        for dof_index = 1:num_degrees_freedom - 1
            temp_sensitivity_parameter_da = [temp_sensitivity_parameter_da, coeff_updates(dof_index*num_harmonics + 1:(dof_index + 1)*num_harmonics, 1:2)];
        end
        
        sensitivity_matrix = zeros(1, 2 * num_degrees_freedom);
        sensitivity_matrix(1,1) = coeff_updates(num_degrees_freedom * num_harmonics + 1, 1);
        sensitivity_matrix(1,3) = coeff_updates(num_degrees_freedom * num_harmonics + 1, 2);
        sensitivity_matrix = [sensitivity_matrix; temp_sensitivity_parameter_da];
        
        % Check if the updated parameter is within the specified range
        if ~coefficients_check(temp_harmonic_coefficients + sensitivity_matrix)
            % If not, increase lambda until the parameter is within range
            lambda_optimal = lambda_optimal * gamma_trust;
            continue;
        end
        
         % Re-compute real_coeff_updates after adjusting lambda
        real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_optimal);

        % Store the initial coeff_updates
        initial_da = real_coeff_updates;

        % Reshape real_coeff_updates to match the harmonic coefficients matrix structure
        coeff_updates = reshape(real_coeff_updates, 2, num_degrees_freedom * num_harmonics + 1)';

        % Extract sensitivity_matrix (formerly sensitivity_parameter_da)
        temp_sensitivity_parameter_da = coeff_updates(1:num_harmonics, 1:2);
        for dof_index = 1:num_degrees_freedom - 1
            temp_sensitivity_parameter_da = [temp_sensitivity_parameter_da, coeff_updates(dof_index * num_harmonics + 1:(dof_index + 1) * num_harmonics, 1:2)];
        end
        
        sensitivity_matrix = zeros(1, 2 * num_degrees_freedom);
        sensitivity_matrix(1,1) = coeff_updates(num_degrees_freedom * num_harmonics + 1, 1);
        sensitivity_matrix(1,3) = coeff_updates(num_degrees_freedom * num_harmonics + 1, 2);
        sensitivity_matrix = [sensitivity_matrix; temp_sensitivity_parameter_da];
        
        % Update the harmonic coefficients matrix
        harmonic_coefficients = temp_harmonic_coefficients + sensitivity_matrix;
        
        % Re-calculate residual using the updated harmonic coefficients matrix
        residual_updates = calculate_residual(harmonic_coefficients);
        temp_residual_vector = -reshape(residual_updates(:,1:num_degrees_freedom), num_degrees_freedom * length(time_data), 1);
        
        % Calculate L_neg_residual (dot product of sensitivity matrix and initial_coeff_updates)
        L_neg_residual = sensitivity_parameter_updates * initial_da - neg_residual;

        % Calculate agreement indicator (rhos)
        rhos = (neg_residual' * neg_residual - temp_residual_vector' * temp_residual_vector) / (neg_residual' * neg_residual - L_neg_residual' * L_neg_residual);
        
        % Check if rhos meets the trust region criteria\
        if rhos >= rho_trust
            break;
        end
        lambda_optimal = lambda_optimal * gamma_trust;
    end
    
    % Calculate tolerance for convergence
    tolerance = norm(coeff_updates) / norm(harmonic_coefficients);

    % Store the harmonic coefficients matrix and lambda_inverse
    harmonic_coefficients_record = [harmonic_coefficients_record, harmonic_coefficients];
    trust_region_record = [trust_region_record; lambda_optimal];

    % Display the current harmonic coefficients matrix
    harmonic_coefficients;
    
    % Check for convergence
    if tolerance <= error_tolerance
        break;
    end
    
    % Store harmonic coefficients matrix at each iteration
    every_a(iteration_index).harmonic_coefficients = harmonic_coefficients;
    iteration_index
end
toc;
%%
% Calculate the final residuals
final_residual = calculate_residual(harmonic_coefficients);

% Plot the residuals
figure;
plot(time_data, final_residual(:,1), 'r-', 'LineWidth', 1.5);
hold on;
plot(time_data, final_residual(:,2), 'b-', 'LineWidth', 1.5);
legend('$$x_1$$', '$$x_2$$', 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

%% Calculate response for a longer time range
time_points_analysis = 0:0.01:10 * 2 * pi / omega;
harmonic_parameters = harmonic_coefficients(2:end, :);

% Initialize displacement and velocity responses
displacement = zeros(num_degrees_freedom, length(time_points_analysis));
velocity = zeros(num_degrees_freedom, length(time_points_analysis));

% Compute displacement and velocity responses
for j = 1:num_degrees_freedom
    for harmonic = 1:num_harmonics
        displacement(j, :) = displacement(j, :) + harmonic_parameters(harmonic, 2 * j - 1) * cos(harmonic * frequency_updated * time_points_analysis) + harmonic_parameters(harmonic, 2 * j) * sin(harmonic * frequency_updated * time_points_analysis);
        velocity(j, :) = velocity(j, :) - frequency_updated * harmonic * harmonic_parameters(harmonic, 2 * j - 1) * sin(harmonic * frequency_updated * time_points_analysis) + frequency_updated * harmonic * harmonic_parameters(harmonic, 2 * j) * cos(harmonic * frequency_updated * time_points_analysis);
    end
    displacement(j, :) = displacement(j, :) + harmonic_coefficients(1, 2 * j - 1);
end

% Plot displacement responses
figure;
plot(time_points_analysis, displacement(1, :), 'r', 'LineWidth', 1.5);
legend('$$TMRM-x_1$$', 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

figure;
plot(time_points_analysis, displacement(2, :), 'b', 'LineWidth', 1.5);
legend('$$TMRM-x_2$$', 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Plot phase diagrams
figure;
plot(displacement(1, fix(7 * 2 * pi / omega / 0.01):end), velocity(1, fix(7 * 2 * pi / omega / 0.01):end), 'r', 'LineWidth', 1.5);
hold on;
plot(displacement(2, fix(7 * 2 * pi / omega / 0.01):end), velocity(2, fix(7 * 2 * pi / omega / 0.01):end), 'b', 'LineWidth', 1.5);
legend('$$TMRM-x_1$$', '$$TMRM-x_2$$', 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

figure;
plot(displacement(1, fix(7 * 2 * pi / omega / 0.01):end), displacement(2, fix(7 * 2 * pi / omega / 0.01):end), 'b', 'LineWidth', 1.5);
legend('$$TMRM-x_2$$', 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

initial_harmonic_coefficients

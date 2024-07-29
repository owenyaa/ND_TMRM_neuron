function residual = calculate_residual(harmonic_coefficients)
    % Global variables for the system parameters
    global coefficient_r1 coefficient_r2 num_degrees_freedom num_harmonics time_data
    global matrix_C matrix_K A omega frequency_updated alpha delta mu
    
    % Extract harmonic parameters
    harmonic_parameters = harmonic_coefficients(2:end,:);
    
    % Initialize displacement and velocity responses
    displacement = zeros(num_degrees_freedom, length(time_data));
    velocity = zeros(num_degrees_freedom, length(time_data));
    
    % Calculate displacement and velocity for each degree of freedom
    for dof = 1:num_degrees_freedom
        for harmonic = 1:num_harmonics
            displacement(dof, :) = displacement(dof, :) + harmonic_parameters(harmonic, 2*dof-1) * cos(harmonic * frequency_updated * time_data) + harmonic_parameters(harmonic, 2*dof) * sin(harmonic * frequency_updated * time_data);
            velocity(dof, :) = velocity(dof, :) - frequency_updated * harmonic * harmonic_parameters(harmonic, 2*dof-1) * sin(harmonic * frequency_updated * time_data) + frequency_updated * harmonic * harmonic_parameters(harmonic, 2*dof) * cos(harmonic * frequency_updated * time_data);
        end
        displacement(dof, :) = displacement(dof, :) + harmonic_coefficients(1, 2*dof-1);
    end
    
    % Nonlinear matrix initialization
    damping_matrix_nonlinear = zeros(num_degrees_freedom, num_degrees_freedom);
    damping_matrix_nonlinear(1,1) = coefficient_r1;
    damping_matrix_nonlinear(2,2) = coefficient_r2;
    
    % External forcing function
    force_matrix = zeros(num_degrees_freedom, length(time_data));
    force_matrix(1, :) = A * cos(omega * time_data);
    force_matrix(2, :) = alpha * mu;
    
    % Calculate residuals for displacement and velocity
    residual(1:num_degrees_freedom, :) = matrix_C * velocity + matrix_K * displacement + damping_matrix_nonlinear * displacement.^3 - force_matrix;
    
    %% Compute sensitivities for harmonic coefficients
    for i = 1:2 * num_harmonics * num_degrees_freedom
        sensitivity_parameter = zeros(2 * num_harmonics * num_degrees_freedom, 1);
        displacement_sensitivity = zeros(num_degrees_freedom, length(time_data));
        velocity_sensitivity = zeros(num_degrees_freedom, length(time_data));
        
        % Set sensitivity parameter
        sensitivity_parameter(i, 1) = 1;
        sensitivity_parameter = reshape(sensitivity_parameter, 2, num_harmonics * num_degrees_freedom)';
        sensitivity_parameters = sensitivity_parameter(1:num_harmonics, 1:2);
        
        for dof_index = 1:num_degrees_freedom-1
            sensitivity_parameters = [sensitivity_parameters, sensitivity_parameter(dof_index*num_harmonics+1:(dof_index+1)*num_harmonics, 1:2)];
        end
        
        % Compute displacement and velocity sensitivities
        for dof = 1:num_degrees_freedom
            for harmonic = 1:num_harmonics
                displacement_sensitivity(dof, :) = displacement_sensitivity(dof, :) + sensitivity_parameters(harmonic, 2*dof-1) * cos(harmonic * frequency_updated * time_data) + sensitivity_parameters(harmonic, 2*dof) * sin(harmonic * frequency_updated * time_data);
                velocity_sensitivity(dof, :) = velocity_sensitivity(dof, :) - frequency_updated * harmonic * sensitivity_parameters(harmonic, 2*dof-1) * sin(harmonic * frequency_updated * time_data) + frequency_updated * harmonic * sensitivity_parameters(harmonic, 2*dof) * cos(harmonic * frequency_updated * time_data);
            end
        end
        
        % Calculate residuals for sensitivities
        residual(num_degrees_freedom*i+1:num_degrees_freedom*(i+1), :) = matrix_C * velocity_sensitivity + matrix_K * displacement_sensitivity + 3 * damping_matrix_nonlinear * (displacement.^2 .* displacement_sensitivity);
    end
    
    %% Compute sensitivities for constant terms
    for dof = 1:num_degrees_freedom
        displacement_sensitivity_constant = zeros(num_degrees_freedom, length(time_data));
        velocity_sensitivity_constant = zeros(num_degrees_freedom, length(time_data));
        displacement_sensitivity_constant(dof, :) = ones(1, length(time_data));
        
        % Calculate residuals for constant sensitivities
        residual(num_degrees_freedom*(2*num_harmonics*num_degrees_freedom + dof)+1:num_degrees_freedom*(2*num_harmonics*num_degrees_freedom + 1 + dof), :) = matrix_C * velocity_sensitivity_constant + matrix_K * displacement_sensitivity_constant + 3 * damping_matrix_nonlinear * (displacement.^2 .* displacement_sensitivity_constant);
    end
    
    residual = residual';
end

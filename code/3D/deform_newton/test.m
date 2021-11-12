function [energy_vector] = test(kernel, make_plots)
addpath(genpath('../'))

% example_elephant_init
% example_large_init
% example_shear_bar_init
% example_twist_bar_init
% example_stretch_bar_init
% example_rand_mesh_init
% example_zero_mesh_init
% example_armadillo_init
% example_dancer_init
% example_botijo_init
% example_dilo_test
% example_homer_init
% example_horse_init
% example_cube_init
% example_wrench_init
% example_statue_init
% example_armadillo_dance_init
% example_homer_bend_init

% Each element of param_groups is a set of hyperparameters defining an
% experiment
param_groups = [struct(), struct()];

% Parameters for default run
param_groups(1).preconditioner = ""; 
param_groups(1).pcg_parameters = struct();
param_groups(1).use_direct = true;

% Parameters for adaptive pcg run
param_groups(2).preconditioner = "incomplete_LU";
% preconditioner = 'diagonal'; 
param_groups(2).pcg_parameters = struct( ...
    'tol', 1e-4, ...
    'restart', 100, ...
    'maxit', 100, ...
    'use_warm_start', false, ...
    'energy_tol', 1e-4, ...
    'line_check_jump', 5);
param_groups(2).use_direct = false; 

energies = cell(length(param_groups), 1);
for i = 1:length(param_groups)
    % Set up the experiment
    if kernel == "elephant"
        example_elephant_init 
    elseif kernel == "shear_bar"
        example_shear_bar_init
    elseif kernel == "botijo"
        example_botijo_init
    elseif kernel == "armadillo"
        example_armadillo_init
    end
    energies{i} = newton_solver(u_n, param_groups(i).preconditioner, param_groups(i).pcg_parameters, param_groups(i).use_direct, kernel, i, make_plots);
end

if make_plots
    energy_figure = figure;
    figure(energy_figure);
    hold on;
    title(['Energy Plot (', kernel, ')']);
    xlabel('Newton Iteration');
    ylabel('Energy Value');
    for j = 1:length(param_groups)
        plot(energies{j}, 'DisplayName', ['Param Group - ', num2str(j)]);
    end
    legend();
    hold off;
    saveas(energy_figure, strcat(kernel, '_energy.fig'))
end 

end


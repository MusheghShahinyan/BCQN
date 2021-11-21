function [] = test(kernel, make_plots, re_run)
%TEST run tests for the given kernel
%
% kernel: Kernel name (e.g. "shear_bar" or "elephant")
% make_plots: A boolean. Whether to produce plots or not
% re_run: The indices of experiments to re_run, as opposed to loading
% the results from a cache file 

addpath(genpath('../'))

if nargin <= 2
  re_run = 0;
end

% Each element of param_groups is a set of hyperparameters defining an
% experiment
param_groups = [struct(), struct(), struct()];

%% Parameters for default run
param_groups(1).name = "Direct Solver";
param_groups(1).short_name = "dirsolve";
param_groups(1).preconditioner = ""; 
param_groups(1).pcg_parameters = struct();
param_groups(1).use_direct = true;

%% Parameters for default fixed threshold pcg run
param_groups(2).preconditioner = "incomplete_LU";
param_groups(2).name = "Iterative Golden (ilu)";
param_groups(2).short_name = "pcg_ilu";
param_groups(2).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 200, ...
    'use_warm_start', false ...
);
param_groups(2).use_direct = false; 

%% Parameters for adaptive pcg run
param_groups(3).preconditioner = "incomplete_LU";
param_groups(3).name = "Adaptive Iterative (ilu)";
param_groups(3).short_name = "apcg_ilu";
param_groups(3).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 200, ...
    'use_warm_start', false, ...
    'energy_tol', 1e-4, ...
    'line_check_jump', 5 ...
);
param_groups(3).use_direct = false; 

function [results] = run_newton_solver(i, u_n)
    %RUN_NEWTON_SOLVER runs newton_solver for the parameter group i

    disp(['Running newton_solver for param_group(', num2str(i), ')']);
    results = newton_solver( ...
        u_n, ...
        param_groups(i).preconditioner, ...
        param_groups(i).pcg_parameters, ...
        param_groups(i).use_direct, ...
        kernel, ...
        i, ...
        make_plots ...
    );
    save_to_cache(kernel, param_groups(i), results)
end

%% Run the experiments 

param_group_results = cell(length(param_groups), 1);

for i = 1:length(param_groups)
    [~, u_n] = initialize_kernel(kernel);  
    
    [success, results] = load_results_from_cache(kernel, re_run, param_groups(i), i);

    if not(success)
        results = run_newton_solver(i, u_n);
    end

    param_group_results{i} = results;
    
end

if make_plots
    close all

    plot_pcg(param_group_results, param_groups, kernel);
    plot_energy(param_group_results, param_groups, kernel);
    plot_energy_grad(param_group_results, param_groups, kernel);
    plot_est_grad(param_group_results, param_groups, kernel);
    % plot_tsne(param_group_results, param_groups, kernel, u_n);

end 

end


function [] = test(kernel, make_plots, re_run, render)
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
param_groups = [];

%% Parameters for default run
param_groups(1).name = "Direct Solver";
param_groups(1).short_name = "dirsolve";
param_groups(1).preconditioner = ""; 
param_groups(1).pcg_parameters = struct();
param_groups(1).use_direct = true;
param_groups(1).stop_criterion = "";
% param_groups(1).stop_criterion = "zeroth_third";

%% Parameters for default fixed threshold pcg run
param_groups(2).preconditioner = "incomplete_LU";
param_groups(2).name = "Iterative Golden (ilu, cg 10)";
param_groups(2).short_name = "pcg_ilu_cg_10";
param_groups(2).use_direct = false; 
param_groups(2).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 10, ...
    'use_warm_start', true ...
);
param_groups(2).stop_criterion = "fixed_steps50";

%% Parameters for default fixed threshold pcg run
param_groups(3).preconditioner = "incomplete_LU";
param_groups(3).name = "Iterative Golden (ilu, cg 50)";
param_groups(3).short_name = "pcg_ilu_cg_50";
param_groups(3).use_direct = false; 
param_groups(3).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 50, ...
    'use_warm_start', true ...
);
param_groups(3).stop_criterion = "";

%% Parameters for default fixed threshold pcg run
param_groups(4).preconditioner = "incomplete_LU";
param_groups(4).name = "Iterative Golden (ilu, cg 100)";
param_groups(4).short_name = "pcg_ilu_cg_100";
param_groups(4).use_direct = false; 
param_groups(4).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 100, ...
    'use_warm_start', true ...
);
param_groups(4).stop_criterion = "";


%% Parameters for default fixed threshold pcg run
param_groups(5).preconditioner = "incomplete_LU";
param_groups(5).name = "Iterative Golden (ilu, cg 250)";
param_groups(5).short_name = "pcg_ilu_cg_250";
param_groups(5).use_direct = false; 
param_groups(5).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 250, ...
    'use_warm_start', true ...
);
param_groups(5).stop_criterion = "";


%% Parameters for adaptive pcg run
param_groups(6).preconditioner = "incomplete_LU";
param_groups(6).name = "Adaptive Iterative (grad check 10)";
param_groups(6).short_name = "apcg_ilu_grad_10";
param_groups(6).use_direct = false; 
param_groups(6).use_custom_pcg = true;
param_groups(6).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 10, ...
    'use_warm_start', true, ...
    'check_grad', true, ...
    'check_grad_jump', 7 ...
    ... %'stopping_pairs', [0 50] ...
);
param_groups(6).stop_criterion = "";

%% Parameters for adaptive pcg run
param_groups(7).preconditioner = "incomplete_LU";
param_groups(7).name = "Adaptive Iterative (grad check 50)";
param_groups(7).short_name = "apcg_ilu_grad_50";
param_groups(7).use_direct = false; 
param_groups(7).use_custom_pcg = true;
param_groups(7).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 50, ...
    'use_warm_start', true, ...
    'check_grad', true, ...
    'check_grad_jump', 7 ...
    ... %'stopping_pairs', [0 50] ...
);
param_groups(7).stop_criterion = "";

%% Parameters for adaptive pcg run
param_groups(8).preconditioner = "incomplete_LU";
param_groups(8).name = "Adaptive Iterative (grad check 100)";
param_groups(8).short_name = "apcg_ilu_grad_100";
param_groups(8).use_direct = false; 
param_groups(8).use_custom_pcg = true;
param_groups(8).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 100, ...
    'use_warm_start', true, ...
    'check_grad', true, ...
    'check_grad_jump', 7 ...
    ... %'stopping_pairs', [0 50] ...
);
param_groups(8).stop_criterion = "";

%% Parameters for adaptive pcg run
param_groups(9).preconditioner = "incomplete_LU";
param_groups(9).name = "Adaptive Iterative (grad check 250)";
param_groups(9).short_name = "apcg_ilu_grad_250";
param_groups(9).use_direct = false; 
param_groups(9).use_custom_pcg = true;
param_groups(9).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 250, ...
    'use_warm_start', true, ...
    'check_grad', true, ...
    'check_grad_jump', 7 ...
);
%}

function [results] = run_newton_solver(i, u_n)
    %RUN_NEWTON_SOLVER runs newton_solver for the parameter group i

    disp(['Running newton_solver for param_group(', num2str(i), ')']);
    results = newton_solver( ...
        u_n, ...
        param_groups(i).preconditioner, ...
        param_groups(i).pcg_parameters, ...
        param_groups(i).use_direct, ...
        param_groups(i).use_custom_pcg, ...
        param_groups(i).stop_criterion, ...
        kernel, ...
        i, ...
        make_plots ...
    );
    save_to_cache(kernel, param_groups(i), results)
end

%% Run the experiments 

param_group_results = cell(length(param_groups), 1);  
kernel_initialized = false;

for i = 1:length(param_groups)
    
    [success, results] = load_results_from_cache(kernel, re_run, param_groups(i), i);

    if not(success)
        [~, u_n] = initialize_kernel(kernel);
        kernel_initialized = true;
        results = run_newton_solver(i, u_n);
    end

    %{
    if render(1) == i
        if not(kernel_initialized)
           [~, ~] = initialize_kernel(kernel); 
        end
        size(results.uns)
        opengl_render_result( results.uns(:, render(2)), 0 )
    end
    %}
    
    param_group_results{i} = results;
end

if make_plots
    close all

     %plot_residual_comparison(param_group_results, param_groups, kernel, [5 9]);
     cumulative_cg(param_group_results, param_groups, kernel);
     % plot_energy(param_group_results, param_groups, kernel);
     % plot_pcg(param_group_results, param_groups, kernel);
     % plot_energy_grad(param_group_results, param_groups, kernel);
     % plot_est_grad(param_group_results, param_groups, kernel);
     % plot_pcg_grad(param_group_results, param_groups, kernel);

   	 % plot_tsne(param_group_results, param_groups, kernel, u_n);

end 

end


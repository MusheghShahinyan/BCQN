function [] = test(kernel, make_plots, re_run)

addpath(genpath('../'))

if nargin <= 2
  re_run = 0;
end

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
param_groups = [struct(), struct(), struct(), struct()];

%% Parameters for default run
param_groups(1).name = "Direct Solver";
param_groups(1).preconditioner = ""; 
param_groups(1).pcg_parameters = struct();
param_groups(1).use_direct = true;

%% Parameters for default fixed threshold pcg run
param_groups(2).preconditioner = "incomplete_LU";
param_groups(2).name = "Iterative Golden (ilu)";
% preconditioner = 'diagonal'; 
param_groups(2).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 400, ...
    'use_warm_start', false);
param_groups(2).use_direct = false; 

%% Parameters for adaptive pcg run
param_groups(3).preconditioner = "incomplete_LU";
% preconditioner = 'diagonal'; 
param_groups(3).name = "Adaptive Iterative (ilu)";
param_groups(3).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 400, ...
    'use_warm_start', false, ...
    'energy_tol', 1e-4, ...
    'line_check_jump', 5);
param_groups(3).use_direct = false; 

%% Parameters for adaptive pcg run
param_groups(4).preconditioner = "incomplete_LU";
% preconditioner = 'diagonal'; 
param_groups(4).name = "Adaptive Iterative (neg. allowed, ilu)";
param_groups(4).pcg_parameters = struct( ...
    'tol', 1e-8, ...
    'restart', 100, ...
    'maxit', 400, ...
    'use_warm_start', false, ...
    'energy_tol', 1e-4, ...
    'line_check_jump', 5, ...
    'allow_negative_energy_delta', true);
param_groups(4).use_direct = false;  

%%

param_group_results = cell(length(param_groups), 1);

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
    elseif kernel == "wrench"
        example_wrench_init
    elseif kernel == "homer" 
        example_homer_bend_init
    elseif kernel == "horse"
        example_horse_init
    elseif kernel == "cube"  
        example_cube_init
    elseif kernel == "dancer"    
        example_dancer_init
    elseif kernel == "twist"  
        example_twist_bar_init
    end
    
    cache_file = sprintf('run_cache_%s_%s.mat', kernel, DataHash(param_groups(i)));
    
    if isfile(cache_file) && ...
            not((isvector(re_run) && any(re_run == i)) ...
            || (isstring(re_run) && re_run == "all")) %...
            %|| (isnumeric(re_run) && re_run == i))
        disp(['Using cache for param_group(', num2str(i), ')']);
        
        load(cache_file, 'results');
        param_group_results{i} = results;
    else
        disp(['Running param_group(', num2str(i), ')']);
        results = newton_solver(u_n, param_groups(i).preconditioner, param_groups(i).pcg_parameters, param_groups(i).use_direct, kernel, i, make_plots);
        params = param_groups(i);

        save(cache_file, 'results', 'params');
        param_group_results{i} = results;
    end
end

if make_plots
    close all
    
    figure; energy_axes = axes; hold(energy_axes, "on");
    figure; energy_axes_grad = axes; hold(energy_axes_grad, "on");
    figure; est_grad = axes; hold(est_grad, "on");
    figure; eta_axes = axes; hold(eta_axes, "on"); 
    figure; rnorm_axes = axes; hold(eta_axes, "on"); 
    
    title(est_grad, ['Est. Grad', kernel]);

    title(energy_axes_grad, ['Energy Plot Grad', kernel]);
    xlabel(energy_axes_grad, 'Newton Iteration');
    ylabel(energy_axes_grad, 'Energy Value');
    
    title(energy_axes, ['Energy Plot ', kernel,]);
    xlabel(energy_axes, 'Newton Iteration');
    ylabel(energy_axes, 'Energy Value');
    
    title(eta_axes, ['Eta Plot (', kernel, ')']);
    xlabel(eta_axes, 'Newton Iteration');
    ylabel(eta_axes, 'Eta Value');
    
    est_grads = zeros(length(param_groups), 100);
    act_grads = zeros(length(param_groups), 100);
    labels = cell(1, length(param_groups)*2);
    max_timestep = 0;
    
    for j = 1:length(param_groups)
        results = param_group_results{j};
        pcg_parameters = param_groups(j).pcg_parameters;
        
        size(results.grads_pre_ls)
        estimated_grad = zeros(size(results.grads_pre_ls));
        
        for i = 1:length(results.grads_pre_ls(1,:))
            H = results.hessians_pre_ls{i};
            grad = results.grads_pre_ls(:, i);
            p = results.search_directions(:, i);
            estimated_grad(:, i) = grad + H*p;
        end
        
        labels{j} = param_groups(j).name + " est. grad";
        labels{length(param_groups) + j} = ...
            param_groups(j).name + " act. grad";
        
        est_grads(j, 1:length(estimated_grad(1,:))) = ...
            vecnorm(estimated_grad);
        act_grads(j, 1:length(results.grads_pre_ls(1,:))) = ...
            vecnorm(results.grads_pre_ls);
        
        % Shift x axis to align with newton steps
        disp(param_groups(j).name);
        
        % Shift x axis to align with newton steps + starting point
        x_axis_energies = linspace(-1, length(results.energies) - 2, length(results.energies));
        % Shift x axis to align with newton steps
        x_axis = linspace(0, length(results.etas) - 1, length(results.etas));
        
        plot(energy_axes, x_axis_energies, results.energies, 'DisplayName', param_groups(j).name);

        if max(max_timestep, length(results.bs(1,:)) > max_timestep)
            max_timestep = max(max_timestep, length(results.bs(1,:)));
        end
        
        plot(energy_axes_grad, x_axis_energies, results.energies, 'DisplayName', param_groups(j).name);
        yyaxis(energy_axes_grad, 'right')
        plot(energy_axes_grad, x_axis, vecnorm(results.bs), 'DisplayName', param_groups(j).name);
        set(energy_axes_grad,'Yscale','log');
        yyaxis(energy_axes_grad, 'left')
        
        if isfield(results, 'num_iter')
            yyaxis(energy_axes, 'right')
            ylim([0 250])
            plot(energy_axes, x_axis, results.num_iter, 'DisplayName', param_groups(j).name);
            yyaxis(energy_axes, 'left')
        end
        
        set(findall(energy_axes,'YAxisLocation','left'),'Yscale','log');
        plot(eta_axes, results.etas, 'DisplayName', param_groups(j).name);

        adaptive_pcg = not(param_groups(j).use_direct) ...
            && isfield(pcg_parameters, 'energy_tol') ...
            && isfield(pcg_parameters, 'line_check_jump');
        
        if adaptive_pcg
            disp("plotting pcg")
            figure; pcg_1 = subplot(3,1,1); pcg_2 = subplot(3,1,2); pcg_3 = subplot(3,1,3); 
            hold(pcg_1, "on"); hold(pcg_2, "on"); hold(pcg_3, "on"); 
            
            for iter = 1:length(results.resvecs)
                energyvec = results.energyvecs{iter};
                plot(pcg_1, (energyvec - min(energyvec)) / (energyvec(1) - min(energyvec)), 'DisplayName', ['iter ', num2str(iter - 1)]);
                plot(pcg_2, results.resvecs{iter}, 'DisplayName', ['iter ', num2str(iter - 1), 'res']);
                plot(pcg_3, results.anglesvecs{iter}, 'DisplayName', ['iter ', num2str(iter - 1), 'res']);
            end
            
            set(pcg_1,'Yscale','linear');
            legend(pcg_1);
            hold(pcg_1, "off");
            
            set(pcg_2,'Yscale','log');
            legend(pcg_2);
            hold(pcg_2, "off");
            
            set(pcg_3,'Yscale','linear');
            legend(pcg_3);
            hold(pcg_3, "off");
        end
        
    end
    
    % plot grad differences
    est_grads = est_grads(:, 1:max_timestep)
    act_grads = act_grads(:, 1:max_timestep)
    
    % shift est. grads since they estimate the gradient of the next
    %   timestep
    est_grads = [zeros(length(est_grads(:,1)), 1) est_grads]
    act_grads = [act_grads zeros(length(est_grads(:,1)), 1)]
    
    len = length(act_grads(1,:));
    x_axis = linspace(0, len - 1, len);
    bar(est_grad, x_axis, [est_grads; act_grads])
    set(est_grad,'Yscale','log');
    legend(est_grad, labels);
    
    set(energy_axes_grad,'Yscale','log');
    
    legend(energy_axes);
    legend(energy_axes_grad);
    legend(eta_axes);
    legend(est_grad);

    hold(energy_axes, "off");
    hold(energy_axes_grad, "off");
    hold(eta_axes, "off");
    hold(est_grad, "off");

    saveas(energy_axes, strcat(kernel, '_energy.fig'));
    saveas(eta_axes, strcat(kernel, '_eta.fig'));
end 

end


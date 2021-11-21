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
        
        figure; tsne_axes = axes; hold(tsne_axes, "on");
    
        figure; energy_axes = axes; hold(energy_axes, "on");
        figure; energy_axes_grad = axes; hold(energy_axes_grad, "on");
        figure; est_grad = axes; hold(est_grad, "on");
        figure; eta_axes = axes; hold(eta_axes, "on");
        
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
    
        total_guesses = 0;
        for j = 1:length(param_groups)
            results = param_group_results{j};
            total_guesses = total_guesses + size(results.guesses, 1);
        end
        tsne_X = zeros(total_guesses, length(u_n));
    
        cur_tsne_count = 1; 
    
        for j = 1:length(param_groups)
            results = param_group_results{j};
            pcg_parameters = param_groups(j).pcg_parameters;
    
            tsne_X(cur_tsne_count:(cur_tsne_count+size(results.guesses, 1)-1), :) = results.guesses;
            
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
                plot(energy_axes, x_axis, results.num_iter, 'DisplayName', strcat(param_groups(j).name, '- iters'));
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
                title(pcg_1, "Energy");
                hold(pcg_1, "off");
                
                set(pcg_2,'Yscale','log');
                legend(pcg_2);
                title(pcg_2, "Residuals")
                hold(pcg_2, "off");
                
                set(pcg_3,'Yscale','linear');
                legend(pcg_3);
                title(pcg_3, "Angles (w.r.t -grad)")
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
       
        % tsne_X tsne_axes 
        % tsne_Y = tsne(tsne_X, 'LearnRate',20000);
        % gscatter(tsne_axes, tsne_Y(:, 1), tsne_Y(:, 2));
        %legend(tsne_axes);
        title(tsne_axes, "Search Path Embeddign");
        hold(tsne_axes, "off");
    
        legend(energy_axes);
        legend(energy_axes_grad);
        legend(eta_axes);
        legend(est_grad);
    
        hold(energy_axes, "off");
        hold(energy_axes_grad, "off");
        hold(eta_axes, "off");
        hold(est_grad, "off");
    
        saveas(energy_axes, strcat('figures/', kernel, '_energy.fig'));
        saveas(eta_axes, strcat('figures/', kernel, '_eta.fig'));
    end 

end


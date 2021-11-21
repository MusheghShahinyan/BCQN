function [ results ] = newton_solver( un, preconditioner, pcg_parameters, use_direct, kernel, param_group_id, make_plots)
%NEWTON_SOLVER Summary of this function goes here
%   Detailed explanation goes here

u = un;
p = 0 * u;

L = get_precomputed_split_info();

order = get_ordering(u);

'newton'

x0 = zeros(length(u), 1);
energy_vector = zeros(2000, 1);
eta_vector = zeros(2000, 1);
final_angle_from_grad = zeros(2000, length(u));

grads_pre_ls = zeros(length(u), 2000);
b_vector = zeros(length(u), 2000);
search_directions = zeros(length(u), 2000);
hessians = cell(2000, 1);
hessians_pre_ls = cell(2000, 1);

energy_vector(1) = energy_value(un);
results = struct();

results.guesses = zeros(20, length(u));
results.guesses(1, :) = x0; 

adaptive_pcg = not(use_direct) && isfield(pcg_parameters, 'energy_tol') && isfield(pcg_parameters, 'line_check_jump');

if not(use_direct)
   results.num_iter = zeros(2000, 1);
   results.resvecs = cell(2000, 1);
   results.normalized_resvecs = cell(2000, 1);
   results.relres = zeros(2000, 1);

   if adaptive_pcg
      results.energyvecs = cell(2000, 1);
      results.anglesvecs = cell(2000, 1);
   end
end


% if make_plots
%     figure; 
%     pcg_axes_1 = subplot(2,1,1); hold(pcg_axes_1, "on"); 
%     pcg_axes_2 = subplot(2,1,2); hold(pcg_axes_2, "on"); 
% 
%     set(pcg_axes_1, 'YScale', 'log');
%     title(pcg_axes_1, strcat('PCG Curves. kernel: ', kernel, '. Param Group:', num2str(param_group_id), '.'))
%     ylabel(pcg_axes_1, 'normed res (log)');
%     xlabel(pcg_axes_1, 'iteration');
% end

for i = 0 : 2000
    
    [grad, H] = grad_hessian_function(u, 0);

    if preconditioner == "incomplete_LU"
        [L,U] = ilu(H, struct('type','nofill'));
    elseif preconditioner == "diagonal"
        L = spdiags(1./diag(H), 0, length(H),length(H));
        U = eye(size(L, 1));
    end
      
    if use_direct 
        % Original Implementation
        p(order) = -1.0 * H(order, order) \ grad(order);
        iterations = -1;
        flag = -1; 
        relres = norm((H * p) + grad) / norm(grad);
    else   
        % Preconditioned CG 

        % Warm start
        if pcg_parameters.use_warm_start
            x0 = p;
        else
            x0 = 0 * p;
        end 

        results.energies = energy_vector;
        results.etas = results;
        
        if adaptive_pcg
            %pcg_parameters.ignore_stop = (i == 1);
            [p, flag, relres, iterations, resvec, energyvec, anglesvec] ...
                = pcg_copy(H, -1.0 * grad, pcg_parameters.tol, ...
                pcg_parameters.maxit, L,U, x0, u, pcg_parameters);
            results.energyvecs{i+1} = energyvec;
            results.anglesvecs{i+1} = anglesvec;
        else
            [p, flag, relres, iterations, resvec] = pcg(H, -1.0 * grad, pcg_parameters.tol, pcg_parameters.maxit, L, U, x0);
        end
            
        results.num_iter(i+1) = length(resvec);
        results.relres(i+1) = relres;
        results.resvecs{i+1} = resvec;
        results.normalized_resvecs{i+1} = resvec / norm(-1.0 * grad);
    end 
    
    
%     % Plot PCG residual progress
%     if not(use_direct) && make_plots
%         %plot(pcg_axes, resvec / norm(-1.0 * grad), 'DisplayName', ['iter ', num2str(i)]);
%         disp("plotting");
%         plot(pcg_axes_1, energy / energy(1), 'DisplayName', ['iter ', num2str(i)]);
%         plot(pcg_axes_2, resvec, 'DisplayName', ['iter ', num2str(i), 'res']);
% 
%     end


    [grads_pre_ls(:, i+1), hessians_pre_ls{i+1}] = grad_hessian_function(u + p, 0);
    search_directions(:, i+1) = p;
    hessians{i + 1} = H;
    
    [un, alp, balp] = line_check_search(p, u, grad);
    energy = energy_value(un);
    
    b_vector(:, i+1) = -1.0 * grad;
    energy_vector(i+2) = energy;
    results.guesses(i+2, :) = un;
    eta_vector(i+1) = relres; 

    disp(['= Newton ', num2str(i), ' => ', ...
        ' alp: ', num2str(alp), ...
        ' balp: ', num2str(balp), ...
        ' num iter: ', num2str(iterations), ' relres: ', num2str(relres), ' energy(un): ', num2str(energy), ...
        ' pcg flag: ', num2str(flag)])
    
    if stop_check(un, u, grad)  
        break;
    end  
    
    u = un;
end

%energy_vector(i + 2) = energy_value(un);
energy_vector = energy_vector(1:(i+2), :);
eta_vector = eta_vector(1:(i+1), :);
final_angle_from_grad = final_angle_from_grad(1:(i+1), :);

results.energies = energy_vector;
results.etas = eta_vector;
results.bs = b_vector(1:(i+1), :); % lhs in Ax = b
results.guesses = results.guesses(1:(i+2), :);
results.final_angle_from_grad = final_angle_from_grad;
results.bs = b_vector(:, 1:(i+1)); % lhs in Ax = b
results.hessians_pre_ls = hessians_pre_ls(1:(i+1));
results.grads_pre_ls = grads_pre_ls(:, 1:(i+1));
results.search_directions = search_directions(:, 1:i+1);

if not(use_direct)
   results.num_iter = results.num_iter(1:(i+1), :);
   results.resvecs = results.resvecs(1:(i+1), :);
   results.normalized_resvecs = results.normalized_resvecs(1:(i+1), :);
   results.relres = results.relres(1:(i+1), :);
   
   if adaptive_pcg
      results.energyvecs = results.energyvecs(1:(i+1), :);
      results.anglesvecs = results.anglesvecs(1:(i+1), :);
   end
end

% if make_plots
%     set(pcg_axes_1,'Yscale','linear');
%     legend(pcg_axes_1);
%     hold(pcg_axes_1, "off");
%     saveas(pcg_axes_1, strcat(kernel, '_', num2str(param_group_id), '_pcg.fig'))
% end

end


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
b_vector = zeros(2000, length(u));
final_angle_from_grad = zeros(2000, length(u));

prev_energy = 0;

results = struct();

adaptive_pcg = not(use_direct) && isfield(pcg_parameters, 'energy_tol') && isfield(pcg_parameters, 'line_check_jump');
if not(use_direct)
   results.num_iter = zeros(2000, 1);
   results.resvecs = cell(2000, 1);
   results.normalized_resvecs = cell(2000, 1);
   results.rnorm = zeros(2000, 1);
   results.relative_rnorm = zeros(2000, 1);

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
        rnorm = norm((H * p) + grad);
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
            [p, flag, rnorm, iterations, resvec, energyvec, anglesvec] = pcg_copy(H, -1.0 * grad, pcg_parameters.tol, pcg_parameters.maxit, L,U, x0, u, pcg_parameters.energy_tol, pcg_parameters.line_check_jump);
            results.energyvecs{i+1} = energyvec;
            results.anglesvecs{i+1} = anglesvec;
        else
            [p, flag, rnorm, iterations, resvec] = pcg(H, -1.0 * grad, pcg_parameters.tol, pcg_parameters.maxit, L, U, x0);
        end
            
        results.num_iter(i+1) = length(resvec);
        results.rnorm(i+1) = rnorm;
        results.relative_rnorm(i+1) = rnorm / norm(-1.0 * grad);
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

    energy = energy_value(un);
    un = line_check_search(p, u, grad);
    
    disp(['= Newton ', num2str(i), ' => ', ...
        'num iter: ', num2str(iterations), ' normed_res: ', num2str(rnorm), ' energy(un): ', num2str(energy), ...
        ' pcg flag: ', num2str(flag)])
    
    if stop_check(un, u, grad)  
        break;
    end  
    
    u = un;
    b_vector(i+1, :) = -1.0 * grad;
    energy_vector(i+1) = energy;
    eta_vector(i+1) = rnorm; 

end

energy_vector = energy_vector(1:(i+1), :);
eta_vector = eta_vector(1:(i+1), :);
b_vector = b_vector(1:(i+1), :);
final_angle_from_grad = final_angle_from_grad(1:(i+1), :);

results.energies = energy_vector;
results.etas = eta_vector;
results.bs = b_vector(1:(i+1), :); % lhs in Ax = b
results.final_angle_from_grad = final_angle_from_grad;

if not(use_direct)
   results.num_iter = results.num_iter(1:(i+1), :);
   results.resvecs = results.resvecs(1:(i+1), :);
   results.normalized_resvecs = results.normalized_resvecs(1:(i+1), :);
   results.rnorm = results.rnorm(1:(i+1), :);
   results.relative_rnorm = results.relative_rnorm(1:(i+1), :);
   
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


function [ results ] = newton_solver( un, preconditioner, pcg_parameters, use_direct, kernel, param_group_id, make_plots)
%NEWTON_SOLVER Summary of this function goes here
%   Detailed explanation goes here

u = un;
p = 0 * u;

L = get_precomputed_split_info();

order = get_ordering(u);

'newton'

x0 = zeros(length(u), 1);

results = struct();

results.final_angle_from_grad = zeros(2000, length(u));
results.grads_pre_ls = zeros(length(u), 2000);
results.b_vector = zeros(length(u), 2000);
results.search_directions = zeros(length(u), 2000);
results.hessians = cell(2000, 1);
results.hessians_pre_ls = cell(2000, 1);
results.guesses = zeros(20, length(u));
results.guesses(1, :) = x0; 
results.energies = zeros(2000, 1);
results.energies(1) = energy_value(un);

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
    
    
    [ ...
        results.grads_pre_ls(:, i+1), ...
        results.hessians_pre_ls{i+1} ...
    ] = grad_hessian_function(u + p, 0);
    results.search_directions(:, i+1) = p;
    results.hessians{i + 1} = H;
    
    [un, alp, balp] = line_check_search(p, u, grad);
    energy = energy_value(un);
    
    results.b_vector(:, i+1) = -1.0 * grad;
    results.energies(i+2) = energy;
    results.guesses(i+2, :) = un;

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

results.final_angle_from_grad = results.final_angle_from_grad(1:(i+1), :);

results.energies = results.energies(1:(i+2), :);
results.bs = results.b_vector(1:(i+1), :); % lhs in Ax = b
results.guesses = results.guesses(1:(i+2), :);
results.bs = results.b_vector(:, 1:(i+1)); % lhs in Ax = b
results.hessians_pre_ls = results.hessians_pre_ls(1:(i+1));
results.grads_pre_ls = results.grads_pre_ls(:, 1:(i+1));
results.search_directions = results.search_directions(:, 1:i+1);

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

end


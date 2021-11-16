function [ energy_vector, eta_vector ] = newton_solver( un, preconditioner, pcg_parameters, use_direct, kernel, param_group_id, make_plots)
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


if make_plots
    figure; pcg_axes = axes; hold(pcg_axes, "on"); 
    set(pcg_axes, 'YScale', 'log');
    title(pcg_axes, strcat('PCG Curves. kernel: ', kernel, '. Param Group:', num2str(param_group_id), '.'))
    ylabel(pcg_axes, 'normed res (log)');
    xlabel(pcg_axes, 'iteration');
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
        rnorm = norm((H * p) + grad);
    else   
        % Preconditioned CG 

        % Warm start
        if pcg_parameters.use_warm_start
            x0 = p;
        else
            x0 = 0 * p;
        end 

        [p, flag, rnorm, iterations, resvec] = pcg_copy(H, -1.0 * grad, pcg_parameters.tol, pcg_parameters.maxit, L,U, x0, u, pcg_parameters.energy_tol, pcg_parameters.line_check_jump);
    end 
    
    
    % Plot PCG residual progress
    if ~use_direct && make_plots
        plot(pcg_axes, resvec / norm(-1.0 * grad), 'DisplayName', ['iter ', num2str(i)]);
    end

    energy = energy_value(un);
    un = line_check_search(p, u, grad);
        
    disp(['= Newton ', num2str(i), ' => ', ...
        'num iter: ', num2str(iterations), ' normed_res: ', num2str(rnorm), ' energy(un): ', num2str(energy), ...
        ' pcg flag: ', num2str(flag)])
    
    if stop_check(un, u, grad)  
        break;
    end  
    
    u = un;
    energy_vector(i+1) = energy;
    eta_vector(i+1) = rnorm; 

end

energy_vector = energy_vector(1:(i+1), :);
eta_vector = eta_vector(1:(i+1), :);

if make_plots
    legend(pcg_axes);
    hold(pcg_axes, "off");
    saveas(pcg_axes, strcat(kernel, '_', num2str(param_group_id), '_pcg.fig'))
end

end


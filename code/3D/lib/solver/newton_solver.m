function [ ] = newton_solver( un, preconditioner, pcg_parameters, use_direct, kernel)
%NEWTON_SOLVER Summary of this function goes here
%   Detailed explanation goes here

u = un;
p = 0 * u;

L = get_precomputed_split_info();

order = get_ordering(u);

'newton'

x0 = zeros(length(u), 1);

pcg_figure = figure;

hold on
figure(pcg_figure);
set(gca, 'YScale', 'log');
title(strcat(kernel, ': PCG Curves'))
ylabel('normed res (log)');
xlabel('iteration');
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
    else   
        % Preconditioned CG 

        % Warm start
        if pcg_parameters.use_warm_start
            x0 = p;
        else
            x0 = 0 * p;
        end 

        [p, d, rnorm, iterations, resvec] = pcg_copy(H, -1.0 * grad, pcg_parameters.tol, pcg_parameters.maxit, L,U, x0);
    end 
    
    
    if ~use_direct
        % Plot PCG residual progress
        figure(pcg_figure);
        plot(resvec / norm(-1.0 * grad), 'DisplayName', ['iter ', num2str(i)]);
    end

    un = line_check_search(p, u, grad);
        
    disp(['= Newton ', num2str(i), ' => ', ...
        'num iter: ', num2str(iterations), ' res: ', num2str(rnorm), ' un: ', num2str(energy_value(un))])
    
    if stop_check(un, u, grad)  
        break;  
    end  
    
    u = un;

end

figure(pcg_figure);
legend();
hold off

saveas(pcg_figure, strcat(kernel, '.fig'))
end


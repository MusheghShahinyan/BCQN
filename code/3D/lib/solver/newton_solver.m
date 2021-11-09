function [ ] = newton_solver( un )
%NEWTON_SOLVER Summary of this function goes here
%   Detailed explanation goes here

u = un;
p = 0 * u;

L = get_precomputed_split_info();

order = get_ordering(u);

'newton'

x0 = zeros(length(u), 1);

hold on
for i = 0 : 2000
    
    [grad, H] = grad_hessian_function(u, 0);

% Incomplete LU factorization, preconditioner
    [L,U] = ilu(H, struct('type','nofill'));

% Diagonal Preconditioner
    % M = spdiags(1./diag(H), 0, length(H),length(H));
     
% Preconditioned CG

    tol = 1e-4;
    restart = 100;
    maxit = 100;
    iterations = 0;
    rnorm = 0;
    
    [p, d, rnorm, iterations, resvec] = pcg(H, -1.0 * grad, tol, maxit, L,U);
     
% Plot PCG residual progress 
    plot(resvec / norm(-1.0 * grad), 'DisplayName', ['iter ', num2str(i)]);
    set(gca, 'YScale', 'log')
     
% Original Implementation
%      p(order) = -1.0 * H(order, order) \ grad(order);
    
    un = line_check_search(p, u, grad);
    
% Warm start
    x0 = p;
    
    disp(['= Newton ', num2str(i), ' => ', ...
        'num iter: ', num2str(iterations), ' res: ', num2str(rnorm), ' un: ', num2str(energy_value(un))])
    
    if stop_check(un, u, grad)
        
        break;
        
    end  
    
    u = un;

end

legend();
end


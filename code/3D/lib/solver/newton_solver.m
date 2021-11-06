function [ ] = newton_solver( un )
%NEWTON_SOLVER Summary of this function goes here
%   Detailed explanation goes here

u = un;
p = 0 * u;

L = get_precomputed_split_info();

order = get_ordering(u);

'newton'

for i = 0 : 2000
 
    
    %plot_result( u, i )

    [grad, H] = grad_hessian_function(u, 0);

    x0 = zeros(length(grad), 1);
    [p, iterations, rnorm] = diag_conjgrad(H, grad, x0);
    p(order) = -1.0 * H(order, order) \ grad(order);

    un = line_check_search(p, u, grad);
    
    disp(['= Newton ', num2str(i), ' => ', ...
        'num iter: ', num2str(iterations), ' res: ', num2str(rnorm)])
    
    if stop_check(un, u, grad)
        
        break;
        
    end  
    
    u = un;

end

end


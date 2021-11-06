function [x, i, rnorm] = diag_conjgrad(A, b, x)
    % Based on
    %  slide 28
    %  https://stanford.edu/class/ee364b/lectures/conj_grad_slides.pdf
    
    r = b - A * x;
    M = diag(1./diag(A)); % inverse of the diagonal of A
    z = M*r;
    
    p = r;
    rho_k = r.' * z;
    
    %disp(['0: residual ', num2str(sqrt(r' * r))]) 
    for i = 1:length(b)
        
        w = A * p;
        alpha = rho_k / (w.' * p);
        x = x + alpha * p;
        r = r - alpha * w;
        
        rnorm = sqrt(r.' * r);
        
        disp([num2str(i), ': residual ', num2str(rnorm)]) 
        
        if rnorm < 1e-5
            break;
        end
        
        z = M*r;
        rho_k_plus_1 = z.' * r;
        p = z + (rho_k_plus_1 / rho_k) * p;
        
        rho_k = rho_k_plus_1;
    end
end
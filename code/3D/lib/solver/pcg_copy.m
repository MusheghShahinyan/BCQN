function [x, flag, relres, iter, resvec, energyvec, anglesvec, xvec, gradvec, estgradvec] = ...
    pcg_copy(A,b, tol, maxit, M1, M2, x0, line_check_u, custom_params, opts1, opts2, varargin)
%PCG   Preconditioned Conjugate Gradients Method.
%   X = PCG(A,B) attempts to solve the system of linear equations A*X=B for
%   X. The N-by-N coefficient matrix A must be symmetric and positive
%   definite and the right hand side column vector B must have length N.
%
%   X = PCG(AFUN,B) accepts a function handle AFUN instead of the matrix A.
%   AFUN(X) accepts a vector input X and returns the matrix-vector product
%   A*X. In all of the following syntaxes, you can replace A by AFUN.
%
%   X = PCG(A,B,TOL) specifies the tolerance of the method. If TOL is []
%   then PCG uses the default, 1e-6.
%
%   X = PCG(A,B,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then PCG uses the default, min(N,20).
%
%   X = PCG(A,B,TOL,MAXIT,M) and X = PCG(A,B,TOL,MAXIT,M1,M2) use symmetric
%   positive definite preconditioner M or M=M1*M2 and effectively solve the
%   system inv(M)*A*X = inv(M)*B for X. If M is [] then a preconditioner
%   is not applied. M may be a function handle MFUN returning M\X.
%
%   X = PCG(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0 is
%   [] then PCG uses the default, an all zero vector (or a random vector if
%   B is empty).
%
%   X = PCG(A,B,TOL,MAXIT,M1,M2,X0,'flex') changes the stadard PCG into the
%   flexibble PCG. The latter is a bit more expensive, but it works in some
%   cases where the standard PCG fails, e.g., if M is not fixed symmetric
%   positive definite.
%
%   X = PCG(A,B,TOL,MAXIT,M1,M2,X0,'null') can be used, if B is a
%   zero vector, to force the PCG to attempt to calculate the nontrivial
%   solution of the homegeneous system of linear equations A*X=0, where A
%   must be symmetric and positive semi-definite. Without the 'null' option,
%   if B is zero, the code immideately returns the trivial solution X=0.
%   The 'flex' and 'null' options can be used together.
%
%   [X,FLAG] = PCG(A,B,...) also returns a convergence FLAG:
%    0 PCG converged to the desired tolerance TOL within MAXIT iterations
%    1 PCG iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 PCG stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during PCG became too
%      small or too large to continue computing.
%
%   [X,FLAG,RELRES] = PCG(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If B is zero, then RELRES is set to 0, unless
%   the 'null' option is set, where the relative residual is defined as
%   NORM(A*X)/NORM(X). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = PCG(A,B,...) also returns the iteration number
%   at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = PCG(A,B,...) also returns a vector of the
%   estimated residual norms at each iteration including NORM(B-A*X0) if
%   the 'null' option in absent. With the 'null' option, the residuals and
%   the relative residuals are defined to be the same.
%
% % Example:
% clear all; n = 100; A = spdiags([1:n]',0,n,n);
% tol = 1e-6;  maxit = 20;
% M = A;
% A=A+sprandsym(n,.1); b = A*ones(n,1);
% x = pcg(A,b,tol,maxit,M);
% % In the last line, one can use a matrix-vector product function as well:
% afun = @(x)A*x; x = pcg(@(x)afun(x),b,tol,maxit,M);
%
% %  Example of the 'flex' option is the same as above, except for M:
% clear all; n = 100; A = spdiags([1:n]',0,n,n);
% tol = 1e-6;  maxit = 20;
% M = A; M(1,2) = 2; % M is no longer symmetric 
% A=A+sprandsym(n,.1); b = A*ones(n,1);
% [x,~,~,~,rv] = pcg(A,b,tol,maxit,M);
% [xf,~,~,~,rvf] = pcg(A,b,tol,maxit,M,[],[],'flex');
% semilogy(rv); hold on; semilogy(rvf,'--rs');
% title('Standard vs. Flexible PCG. Nonsymmetric preconditioning.')
%
% %  Example of the 'null' option with and without the 'flex' option:
% clear all; close all; n = 100; b = zeros(n,1); 
% v = (0:n-1)'; A = spdiags([v 2*v v-1],-1:1,n,n); % symmetric semidefinite
% % A has a null-space spanned by the first coordinate vector 
% tol = 1e-6;  maxit = 50; v = (1:n)'; 
% M = spdiags([v 2*v v-1],-1:1,n,n); % M is SPD
% [xn,~,~,~,rvn] = pcg(A,b,tol,maxit,M,[],[],'null');
% [xnf,~,~,~,rvnf] = pcg(A,b,tol,maxit,M,[],[],'null','flex');
% figure(665); semilogy(rvn); hold on; semilogy(rvnf,'--rs'); 
% title('Standard vs. Flexible PCG null. SPD preconditioning.')
% M = spdiags([v 2*v v],-1:1,n,n);   % M is non-symmetric
% [xn,~,~,~,rvn] = pcg(A,b,tol,maxit,M,[],[],'null');
% [xnf,~,~,~,rvnf] = pcg(A,b,tol,maxit,M,[],[],'null','flex');
% figure(667); semilogy(rvn); hold on; semilogy(rvnf,'--rs');
% title('Standard vs. Flexible PCG null. Nonsymmetric preconditioning.')
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, MINRES, QMR,
%   SYMMLQ, TFQMR, CHOLINC, FUNCTION_HANDLE.
% This is an updated for MATLAB R2011a Revision 2.0 of the code
% http://www.mathworks.com/matlabcentral/fileexchange/50-pcgnull-m
% by Andrew Knyazev, andrew.knyazev@na-net.ornl.gov
% Updated by Andrew Knyazev to implement the new 'flex' option and to make 
% the code PCG compatible so that it can be used as a PCG replacement. 
% Also included new examples in the header.  
%
% This Revision 2.0 is subject to the conditiones of original Revision 1.0:
% "This is a modified version of PCG, Revision 1.6, 1996, by Penny Anderson,
% modified with the permission of The MathWorks, Inc., the copyright owner.
% This Revision 1.0 may not be used with any products other than products of
% The MathWorks, Inc., nor may it be used in or as part of another computer program.
% MATLAB is a registered trademark of The MathWorks, Inc."
if (nargin < 2)
    error('MATLAB:pcg:NotEnoughInputs', 'Not enough input arguments.');
end

if (nargin >= 9)
    if isfield(custom_params, 'energy_tol') && isfield(custom_params, 'line_check_jump')
        energy_tol = custom_params.energy_tol;
        line_check_jump = custom_params.line_check_jump;
        check_energy = true;
    else
        check_energy = false;
    end 

    allow_negative_energy_delta = isfield(custom_params, 'allow_negative_energy_delta') ...
        && custom_params.allow_negative_energy_delta;
    
    store_grad = isfield(custom_params, 'calc_grad') ...
        && custom_params.calc_grad;
    
    check_stopping_pairs = isfield(custom_params, 'stopping_pairs');
        
    check_grad = isfield(custom_params, 'check_grad') ...
        && custom_params.check_grad;
    
    sgd_fallback = isfield(custom_params, 'sgd_fallback') ...
        && custom_params.sgd_fallback;
    
    running_approximation_score = isfield(custom_params, 'running_approximation_score') ...
        && custom_params.running_approximation_score;
end



% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        error('MATLAB:pcg:NonSquareMatrix', 'Matrix must be square.');
    end
    if ~isequal(size(b),[m,1])
        error('MATLAB:pcg:RSHsizeMatchCoeffMatrix', ...
            ['Right hand side must be a column vector of' ...
            ' length %d to match the coefficient matrix.'],m);
    end
else
    m = size(b,1);
    n = m;
    if ~iscolumn(b)
        error('MATLAB:pcg:RSHnotColumn',...
            'Right hand side must be a column vector.');
    end
end
% Assign default values to unspecified parameters
if (nargin < 3) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol <= eps
    warning('MATLAB:pcg:tooSmallTolerance', ...
        strcat('Input tol is smaller than eps and may not be achieved',...
        ' by PCG\n','         Try to use a bigger tolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning('MATLAB:pcg:tooBigTolerance', ...
        strcat('Input tol is bigger than 1 \n',...
        '         Try to use a smaller tolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 4) || isempty(maxit)
    maxit = min(n,20);
end
if ((nargin >= 5) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[m,m])
            error('MATLAB:pcg:WrongPrecondSize', ...
                ['Preconditioner must be a square matrix' ...
                ' of size %d to match the problem size.'],m);
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end
if ((nargin >= 6) && ~isempty(M2))
    existM2 = 1;
    [m2type,m2fun,m2fcnstr] = iterchk(M2);
    if strcmp(m2type,'matrix')
        if ~isequal(size(M2),[m,m])
            error('MATLAB:pcg:WrongPrecondSize', ...
                ['Preconditioner must be a square matrix' ...
                ' of size %d to match the problem size.'],m);
        end
    end
else
    existM2 = 0;
    m2type = 'matrix';
end
flexible = 0; nullPCG=0; %the default
if ((nargin >= 11))
    if strcmp(opts1,'flex')
        flexible = 1; %Flexible PCG
    elseif strcmp(opts1,'null')
        nullPCG = 1; %nullPCG
    end
end
if ((nargin >= 11))
    if strcmp(opts2,'flex')
        flexible = 1; %Flexible PCG
    elseif strcmp(opts2,'null')
        nullPCG = 1; %nullPCG
    end
end
if ((nargin >= 7) && ~isempty(x0))
    if ~isequal(size(x0),[n,1])
        error('MATLAB:pcg:WrongInitGuessSize', ...
            ['Initial guess must be a column vector of' ...
            ' length %d to match the problem size.'],n);
    elseif norm(x0)==0 && nullPCG
        error('MATLAB:pcg:WrongInitGuessSizeNullPCG', ...
            ['Initial nullPCG guess must be a nonzero vector.'],n);
    else
        if ~nullPCG
            x = x0;
        else
            x = x0/norm(x0); %normalize nullPCG initial guess
        end
    end
else
    if ~nullPCG
        x = zeros(n,1);
    else
        x0 = randn(n,1); %nullPCG does not allow zero initial guess
        x = x0/norm(x0); %better make it normalized
    end
end
if ((nargin > 11) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error('MATLAB:pcg:TooManyInputs', 'Too many input arguments.');
end
% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % Norm of rhs vector, b
if (n2b == 0) && ~nullPCG          % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    iter = 0;                      % no iterations need be performed
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('pcg',tol,maxit,0,flag,iter,NaN);
    end
    return
end
% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
if ~nullPCG
    tolb = tol * n2b;                  % Relative tolerance
else
    tolb = tol;
end
r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
rp = zeros(size(r));
normr = norm(r); % Norm of residual
normr_act = normr;
if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    if ~nullPCG
        relres = normr / n2b;
    else
        relres = normr;
    end
    iter = 0;
    resvec = normr;
    if (nargout < 2)
        itermsg('pcg',tol,maxit,0,flag,iter,relres);
    end
    return
end

energyvec = zeros(maxit+1,1);
anglesvec = zeros(maxit+1,1);
xvec = zeros(length(x), maxit+1);
estgradvec = zeros(maxit+1, length(x)); % TODO be consistent with xvec and store as columns

grad_ii = 1;
grads_for_checking = zeros(maxit+1, 1);
gradnorm = zeros(maxit+1,1);
approxdelta = zeros(maxit+1,1);
resvec = zeros(maxit+1,1);         % Preallocate vector for norm of residuals
resvec(1,:) = normr;               % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of minimum residual
rho = 1;
stag = 0;                          % stagnation of the method
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;

un = line_check_search(b, line_check_u, -1.0 * b);
energyvec(1) = energy_value(un);
anglesvec(1) = 0;

prev_energy = energy_value(un);
prev_x = x;

xvec(:, 1) = x;
if store_grad || check_grad
   grad = grad_function(un); 
end

if store_grad
    [gradvec(1, :)] = grad;
    estgradvec(1, :) = b;
end

if check_grad
    gradnorm(1) = norm(grad) / norm(b);
    approxdelta(1) = 0; % newton approximation of gradient is perfect at the current location
    prev_grad_mean = gradnorm(1);
end

fprintf('\n');
logstringlen = 0;

% loop over maxit iterations (unless convergence or failure)
for ii = 1 : maxit
    if existM1
        y = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
        if ~all(isfinite(y))
            flag = 2;
            break
        end
    else % no preconditioner
        y = r;
    end
    
    if existM2
        z = iterapp('mldivide',m2fun,m2type,m2fcnstr,y,varargin{:});
        if ~all(isfinite(z))
            flag = 2;
            break
        end
    else % no preconditioner
        z = y;
    end
    
    rho1 = rho;
    rho = r' * z;
    if ((rho == 0) || isinf(rho))
        flag = 4;
        break
    end
    if (ii == 1)
        p = z;
    else
        if flexible
            beta = (r-rp)' * z / rho1;  %flexible
        else
            beta = rho / rho1; %standard
        end
        if ((beta == 0) || isinf(beta))
            flag = 4;
            break
        end
        p = z + beta * p;
    end
    q = iterapp('mtimes',afun,atype,afcnstr,p,varargin{:});
    pq = p' * q;
    if ((pq <= 0) || isinf(pq))
        flag = 4;
        break
    else
        alpha = rho / pq;
    end
    if isinf(alpha)
        flag = 4;
        break
    end
    
    % Check for stagnation of the method
    if (norm(p)*abs(alpha) < eps*norm(x))
        stag = stag + 1;
    else
        stag = 0;
    end
    
    x = x + alpha * p;             % form new iterate
    if flexible
        rp = r;
    end
    r = r - alpha * q;
    if ~nullPCG
        normr = norm(r);
    else
        normr = norm(r)/norm(x);
    end
    normr_act = normr;
    resvec(ii+1,1) = normr;
    
    un = line_check_search(x, line_check_u, -1.0 * b);
    energyvec(ii+1) = energy_value(un);
    
    CosTheta = max(min(dot(x,b)/(norm(x)*norm(b)),1),-1);
    anglesvec(ii+1) = real(acosd(CosTheta));
    
    xvec(:, ii+1) = x;
    
    if store_grad || check_grad
        grad = grad_function(un);
    end
    
    if store_grad
       [gradvec(ii+1, :)] = grad;
       estgradvec(ii+1, :) = r;
    end
    
    logstring = strcat("cg iter ", num2str(ii));
    
    % gradient stopping condition
    if running_approximation_score
        gradnorm(ii+1) = norm(grad) / norm(b);
        approxdelta(ii + 1) = norm(grad - r) / norm(b);
        
        % moving average window of width 3
        est = mean(approxdelta(max(1, ii-3):ii+1));        
        
        logstring = strcat("cg iter ", num2str(ii), " est ",...
                            num2str(est), " maxit ", num2str((1 - est) * maxit), ...
                            " thresh ", num2str((1 + (0.001/tolb) * est) * tolb),...
                            " normr ", num2str(normr));
        
        if ii >= 20 && (normr / norm(b)) < 0.01
            if sgd_fallback && est > 2
                flag = 11;
                xmin = b;
                imin = 0;
                break;
            end
            
            if ii >= (1 - est) * maxit 
                flag = 9;
                
                if (normr < normrmin)      % update minimal norm quantities
                    xmin = x;
                    imin = ii;
                end
                
                break;
            end
            
            if normr < (1 + (0.001/tolb) * est) * tolb
                flag = 10;
                
                if (normr < normrmin)      % update minimal norm quantities
                    xmin = x;
                    imin = ii;
                end
                
                break;
            end
        end
    end 
    
    if check_grad 
        gradnorm(ii+1) = norm(grad) / norm(b);
        
        if mod(ii, custom_params.check_grad_jump) == 0 
            grads_for_checking(grad_ii) = gradnorm(ii+1);
            grad_ii = grad_ii + 1;

            % moving average window of width 5
            %grad_mean = mean(gradnorm(max(1, ii-5):ii+1));
            grad_mean = mean(grads_for_checking(max(1, grad_ii-7):(grad_ii - 1)));

            % We must have atleast 10 iterations so that we don't trigger
            %  on an initial hump, also make sure we actually make progress 
            %  (the norm derease below the starting norm)
            if ii >= 10 && grad_mean > prev_grad_mean && grad_mean < gradnorm(1)
                [~, idx] = min(gradnorm(1:ii+1));
                xmin = xvec(:, idx);
                imin = idx - 1;

                flag = 8;
                break;
            end
         
            prev_grad_mean = grad_mean;
        end
    end
    
    
    % stopping pairs stopping condition
    if check_stopping_pairs
        stop = false;
        for pair_idx = 1:size(custom_params.stopping_pairs, 1)
            if (custom_params.stopping_pairs(pair_idx, 1) == custom_params.newton_iter ...
                  && custom_params.stopping_pairs(pair_idx, 2) == ii)
                stop = true;
                break;
            end
        end
        
        if stop
            flag = 7;
            xmin = x;
            imin = ii;
            break;
        end
    end
     
    % line search energy stopping condition, check outer line search 
    if check_energy 
        if mod(ii, line_check_jump) == 0
            un = line_check_search(x, line_check_u, -1.0 * b);
            energy = energy_value(un);

            energy_delta = (prev_energy - energy) / prev_energy;
            if not(allow_negative_energy_delta) && (energy_delta < 0)
                flag = 5;
                x = prev_x;
                xmin = x;
                imin = ii - line_check_jump;
                break;
            else
                if allow_negative_energy_delta
                    energy_delta = abs(energy_delta);
                end 
                
                if (energy_delta < energy_tol)
                    flag = 6;
                    break;
                end 
            end
            prev_x = x;
            prev_energy = energy;
        end
    end
    
    fprintf(repmat('\b', 1, logstringlen));
    logstringlen = strlength(logstring);
    fprintf(logstring)
    
    
        
    % check for convergence (default implementation)
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
        normr_act = norm(r);
        resvec(ii+1,1) = normr_act;
        if (normr_act <= tolb)
            flag = 0;
            iter = ii;
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning('MATLAB:pcg:tooSmallTolerance', ...
                        strcat('Input tol may be smaller than eps*cond(A)',...
                        ' and might not be achieved by PCG\n',...
                        '         Try to use a bigger tolerance'));
                end
                flag = 3;
                iter = ii;
                break;
            end
        end
    end
    if (normr_act < normrmin)      % update minimal norm quantities
        normrmin = normr_act;
        xmin = x;
        imin = ii;
    end
    if stag >= maxstagsteps
        flag = 3;
        break;
    end
end                                % for ii = 1 : maxit
% returned solution is first with minimal residual
if (flag == 0)
    if ~nullPCG
        relres = normr_act / n2b;
    else
        relres = normr_act;
    end
else
    r_comp = b - iterapp('mtimes',afun,atype,afcnstr,xmin,varargin{:});
    if norm(r_comp) <= normr_act
        x = xmin;
        iter = imin;
        if ~nullPCG
            relres = norm(r_comp) / n2b;
        else
            relres = norm(r_comp);
        end
    else
        iter = ii;
        if ~nullPCG
            relres = normr_act / n2b;
        else
            relres = normr_act;
        end
    end
end
% truncate the zeros from resvec and energyvec
if ((flag <= 1) || (flag == 3))
    resvec = resvec(1:ii+1,:);
    energyvec = energyvec(1:ii+1,:);
    anglesvec = anglesvec(1:ii+1,:);
    gradvec = gradvec(1:ii+1,:);
    estgradvec = estgradvec(1:ii+1,:);
    xvec = xvec(:, 1:ii+1);
else
    resvec = resvec(1:ii,:);
    energyvec = energyvec(1:ii,:);
    anglesvec = anglesvec(1:ii,:);
    gradvec = gradvec(1:ii+1,:);
    estgradvec = estgradvec(1:ii+1,:);
    xvec = xvec(:, 1:ii+1);
end
% only display a message if the output flag is not used
if (nargout < 2)
    itermsg('pcg',tol,maxit,ii,flag,iter,relres);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                            %%%
%%%      Partial Conjugate Gradient Method     %%%
%%%            (P) minimize z=f(x)             %%%
%%%                                            %%%
%%%    arguments (*required)                   %%%
%%%    *f  - symbolic scalar function          %%%
%%%    *x0 - numeric initial location          %%%
%%%    *m  - restart after m<n iterations      %%%
%%%          m=0 -> gradient descent           %%%
%%%          m=n -> conjugate gradient         %%%
%%%    lineSearch - "A" for Armijo's rule      %%%
%%%                 "W" for Wolfe's rule       %%%
%%%                 "G" for Goldstein's rule   %%%
%%%    method - "P" for Polak-Ribiere          %%%
%%%             "F" for Fletcher-Reeves        %%%
%%%    tolerance - convergence tolerance       %%%
%%%    maxIt - maximum iterations              %%%
%%%                                            %%%
%%%    returns                                 %%%
%%%    x_star - optimal location               %%%
%%%    z_star - optimal objective value        %%%
%%%    x_values - matrix where kth location    %%%
%%%                is kth column               %%%
%%%    iterations - number of iterations       %%%
%%%                                            %%%
%%%    Author: Harrison Watts                  %%%
%%%                                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_star, z_star, x_values, iterations] = pcgrad(f,x0,m,lineSearch,method,tolerance,maxIt)
 
    %%% optional arguments %%%
    if ~exist('lineSearch')  %
        lineSearch = 'A';    %
    end                      %
                             %
    if ~exist('method')      %
        method = 'F';        %
    end                      %
                             %
    if ~exist('tolerance')   %
        tolerance = 1e-12;   %
    end                      %
                             %
    if ~exist('maxIt')       %
        maxIt = 50;          %
    end                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % dimension
    n = size(x0,2);
    
    % initial x vector
    x = zeros(n,m+2);
    
    % initial y vector
    y = zeros(n,2);
    
    % initial location
    y(:,1)=x0;
    
    y(:,2)= y(:,1);
        
    x_values = nan(n,maxIt);
    
    x_values(:,1) = x0;
    
    % initial gradient vector
    g = zeros(n,m+2);
    
    % convert location for symbolic use
    x0c = num2cell(x0);
    
    % gradient 
    grad = gradient(f);
    
    g(:,1) = grad(x0c{:});
    
    % initial step size vector
    alpha = zeros(m+1,1);
    
    % initial beta vector
    beta = alpha;

    % loop until convergence
    converged = 0;
    
    iterations = 1;
    
    while (~converged) && (iterations <=maxIt)
    
        % restart from last step
        x(:,1) = y(:,2);
        
        y(:,1) = x(:,1);
        
        x0c=num2cell(x(:,1));
        
        % restart with steepest descent direction
        d(:,1) = -grad(x0c{:});
        
        % restart after m iterations of conjugate gradient
        for k=1:m+1
            
            if lineSearch == 'G'
                % use Goldstein's rule to line search
                alpha(k) = goldstein(f,x(:,k),d(:,k));
                
            elseif lineSearch == 'W'
                % use Wolfe's rule to line search
                alpha(k) = wolfe(f,x(:,k),d(:,k));
                
            else 
                % use Armijo's rule to line search
                alpha(k) = armijo(f,x(:,k),d(:,k));
                
            end
            
            % new partial location
            x(:,k+1)=x(:,k)+alpha(k)*d(:,k);
            
            xc = num2cell(x(:,k+1));
            
            % gradient at new location
            g(:,k+1) = grad(xc{:});
            
            % Polak-Ribiere
            if method == 'P'
                
                beta(k) = ((g(:,k+1)-g(:,k))'*g(:,k+1))/(g(:,k)'*g(:,k));
                
            % Fletcher-Reeves    
            else
            
                beta(k) = (g(:,k+1)'*g(:,k+1))/(g(:,k)'*g(:,k));
            
            end
            
            x_values(:,iterations) = x(:,m+2);
            
            % new conjugate direction
            d(:,k+1)=double(-g(:,k+1)+beta(k)*d(:,k));
            
            
        end

        % new location
        y(:,2) = double(x(:,m+2));
        
        %Cauchy criterion
        converged = norm(y(:,2)-y(:,1))<=tolerance;
        
        iterations = iterations+1;
    
    end
    
    % optimal location
    y_star = num2cell(x(:,m+2));
    
    x_star = double(x(:,m+2));
    
    % optimal objective function value
    z_star = double(f(y_star{:}));
    
    % all location values
    no_index = isnan(x_values)==zeros(n,1);
    x_values = x_values(:,no_index(1,:));
    
    % number of iterations
    iterations = iterations;


end
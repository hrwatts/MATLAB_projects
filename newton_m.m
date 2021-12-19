%% Newton's method
% solve an equation f(x) = 0
% x0 = initial guess
% number of iterations
function [solution, x] = newton_m(f,fprime,x0,N)
    x = x0;
    for n = 2:N
        x(n)=x(n-1)-f(x(n-1))./fprime(x(n-1));
        if abs(x(n)-x(n-1))<=2.2204e-16
            break;
        end
    solution = x(end);
end
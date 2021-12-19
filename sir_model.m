%% sir model ode
% x = [s, i, r]
% a = contact rate
% k = recovery rate
function dxdt = sir_model(t,x,a,k)
    dxdt = [-a*x(1)*x(2); a*x(1)*x(2)-k*x(2); k*x(2)];
end
%% Fit SIR curve
% data = dates,cases,deaths in .csv format
% N = population
function [tspan, curve, k, a] = sir_fit(file,N)
% weekly data
weekdata=week_tab(file);

% 35 initial cases (week 9)
I0 = weekdata.cases(weekdata.cases>1);
I0 = I0(1);
S0 = N-I0;

% peak
[Ihat that] =max(weekdata.cases);


%% Basic Reproduction Number (R0)
% proportions
i0 = I0/N; s0 = S0/N; ihat = Ihat/N;
% basic reproduction number
syms R
% integration constant
K = @(R) i0+s0-R*log(s0); 

% Newton's Method
% i(R0)
iR0 = @(R) -R+R.*log(R)+K(R)-ihat;

% i'(R0)
ipR0 = matlabFunction(diff(iR0(R)));

% initial guess R00
R00 = .5;
it = 10; %iterations
R0 = newton_m(iR0,ipR0,R00,it);

% read dates
dates = readvars(file);
ndays = length(dates);
week = weekdata.week;
tspan = linspace(0,week(end),ndays);
tspan1 = linspace(0,week(end)+12,ndays+12*7);

%% recovery time: quarter-day to 30 days
that1 = []; a = []; k = []; index =0;
    for k1 = 14:-.125:.025
        index=index+1;
        k(index) = 7/k1;
        a(index) = k(index)/R0;

        % solve using ode45 % V = [s, i, r]
        [~,V] = ode45(@(t,x) sir_model(t,x,a(index),k(index)), tspan1, [s0; i0; 0]);
        s = V(:,1); i = V(:,2); r = V(:,3);
        [~, that1(index)]=max(i);
    end
    [~,that1]= min(abs(that-that1./7));
    k = k(that1); a = a(that1);
    [tspan,curve] = ode45(@(t,x) sir_model(t,x,a,k), tspan, [s0; i0; 0]);
    tspan = tspan; curve = curve*N;
end
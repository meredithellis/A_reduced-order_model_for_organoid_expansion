function [t, c, w] = metabolite_concentration_LW_N0_flow(proliferation_rate, uptake_rate, N0, flow_scale, W_tol_dim)
% (proliferation_rate, uptake_rate) where proliferation rate is either
% 'low','mid','high' or can take the number of days over which the cells
% double

% timescale is 24 hours
%% CXP1 parameters
% CXP1_parameters_TDgrowth_casestudy('proliferation_rate', 'uptake_rate')

if nargin<5
    W_tol_dim=0.252;
end

[timescale, Lengthscale, flow_scale, metab_coeff, mu, Ubar, Hh, Kc,Kw, rho, sigma,dcm,dch,dwm,dwh, lambdaM, lambdaH, C_uptake_rate,p,W_tol] = ...
            CXP1_parameters_TDgrowth_casestudy_vary_n0_flow(proliferation_rate, uptake_rate, N0, flow_scale,W_tol_dim);


%% problem parameters

alphaC =1+ Kc*(Hh/(1-Hh));
alphaW=1+ Kw*(Hh/(1-Hh));
beta = mu*Ubar;
dc = dcm+ dch*Kc*(Hh/(1-Hh));
dw=dwm + dwh*Kw*(Hh/(1-Hh));
lambda=lambdaM+ lambdaH*Kw*(Hh/(1-Hh));

gamma = rho*Kc*(Hh/(1-Hh));
eta= sigma*Kw*(Hh/(1-Hh));

%% Problem set-up
% Create an interval of the space domain...
dom = [0 1];
%...and specify a sampling of the time domain:
timestep=0.05/24;
t = 0:timestep:7;

% Make the right-hand side of the PDE.
pdefun = @(t,x,c,w) [(-beta*diff(c)+dc*diff(c,2)-gamma*c*exp(p*t))/alphaC ; (-beta*diff(w)+dw*diff(w,2)+eta*c*exp(p*t)-lambda*w)/alphaW];

% Assign boundary conditions.
bc.left = @(t,c,w) [beta*c-dc*diff(c)-beta; beta*w-dw*diff(w)];
bc.right = @(t,c,w) [-dc*diff(c); -dw*diff(w)];

% Construct a chebfun of the space variable on the domain,
x = chebfun(@(x) x, dom);
% and of the initial conditions.
c0 = chebfun(1/alphaC,dom);
%c0 = chebfun(c0,"trunc",30);
w0 = chebfun(0,dom);
sol0 = [c0, w0];

%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-6, 'Ylim', [0.5,1],"AdjustBCs", false);

%% Call pde15s to solve the problem.
tic
[t, c, w] = pde15s(pdefun, t, sol0, bc, opts);
toc


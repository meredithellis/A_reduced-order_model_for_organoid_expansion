function [t, X, c, w] = metabolite_concentration_SL(proliferation_rate, uptake_rate)
% (proliferation_rate, uptake_rate) where proliferation rate is either
% 'low','mid','high' or can take the number of days over which the cells
% double

% timescale is 24 hours
%% CXP1 parameters

[timescale, Lengthscale, flow_scale, metab_coeff, mu, Ubar, Hh, Kc,Kw, rho, sigma,dcm,dch,dwm,dwh, lambdaM, lambdaH, C_uptake_rate,p] = ...
            CXP1_parameters_TDgrowth_SL(proliferation_rate, uptake_rate);


%% problem parameters

alphaC =1+ Kc*(Hh/(1-Hh));
alphaW=1+ Kw*(Hh/(1-Hh));
beta = mu*Ubar;
dc = dcm+ dch*Kc*(Hh/(1-Hh));
dw=dwm + dwh*Kw*(Hh/(1-Hh));
lambda=lambdaM+ lambdaH*Kw*(Hh/(1-Hh));

gamma = rho*Kc*(Hh/(1-Hh));
eta= sigma*Kw*(Hh/(1-Hh));
%%
s_values   = linspace(-7, 1, 801);
tau = linspace(0, 1/beta, 1401); %1/beta
[ss, tautau] = meshgrid(s_values,tau);
%%
Time_span= tau;
% initial condition
W0 = zeros(length(s_values));
W0_r1 = zeros(length(s_values(s_values>=0)),1);
W0_r2 = zeros(length(s_values(s_values<0)),1);

[tau_r1, W_r1] = ode45( @(tau_r1,W_r1) W_r1_cr_rhs(tau_r1,W_r1,s_values,proliferation_rate,uptake_rate), Time_span, W0_r1);
[tau_r2, W_r2] = ode45( @(tau_r2,W_r2) W_r2_cr_rhs(tau_r2,W_r2,s_values,proliferation_rate,uptake_rate), Time_span, W0_r2);

W_combine = [W_r2 W_r1];
%%
for i=1:length(s_values)
    for j= 1:length(tau)
            s_val = s_values(i);
            tau_val = tau(j);
            if s_val >=0 %in region 1
                t_values(i,j) = alphaW*tau_val;
                X_values(i,j) = beta*tau_val + s_val;
                C_solution(i,j) = C_r1_cr(s_val,tau_val, proliferation_rate, uptake_rate);
                
            elseif s_val < 0
                t_values(i,j) = alphaW*tau_val - s_val;
                X_values(i,j) = beta*tau_val ;
                C_solution(i,j) = C_r2_cr(s_val,tau_val, proliferation_rate, uptake_rate);
                
            end
    end
end

%%
t=t_values;
X = X_values;
c = C_solution;
w = W_combine';
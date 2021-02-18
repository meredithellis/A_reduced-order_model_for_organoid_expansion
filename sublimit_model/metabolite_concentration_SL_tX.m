function [C, W] = metabolite_concentration_SL_tX(t, X, proliferation_rate, uptake_rate)

[~, ~, ~, ~, mu, Ubar, Hh, Kc,Kw, rho, sigma,dcm,dch,dwm,dwh, lambdaM, lambdaH, ~,~] = ...
            CXP1_parameters_TDgrowth_SL(proliferation_rate, uptake_rate);


%% problem parameters

alphaC =1+ Kc*(Hh/(1-Hh));
alphaW=1+ Kw*(Hh/(1-Hh));
beta = mu*Ubar;

gamma = rho*Kc*(Hh/(1-Hh));
eta= sigma*Kw*(Hh/(1-Hh));
%% COMPUTE (s,tau) FROM (t,X)
s_val = X - beta*t/alphaC;

if s_val >= 0 % REGION 1
    tau_val = t/alphaC;
    tau_space = linspace(0, tau_val, 1000); 
    W0_r1=0;
    
    if t==0
        W=W0_r1;
        
    else
    [tau_r1, W_r1] = ode45( @(tau_r1,W_r1) W_r1_cr_rhs(tau_r1,W_r1,s_val,proliferation_rate,uptake_rate), tau_space, W0_r1);

    W = W_r1(end); % evaluate at the time required
    end
    C = C_r1_cr(s_val,tau_val, proliferation_rate, uptake_rate);
    
    
elseif s_val < 0 % REGION 2
    s_val = alphaC*X/beta - t;
    tau_val = X/beta;
    tau_space = linspace(0, tau_val, 1001); 
    W0_r2=0;
    
    if X==0
        W=W0_r2;
    else
    [tau_r2, W_r2] = ode45( @(tau_r2,W_r2) W_r2_cr_rhs(tau_r2,W_r2,s_val,proliferation_rate,uptake_rate), tau_space, W0_r2);

    W = W_r2(end); % evaluate at the time required
    end
    C = C_r2_cr(s_val,tau_val, proliferation_rate, uptake_rate);
end

%%

%%


end
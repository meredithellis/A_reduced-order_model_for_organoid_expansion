function dW = W_r1_cr_rhs(tau_r1, W_r1,s_values, proliferation_rate, uptake_rate)

%% problem parameters
[timescale, Lengthscale, flow_scale, metab_coeff, mu, Ubar, Hh, Kc,Kw, rho, sigma,dcm,dch,dwm,dwh, lambdaM, lambdaH, C_uptake_rate,p] = CXP1_parameters_TDgrowth_SL(proliferation_rate, uptake_rate);
alphaC =1+ Kc*(Hh/(1-Hh));
alphaW=1+ Kw*(Hh/(1-Hh));
beta = mu*Ubar;
dc = dcm+ dch*Kc*(Hh/(1-Hh));
dw=dwm + dwh*Kw*(Hh/(1-Hh));
gamma = rho*Kc*(Hh/(1-Hh));
eta= sigma*Kw*(Hh/(1-Hh));
lambda=lambdaM+ lambdaH*Kw*(Hh/(1-Hh));

%%
%s_values =linspace(-beta*7/alphaC, 1, 101);
S_values = s_values(s_values>=0);


for i= 1:length(S_values)
    Si = S_values(i);
    
    dW(i,1) =  eta*C_r1_cr(Si,tau_r1, proliferation_rate, uptake_rate).*exp(p.*alphaW*tau_r1) - lambda*W_r1(i);
end

end
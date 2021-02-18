function f = C_r1_cr(s,tau, proliferation_rate, uptake_rate)

%% CXP1 parameters
% proliferation_rate='mid';
% uptake_rate='mid';
[~, ~, ~, ~, mu, Ubar, Hh, Kc,Kw, rho, sigma,dcm,dch,dwm,dwh, lambdaM, lambdaH, ~,p] = CXP1_parameters_TDgrowth_SL(proliferation_rate, uptake_rate);


%% problem parameters

alphaC =1+ Kc*(Hh/(1-Hh));
alphaW=1+ Kw*(Hh/(1-Hh));
beta = mu*Ubar;
dc = dcm+ dch*Kc*(Hh/(1-Hh));
dw=dwm + dwh*Kw*(Hh/(1-Hh));
gamma = rho*Kc*(Hh/(1-Hh));
eta= sigma*Kw*(Hh/(1-Hh));
lambda=lambdaM+ lambdaH*Kw*(Hh/(1-Hh));

%%


f= (1/alphaC)*exp(gamma/(alphaC*p)).*exp(- (gamma/(alphaC*p)).*exp(p*alphaC*tau))+0.*s;


end
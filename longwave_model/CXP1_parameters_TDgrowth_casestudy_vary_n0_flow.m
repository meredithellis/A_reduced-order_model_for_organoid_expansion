% Call model parameter values for the current CXP1 operating conditions. 
%(Note, am calling parameters used in the model, not the parameters used in non-dimensionalising)

function [timescale, Lengthscale, flow_scale, metab_coeff, mu, Ubar, Hh, Kc,Kw, rho, sigma,dcm,dch,dwm,dwh, lambdaM, lambdaH, nuC,p,W_tol] = ....
    CXP1_parameters_TDgrowth_casestudy_vary_n0_flow(proliferation_rate, uptake_rate, N0, flow_scale,W_tol_dim)
% Call the parameters of the CXP1 system for 'low', 'mid' or 'high',
% (proliferation_rate, uptake_rate)

if nargin<5
    W_tol_dim=0.252;
end

timescale=3600*24; %[s]c x 
%
Lengthscale = (90e-3); % m
%flow_scale = 1e-6; % m per s
metab_coeff = 0.36; % mol per m2

mu=flow_scale*timescale/Lengthscale;
Ubar=1/3;
% hydrogel:media 
Hh=1/3;
% partition coeff
Kc=1;
Kw=1;
% diffusivities
dcm=6e-10*timescale/((Lengthscale)^2);
dch=6e-10*timescale/((Lengthscale)^2);
dwm=1.4e-9*timescale/((Lengthscale)^2);
dwh=1.2e-9*timescale/((Lengthscale)^2);
% natural decay
lambdaM=1e-14*timescale;
lambdaH=1e-14*timescale;
% initial number of cells seeded
% initial cell seeding density
N0; % = 3.3e10; 


if strcmpi(proliferation_rate,'low') == 1 
    p_rate = 1.9e-6; % prolif per second (doubling every 6 days)
elseif strcmpi(proliferation_rate,'mid') == 1 
    p_rate = 3.9e-6; %(doubling every 3 days)
elseif strcmpi(proliferation_rate,'high') == 1 
    p_rate = 1.2e-5; %(doubling every 1 day)
else
    alpha_day = proliferation_rate;
    p_rate = 1/(alpha_day*24*3600); % doubling every alpha_day 
end

if strcmpi(uptake_rate, 'low') == 1
    nuC = 9.4e-18; % uptake rate per second   - m^2 per cell per second
elseif strcmpi(uptake_rate, 'mid') == 1
    nuC = 9.4e-17; % 
elseif strcmpi(uptake_rate, 'high') == 1  
    nuC = 9.4e-16; % 
elseif strcmpi(uptake_rate, 'midx5') ==1
    nuC = 5*9.4e-17;
end

W_tol=W_tol_dim/metab_coeff;
p = p_rate*timescale;
rho=N0*nuC*timescale; 
sigma=2*rho;


%



end
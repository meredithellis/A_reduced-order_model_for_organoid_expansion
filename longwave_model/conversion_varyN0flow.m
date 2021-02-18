function [conversion_nutrient, amount_consumed] = conversion_varyN0flow(t, C, proliferation_rate, uptake_rate, N0, flow_scale, W_tol_dim)

if nargin<7
    W_tol_dim=0.252;
end

[~, ~, ~, ~, mu, Ubar, Hh, ~,~, rho, sigma,dcm,dch,dwm,dwh, ~, ~, ~,p,~] = ....
    CXP1_parameters_TDgrowth_casestudy_vary_n0_flow(proliferation_rate, uptake_rate, N0, flow_scale);

influx = mu*(1-Hh)*Ubar.*ones(1,length(t));
outflux = mu*(1-Hh)*Ubar.*C(1, :);

total_fed_into_system = cumtrapz(t, influx)+2/3;
total_flow_out_of_system =cumtrapz(t, outflux);


left_in_bioreactor = sum(C(:,:)); % calculate the amount of nutrient left in the bioreactor at time t*

amount_consumed= total_fed_into_system - total_flow_out_of_system - left_in_bioreactor;

conversion_nutrient = amount_consumed./total_fed_into_system;

% for plotting later
amount_consumed=amount_consumed';
conversion_nutrient = conversion_nutrient';
end
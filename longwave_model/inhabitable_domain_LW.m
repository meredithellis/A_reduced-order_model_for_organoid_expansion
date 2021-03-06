function [rW, fraction_uninhabitable] = inhabitable_domain_LW(t, W, W_maxtol)
% rW is the position where W=0.7, so X>= rW is above max lactate threshold
if nargin<5
    W_maxtol=0.7;
end

j=0;
%W_maxtol; %=0.7;
%[timescale, Lengthscale, flow_scale, metab_scale, mu, Ubar, Hh, Kc,Kw, rho, sigma,dcm,dch,dwm,dwh, lambdaM, lambdaH, C_uptake_rate,p] = CXP1_parameters_TDgrowth_casestudy(proliferation_rate, uptake_rate);

for i=1:length(t)
    j=j+1;
r =  roots(W(:,i)-W_maxtol);

  if isempty(r)==0   % if there is a root of W-0.7=0
    rW(j,1) = t(i);
    rW(j,2) = min(r);

  elseif max(W(:,i))> W_maxtol % if there is no root, but the max(W)>Wmaxtol, then have W>Wmaxtol
    rW(j,1) = t(i);
    rW(j,2) = 0;
  else                        % no root and max W<Wmaxtol
    rW(j,1) = t(i);
    rW(j,2) = 1;
      
  end
 
  

end

fraction_uninhabitable = ones(length(t),1)- rW(:,2);
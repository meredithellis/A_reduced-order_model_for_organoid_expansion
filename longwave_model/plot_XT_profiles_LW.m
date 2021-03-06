%% Computing metabolite concentrations for longwave approximation model

% DEFINE THE CONTROL PARAMETER VALUES
% N0 is the initial cell seeding density cells m^{-2}
N0=3.3e10;
% flow_scale is the inlet flow rate m/s
flow_scale= 1e-6;
% DEFINE THE CELL TYPE
%proliferation rate 'low' sets P=1/6, 'mid' sets P=1/3, 'high' sets P=1
proliferation_rate='mid'; 
%uptake rate 'low' sets rho=0.027, 'mid' sets rho=0.27, 'high' sets rho=2.7
uptake_rate='mid';
% ratio of maximum tolerated lactate concentration to inlet glucose
% concentration
W_tol_dim=0.7*0.36;
[~, ~, ~, ~, ~, ~, ~, ~,~, ~, ~,~,~,~,~, ~, ~, ~,~,W_tol] = ....
    CXP1_parameters_TDgrowth_casestudy_vary_n0_flow(proliferation_rate, uptake_rate, N0, flow_scale,W_tol_dim);

%% compute C and W for X \in [0,1] and T \in [0,7]
[T, C_LW, W_LW] = metabolite_concentration_LW_N0_flow(proliferation_rate, uptake_rate, N0, flow_scale);

%% plot glucose concentration C(X,T)

surf(C_LW,T)
colormap viridis
colorbar
caxis([0 1])
view([0 90])
title({'Glucose concentration', ['N_0= ' num2str(N0, '%1.1e') ' cell m^{-1}, [u]=' num2str(flow_scale, '%1.e') 'ms^{-1}']},'FontWeight','Normal')
ax = gca;
ax.FontSize = 14;
xlabel('$X$','Interpreter','latex')
ylabel('time, $T$','interpreter','latex')
set(gcf, 'Position',  [100, 100, 600, 500])

%% plot lactate concentration W(X,T)

surf(W_LW,T)
colormap viridis
colorbar
caxis([0 inf]) % we use caxis([0 2.1]) in paper figures
view([0 90])
title({'Lactate concentration', ['N_0= ' num2str(N0, '%1.1e') ' cell m^{-1}, [u]=' num2str(flow_scale, '%1.e') 'ms^{-1}']},'FontWeight','Normal')
ax = gca;
ax.FontSize = 14;
xlabel('$X$','Interpreter','latex')
ylabel('time, $T$','interpreter','latex')
set(gcf, 'Position',  [100, 100, 600, 500])

%% plot fraction of uninhabitable domain
%compute fraction of domain which is uninhabitable as a function of T
[~, fraction_uninhabitable] = inhabitable_domain_LW(T, W_LW, W_tol);
%plot
figure(10)
s_uninhabitable = plot(T, fraction_uninhabitable, 'linewidth',2);
s_uninhabitable.Color='k';
shg
ax = gca;
ax.FontSize = 14;
xlabel('time, $T$','Interpreter','latex','fontsize',22)
ylabel('$P_U$','interpreter','latex','fontsize',22)
title({'fraction of domain which is uninhabitable', ['N_0= ' num2str(N0, '%1.1e') ' cell m^{-1}, [u]=' num2str(flow_scale, '%1.e') 'ms^{-1}']},'FontWeight','Normal')

xlim([0 7])
ylim([0 1])
set(gcf, 'Position',  [100, 100, 600, 500])

%% compute turn off time
[~,T_off] = turn_off_time_LW(T, W_LW, W_tol);

%% plot glucose conversion
[conversion_glucose, ~] = conversion_varyN0flow(T, C_LW, proliferation_rate, uptake_rate, N0, flow_scale);
s_conversion = plot(T, conversion_glucose, 'linewidth',2);
s_conversion.Color='k';
shg
ax = gca;
ax.FontSize = 14;
xlabel('time, $T$','Interpreter','latex','fontsize',22)
ylabel('glucose conversion, $Q$','interpreter','latex','fontsize',22)
title({'glucose conversion', ['N_0= ' num2str(N0, '%1.1e') ' cell m^{-1}, [u]=' num2str(flow_scale, '%1.e') 'ms^{-1}']},'FontWeight','Normal')
xlim([0 7])
ylim([0 1])

hold off
set(gcf, 'Position',  [100, 100, 600, 500])



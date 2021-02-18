% Plotting metabolite concentrations for sublimit model
%%
%proliferation rate 'low' sets P=1/6, 'mid' sets P=1/3, 'high' sets P=1
proliferation_rate='mid';
%uptake rate 'low' sets rho=0.027, 'mid' sets rho=0.27, 'high' sets rho=2.7
uptake_rate='mid';
%% compute C and W for region X \in [0,1] and T \in [0,7]
[T,X, C_SL,W_SL]=metabolite_concentration_SL(proliferation_rate, uptake_rate);

%% To compute C and W for coordinates (T,X)
T_coord=0;
X_coord=0;
[C, W] = metabolite_concentration_SL_tX(T_coord, X_coord, proliferation_rate, uptake_rate);

%%
figure(12)
surf(X, T, C_SL, 'edgecolor','none')
xlim([0 1])
ylim([0 7])
view([0 90])
colormap viridis
colorbar
caxis([0 1])
title({'{\bf\fontsize{14} Glucose concentration - sublimit model}'},'FontWeight','Normal')
ax = gca;
ax.FontSize = 14;
xlabel('$X$','Interpreter','latex')
ylabel('$T$','interpreter','latex')
set(gcf, 'Position',  [100, 100, 600, 500])

%%
figure(13)
surf(X, T, W_SL, 'edgecolor','none')
colormap viridis
xlim([0 1])
ylim([0 7])
view([0 90])
colormap viridis
colorbar
max_W_in_range= W_SL(1,801);
caxis([0 max_W_in_range])
title({'{\bf\fontsize{14} Lactate concentration - sublimit model}'},'FontWeight','Normal')
ax = gca;
ax.FontSize = 14;
xlabel('$X$','Interpreter','latex'), 
ylabel('$T$','interpreter','latex')
set(gcf, 'Position',  [100, 100, 600, 500])


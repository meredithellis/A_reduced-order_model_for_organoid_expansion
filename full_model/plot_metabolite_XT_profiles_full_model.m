%% COMSOL - Z-averaged full solution

load('av_CW_concentration_full_model')

% Create grid for plotting
time=t_full;
xp=(0:0.01:1)';
[xx,TT]=meshgrid(xp,time);

%% Plot C(X,T)
figure(103)
s3=surf(xx,TT,avC_data_full);
s3.EdgeColor = 'none';
%ylim([0 168])
view([0 90])
%yticks([0 1*24 2*24 3*24 4*24 5*24 6*24 7*24])
%yticklabels({'0','1','2','3','4','5','6','7'})
colormap viridis
colorbar
shg
caxis([0 1])
title({'{\bf\fontsize{14} Glucose concentration}';...
    '{\fontsize{12} Z-averaged full solution}'},'FontWeight','Normal')
ax = gca;
ax.FontSize = 14;
xlabel('horizontal position, $X$','Interpreter','latex','fontsize',22), 
ylabel('time, $T$','interpreter','latex','fontsize',22)
set(gcf, 'Position',  [100, 100, 600, 500])


%% Plot W(X,T)
figure(104)
s4=surf(xx,TT,avW_data_full);
s4.EdgeColor = 'none';
%ylim([0 168])
view([0 90])
%yticks([0 1*24 2*24 3*24 4*24 5*24 6*24 7*24])
%yticklabels({'0','1','2','3','4','5','6','7'})
colormap viridis
colorbar
caxis([0 1.75]) %1.7
title({'{\bf\fontsize{14} Lactate concentration}';...
    '{\fontsize{12} Z-averaged full solution }'},'FontWeight','Normal')
ax = gca;
ax.FontSize = 14;
xlabel('horizontal position, $X$','Interpreter','latex','fontsize',22),
ylabel('time, $T$','interpreter','latex','fontsize',22)
set(gcf, 'Position',  [100, 100, 600, 500])
%%

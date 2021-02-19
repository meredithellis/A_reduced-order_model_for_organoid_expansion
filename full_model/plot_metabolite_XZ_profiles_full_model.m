

% Load data for the full solution at times T=1,3,7

load('CWdata_XZ137_full_model')

    data_C=Cdata_XZ137;
    data_W=Wdata_XZ137;
    x = data_C(:,1) ; y = data_C(:,2) ; 
    
    days_considered=[1,3,7];
   %% 
    for k=1:length(days_considered) %looping over the days
        C_val = data_C(:,2+k) ; W_val = data_W(:,2+k);
dt = delaunayTriangulation(x,y) ;
tri = dt.ConnectivityList ;
xi = dt.Points(:,1) ; 
yi = dt.Points(:,2) ; 
F_C = scatteredInterpolant(x,y,C_val);
F_W = scatteredInterpolant(x,y,W_val);
C_val_i = F_C(xi,yi) ;
W_val_i = F_W(xi,yi);
%%%%%
figure(k)
trisurf(tri,xi,yi,C_val_i) 
view(2)
caxis([0 1])
hold on 
yline(1/3,'r','linewidth',2);
hold off
shading interp
colormap viridis
colorbar
shg
ax = gca;
xlabel('horizontal position, $X$','Interpreter','latex', 'fontsize',24)
ylabel('vertical position, $Z$','interpreter','latex','fontsize',24)

% titles
title({sprintf('Glucose concentration at t=%d',days_considered(k))} ,'Fontsize',14); %, 'interpreter','latex');

set(gcf, 'Position',  [100, 100, 600, 500])

%%%%%%%%
figure((length(days_considered)+k))
trisurf(tri,xi,yi,W_val_i) 
view(2)
caxis([0 1.75])
hold on 
yline(1/3,'r','linewidth',2);
hold off
shading interp
colormap viridis
colorbar
shg
title({sprintf('Lactate concentration at t=%d',days_considered(k))} ,'Fontsize',14); %, 'interpreter','latex');

ax = gca;
ax.FontSize = 14;
xlabel('horizontal position, $X$','Interpreter','latex', 'fontsize',24)
ylabel('vertical position, $Z$','interpreter','latex','fontsize',24)
shg
set(gcf, 'Position',  [100, 100, 600, 500])

    end



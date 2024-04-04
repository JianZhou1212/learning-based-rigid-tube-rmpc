clc
clear
close all
load('Results_2.mat');
load('parameters.mat');
N_sam_ini = Results_2.N_sam_ini;
W_Hat_Vol = Results_2.W_Hat_Vol;
W_hat_opt = Results_2.W_hat_opt;

mean_vol = mean(W_Hat_Vol, 2)';
std_vol = std(W_Hat_Vol, 0, 2)';
figure(1)
maxY = max([mean_vol - std_vol; mean_vol + std_vol]); 
minY = min([mean_vol - std_vol; mean_vol + std_vol]); 
yFill = [maxY, fliplr(minY)];
xFill = [(N_sam_ini), fliplr((N_sam_ini))];  
fill(xFill, yFill, [192, 192, 192]/255, 'linestyle', 'None'); 
hold on
plot(N_sam_ini, mean_vol, 'k', 'linewidth', 3);
hold on
plot(N_sam_ini, ones(length(N_sam_ini), 1)*W_hat_opt.volume( ), 'color', [128, 128, 128]/255, 'linestyle', '-.', 'linewidth', 2);
LE = legend('One standard deviation', 'Mean', 'Interpreter','latex','Location','southeast', 'NumColumns',1);
set(LE,'Fontsize',11);
set(gca, 'XScale', 'log');
box on
grid off
ylim([0, 0.06]);
xlim([4, 16384]);
xlabel('$|\mathcal{I}_0^w|$ [-]', 'Interpreter','latex');
ylabel('$\hat{{\rm V}}_0^{\star}$ [-]', 'Interpreter','latex');
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',14);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[3 5 8 8]);
set(gcf, 'PaperSize', [16 7]);
exportgraphics(gcf,'Fig_Vol_Converge.pdf','ContentType','vector');



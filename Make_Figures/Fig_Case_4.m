clc
clear
close all
load('Results_4.mat');
load('parameters.mat');
W_Hat_Vol = Results_4.W_Hat_Vol;
W_hat_opt = Results_4.W_hat_opt;
K_N = Results_4.K_N;
Samples = Results_4.Samples;
steps = 11:K_N + 10;

W_Hat_Vol = W_Hat_Vol';
mean_vol = mean(W_Hat_Vol, 2)';
std_vol = std(W_Hat_Vol, 0, 2)';
figure(1)
maxY = max([mean_vol - std_vol; mean_vol + std_vol]); 
minY = min([mean_vol - std_vol; mean_vol + std_vol]); 
yFill = [maxY, fliplr(minY)];
xFill = [(steps), fliplr((steps))];  
fill(xFill, yFill, [192, 192, 192]/255, 'linestyle', 'None'); 
hold on
plot(steps, mean_vol, 'k', 'linewidth', 3);
hold on
plot(steps, ones(length(steps), 1)*W_hat_opt.volume( ), 'color', [128, 128, 128]/255, 'linestyle', '-.', 'linewidth', 2);
LE = legend('One standard deviation', 'Mean', 'Interpreter','latex','Location','southeast', 'NumColumns',1);
set(LE,'Fontsize',11);
set(gca, 'XScale', 'log');
box on
ylim([0, 0.06]);
xlim([11, 10000]);
xlabel('$|\mathcal{I}_k^w|$ [-]', 'Interpreter','latex');
ylabel('$\hat{{\rm V}}_k^{\star}$ [-]', 'Interpreter','latex');
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',14);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[3 5 8 8]);
set(gcf, 'PaperSize', [16 7]);
exportgraphics(gcf,'Fig_Vol_Converge_Online.pdf','ContentType','vector');



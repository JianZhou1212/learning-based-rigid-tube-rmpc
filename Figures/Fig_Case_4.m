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
%%
plot(Samples{3}(1, :),Samples{3}(2, :), 'k.', 'markersize', 8);
hold on
plot(parameters.opts_ini_set.W_true, 'wire', 1, 'edgecolor', 'm', 'linewidth', 2.5);
%%
W_Hat_Vol = W_Hat_Vol';
mean_vol = mean(W_Hat_Vol, 2)';
std_vol = std(W_Hat_Vol, 0, 2)';
figure(1)
maxY = max([mean_vol - std_vol; mean_vol + std_vol]); 
minY = min([mean_vol - std_vol; mean_vol + std_vol]); 
yFill = [maxY, fliplr(minY)];
xFill = [(steps), fliplr((steps))];  
fill(xFill, yFill, 'c', 'linestyle', 'None'); 
hold on
plot(steps, mean_vol, 'b', 'linewidth', 3);
hold on
plot(steps, ones(length(steps), 1)*W_hat_opt.volume( ), 'k-.', 'linewidth', 3);
LE = legend('One standard deviation', 'Mean', 'Interpreter','latex','Location','southeast', 'NumColumns',1);
set(LE,'Fontsize',12);
set(gca, 'XScale', 'log');
box on
ylim([0, 0.07]);
xlim([11, 10000]);
xlabel('$|\mathcal{I}_k^w|$ [-]', 'Interpreter','latex');
ylabel('${\rm Set \ Volume}$ [-]', 'Interpreter','latex');
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',12);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[3 5 17 5]);
set(gcf, 'PaperSize', [16 7]);
exportgraphics(gcf,'Fig_Vol_Converge_Online.pdf','ContentType','vector');



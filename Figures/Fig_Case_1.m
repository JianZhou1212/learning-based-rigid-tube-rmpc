% Load parameters
clc
clear
close all
load('Results_1.mat');
load('parameters.mat');
N_sam_ini = Results_1.N_sam_ini;
Samples_ini = Results_1.Samples_ini;
Alpha_ini = Results_1.Alpha_ini;
V_ini = Results_1.V_ini;
W_Hat_ini = Results_1.W_Hat_ini;
W = Results_1.W;
W_hat_opt = Results_1.W_hat_opt;
W_true = Results_1.W_true;
F_N_True = Results_1.F_N_True;
F_N_Hat_Opt = Results_1.F_N_Hat_Opt;
F_N_RMPC = Results_1.F_N_RMPC;

Phi = parameters.opts_RMPC.Phi;
S = parameters.opts_RMPC.S;
Nu = parameters.opts_RMPC.Nu;

FeasibleRegion = ComputeFeasibleRegion(parameters.opts_feasible_region);
%%
close all
for i = 1:length(N_sam_ini)
    figure(i)
    plot(W, 'wire', 1, 'edgecolor', 'b', 'linewidth', 2.5);
    hold on
    plot(W_Hat_ini{i}, 'wire', 1, 'edgecolor', 'c', 'linewidth', 2.5);
    hold on
    plot(W_true, 'wire', 1, 'edgecolor', 'm', 'linewidth', 2.5);
    hold on
    plot(W_hat_opt, 'wire', 1, 'edgecolor', [178, 102, 255]/255, 'linewidth', 2.5);
    hold on
    plot(Samples_ini{i}(1, :),Samples_ini{i}(2, :), 'k.','markersize', 8);
    box on
    grid off
    xlim([-0.62, 0.62]);
    ylim([-0.25, 0.25]);
    xlabel('$w_{k, 1}$ [${\rm m}$]', 'Interpreter','latex');
    ylabel('$w_{k, 2}$ [${\rm m/s}$]', 'Interpreter','latex');
    set(gca,'Linewidth',1.5,'GridAlpha',0.5);
    set(gca,'FontName','Times New Roman','FontSize',15);
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gcf,'unit','centimeters','position',[5 5 10 10]);
    set(gcf, 'PaperSize', [16 7]);
    savename = sprintf('Fig_W_hat_vs_W_Samples_%d.pdf', N_sam_ini(i));
    exportgraphics(gcf, savename,'ContentType','vector');
end
%%
close all
for i = 1:length(N_sam_ini)
    figure(i)
    S_hat_ini = Alpha_ini(i)*S + (1 - Alpha_ini(i))*(inv(1 - Phi)*V_ini(:, i));
    hs_UQMPC = FeasibleRegion.ComFeasibleRegion_UQMPC_Hs_ini(S_hat_ini);
    [~, F_N_Hat] = FeasibleRegion.ComFeasibleRegion_UQMPC(S_hat_ini, hs_UQMPC);
    plot(F_N_True, 'wire', 1, 'edgecolor', 'm', 'linewidth', 2.5);
    hold on
    plot(F_N_Hat, 'wire', 1, 'edgecolor', 'c', 'linewidth', 2.5);
    hold on
    plot(F_N_RMPC, 'wire', 1, 'edgecolor', 'b', 'linewidth', 2.5);
    hold on
    plot(F_N_Hat_Opt, 'wire', 1, 'edgecolor', [178, 102, 255]/255, 'linewidth', 2.5);
    box on
    grid off
    xlim([-16, 16]);
    ylim([-8, 8]);
    xlabel('$x_{0,1}$ [${\rm m}$]', 'Interpreter','latex');
    ylabel('$x_{0, 2}$ [${\rm m/s}$]', 'Interpreter','latex');
    set(gca,'Linewidth',1.5,'GridAlpha',0.5);
    set(gca,'FontName','Times New Roman','FontSize',15);
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gcf,'unit','centimeters','position',[5 5 10 10]);
    set(gcf, 'PaperSize', [16 7]);
    savename = sprintf('Fig_F_N_hat_vs_F_N_%d.pdf', N_sam_ini(i));
    exportgraphics(gcf, savename,'ContentType','vector');
end
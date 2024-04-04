% load results
clc
clear
close all
load('parameters.mat');
load('Results_3_small.mat');
load('Results_3_large.mat');
t_small = Results_3_small.t_small;
Control_EV_UQMPC_small = Results_3_small.Control_EV_UQMPC_small;
State_RM_UQMPC_small = Results_3_small.State_RM_UQMPC_small;
Nominal_RM_State_UQMPC_small = Results_3_small.Nominal_RM_State_UQMPC_small;
S_Hat_small = Results_3_small.S_Hat_small;

t_large = Results_3_large.t_large;
Control_EV_UQMPC_large = Results_3_large.Control_EV_UQMPC_large;
State_RM_UQMPC_large = Results_3_large.State_RM_UQMPC_large;
Nominal_RM_State_UQMPC_large = Results_3_large.Nominal_RM_State_UQMPC_large;
S_Hat_large = Results_3_large.S_Hat_large;

acc_bound = parameters.acc_bound;
c_set = [192, 192, 192]/255;
%%
clc
close all
figure(1)
num_steps = 10;
for i = 1:1:num_steps
    Nominal_State = [Nominal_RM_State_UQMPC_small(1, i); Nominal_RM_State_UQMPC_small(2, i)];
    Graphics.show_convex(S_Hat_small{i}  + Nominal_State, c_set, 'FaceAlpha', 0.2);
end
h1 = plot(State_RM_UQMPC_small(1, 1:num_steps), State_RM_UQMPC_small(2, 1:num_steps), 'color', [48, 129, 208]/255, 'linestyle', '-', 'marker', '.', 'linewidth', 1.5, 'markersize', 8);
hold on
h2 = plot(Nominal_RM_State_UQMPC_small(1, 1:num_steps), Nominal_RM_State_UQMPC_small(2, 1:num_steps), 'color', [238, 114, 20]/255, 'linestyle', '-.', 'marker', '.', 'linewidth', 1.5, 'markersize', 8);
box on
grid off
xlim([-13, 1]);
LE = legend([h1, h2], '${\rm True \ State}(|\mathcal{I}_0^w = 100|)$', '${\rm Nom. \ State}(|\mathcal{I}_0^w = 100|)$',  'Interpreter','latex', 'Location','southwest');
set(LE,'Fontsize',8);
xlabel('$x_{k,1}$ [${\rm m}$]', 'Interpreter','latex');
ylabel('$x_{k,2}$ [${\rm m/s}$]', 'Interpreter','latex');
set(LE, 'Fontsize', 12);
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',14);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[5 5 8 8]);
set(gcf, 'PaperSize', [16 7]);
savename = sprintf('Fig_Trajectory_Small_Initial_Sample.pdf');
exportgraphics(gcf, savename,'ContentType','vector');
%%
figure(2)
num_steps = 10;
for i = 1:1:num_steps
    Nominal_State = [Nominal_RM_State_UQMPC_large(1, i); Nominal_RM_State_UQMPC_large(2, i)];
    Graphics.show_convex(S_Hat_large{i}  + Nominal_State,  c_set, 'FaceAlpha', 0.2);
end
hold on
h1 = plot(State_RM_UQMPC_large(1, 1:num_steps), State_RM_UQMPC_large(2, 1:num_steps),  'color', [48, 129, 208]/255, 'linestyle', '-', 'marker', '.', 'linewidth', 1.5, 'markersize', 8);
hold on
h2 = plot(Nominal_RM_State_UQMPC_large(1, 1:num_steps), Nominal_RM_State_UQMPC_large(2, 1:num_steps), 'color', [238, 114, 20]/255, 'linestyle', '-.', 'marker', '.', 'linewidth', 1.5, 'markersize', 8);
box on
grid off
xlim([-13, 1]);
LE = legend([h1, h2], '${\rm True \ State}(|\mathcal{I}_0^w = 2{\rm e}4|)$', '${\rm Nom. \ State}(|\mathcal{I}_0^w = 2{\rm e}4|)$',  'Interpreter','latex', 'Location','southwest');
set(LE,'Fontsize',8);
xlabel('$x_{k,1}$ [${\rm m}$]', 'Interpreter','latex');
ylabel('$x_{k,2}$ [${\rm m/s}$]', 'Interpreter','latex');
set(LE, 'Fontsize', 12);
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',14);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[5 5 8 8]);
set(gcf, 'PaperSize', [16 7]);
savename = sprintf('Fig_Trajectory_Large_Initial_Sample.pdf');
exportgraphics(gcf, savename,'ContentType','vector');

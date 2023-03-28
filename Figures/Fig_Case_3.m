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
%%
figure(1)
h1 = plot(t_small(1:30), Control_EV_UQMPC_small(1:30), 'm', 'linewidth', 3);
hold on
h2 = plot(t_large(1:30), Control_EV_UQMPC_large(1:30), 'b-.', 'linewidth', 3);
hold on
h3 = plot(t_small(1:30),-acc_bound*ones(1, length(Control_EV_UQMPC_small(1:30))), 'k--', 'linewidth', 1.5);
hold on
plot(t_small(1:30),acc_bound*ones(1, length(Control_EV_UQMPC_small(1:30))), 'k--', 'linewidth', 1.5);
ylim([-2.5, 2.5])
box on
grid off
LE = legend([h1, h2, h3], '$|\mathcal{I}_0^w| = 100$','$|\mathcal{I}_0^w| = 20000$', '${\rm Acc. bound}$', 'Interpreter','latex');
xlabel('${\rm Time} $ [${\rm s}$]', 'Interpreter','latex');
ylabel('${\rm Long. \ Acc.}$ [${\rm m/s^2}$]', 'Interpreter','latex');
set(LE, 'Fontsize', 10);
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',12);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[15 5 15 5]);
set(gcf, 'PaperSize', [16 7]);
savename = sprintf('Fig_control_large_small_sample.pdf');
exportgraphics(gcf, savename,'ContentType','vector');
%%
figure(2)
num_steps = 10;
for i = 1:1:num_steps
    Nominal_State = [Nominal_RM_State_UQMPC_small(1, i); Nominal_RM_State_UQMPC_small(2, i)];
    Graphics.show_convex(S_Hat_small{i}  + Nominal_State, 'g', 'FaceAlpha', 0.2);
end
h1 = plot(State_RM_UQMPC_small(1, 1:num_steps), State_RM_UQMPC_small(2, 1:num_steps), 'kd--', 'linewidth', 1.5, 'MarkerFaceColor','k', 'markersize', 4);
hold on
h2 = plot(Nominal_RM_State_UQMPC_small(1, 1:num_steps), Nominal_RM_State_UQMPC_small(2, 1:num_steps), 'r.--','linewidth',1, 'markersize', 18);
box on
grid off
xlim([-13, 1]);
LE = legend([h1, h2], '${\rm True \ State}(|\mathcal{I}_0^w = 100|)$', '${\rm Nominal \ State}(|\mathcal{I}_0^w = 100|)$',  'Interpreter','latex', 'Location','southwest');
xlabel('$x_{k,1}$ [${\rm m}$]', 'Interpreter','latex');
ylabel('$x_{k,2}$ [${\rm m/s}$]', 'Interpreter','latex');
set(LE, 'Fontsize', 12);
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',14);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[5 5 10 10]);
set(gcf, 'PaperSize', [16 7]);
savename = sprintf('Fig_Trajectory_Small_Initial_Sample.pdf');
exportgraphics(gcf, savename,'ContentType','vector');
%%
figure(3)
num_steps = 10;
for i = 1:1:num_steps
    Nominal_State = [Nominal_RM_State_UQMPC_large(1, i); Nominal_RM_State_UQMPC_large(2, i)];
    Graphics.show_convex(S_Hat_large{i}  + Nominal_State, 'g', 'FaceAlpha', 0.2);
end
hold on
h1 = plot(State_RM_UQMPC_large(1, 1:num_steps), State_RM_UQMPC_large(2, 1:num_steps),  'kd--', 'linewidth', 1.5, 'MarkerFaceColor','k', 'markersize', 4);
hold on
h2 = plot(Nominal_RM_State_UQMPC_large(1, 1:num_steps), Nominal_RM_State_UQMPC_large(2, 1:num_steps), 'r.--','linewidth',1, 'markersize', 18);
box on
grid off
xlim([-13, 1]);
LE = legend([h1, h2], '${\rm True \ State}(|\mathcal{I}_0^w = 20000|)$', '${\rm Nominal \ State}(|\mathcal{I}_0^w = 20000|)$',  'Interpreter','latex', 'Location','southwest');
xlabel('$x_{k,1}$ [${\rm m}$]', 'Interpreter','latex');
ylabel('$x_{k,2}$ [${\rm m/s}$]', 'Interpreter','latex');
set(LE, 'Fontsize', 12);
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',14);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[5 5 10 10]);
set(gcf, 'PaperSize', [16 7]);
savename = sprintf('Fig_Trajectory_Large_Initial_Sample.pdf');
exportgraphics(gcf, savename,'ContentType','vector');

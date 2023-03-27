% In this scenario, we try to compare the feasible region of RMPC and
% UQ-RMPC
clc
clear
close all
load('parameters.mat');
opts_ini_set = parameters.opts_ini_set;
opts_Car = parameters.opts_Car;
opts_UQMPC = parameters.opts_UQMPC;
x_des = parameters.x_des;
distance_error = parameters.distance_error;
acc_bound = parameters.acc_bound;
nx = opts_UQMPC.nx;
nu = opts_UQMPC.nu;
nc = opts_UQMPC.nc;
A = opts_UQMPC.A;
B = opts_UQMPC.B;
T = parameters.T;
S = parameters.opts_feasible_region.S;
Phi = parameters.opts_UQMPC.Phi;
%
import casadi.*
% call UQMPC and Car modeling
UQRobustMPC = UQRobustMPC(opts_UQMPC);
Car = ModelingCar(opts_Car);

opts_ini_set.N_pre_sam = 5000; % change this to 10, 100, 500, 2000, 5000, corresponding to |I_0^w|
IniSet = InitialSetComputation(opts_ini_set);

% this checks if the state is close to the bound of initial feasible region
[alpha_ini, v_ini, ~] = IniSet.solve();
S_hat_ini = alpha_ini*S + (1 - alpha_ini)*(inv(1 - Phi)*v_ini);
FeasibleRegion = ComputeFeasibleRegion(parameters.opts_feasible_region);
hs_UQMPC = FeasibleRegion.ComFeasibleRegion_UQMPC_Hs_ini(S_hat_ini);
[~, F_N_Hat] = FeasibleRegion.ComFeasibleRegion_UQMPC(S_hat_ini, hs_UQMPC, parameters.opts_RMPC.Nu);
V = F_N_Hat.V;

% Online iteraction
K_N = 20; % the MPC iterates 20 steps online
x_LV_0 = [100; 10]; % initial state of LV
x_RM_0 = [min(V(:, 1)) + 0.1; max(V(:, 2)) - 0.3]; % initial relative state, this is selected at the bound of feasible region
%%
x_EV_0 = x_LV_0 + x_des + x_RM_0;
N_MC = 300; % number of Monte-Carlo simulations of each |I_0^w|
Control_EV = cell(N_MC, 1);
State_RM = cell(N_MC, 1);
Infeasible_Index = zeros(N_MC, 1);
for i = 1:1:N_MC
    % change the number of samples for initial quantification
    [alpha_ini, v_ini, ~] = IniSet.solve();

    State_LV = ones(nx, K_N + 1);
    State_LV(:, 1) = x_LV_0;

    Control_EV_UQMPC = ones(nu, K_N);
    State_EV_UQMPC = ones(nx, K_N + 1);
    State_RM_UQMPC = ones(nx, K_N + 1);
    State_EV_UQMPC(:, 1) = x_EV_0;
    State_RM_UQMPC(:, 1) = x_EV_0 - x_LV_0 - x_des;
    W_Hat_Vol_loop = ones(1, K_N);
    for k = 1:K_N
        [xi_LV_k, x_LV_k_next] = Car.LVModeling(State_LV(:, k));
        State_LV(:, k + 1) = x_LV_k_next;

        if k == 1
            w_new = [0; 0];
            alpha_before = alpha_ini;
            v_before = v_ini;
        else
            w_new = State_RM_UQMPC(:, k) - A*State_RM_UQMPC(:, k-1) - B*Control_EV_UQMPC(:, k-1);
            alpha_before = alpha_k;
            v_before = v_k;
        end
        [s_k_uqmpc, u_EV_k_uqmpc, alpha_k, v_k, hs_uqmpc, S_hat, W_hat] = UQRobustMPC.solve(State_RM_UQMPC(:, k), alpha_before, v_before, w_new);
        [xi_EV_k_uqmpc, x_EV_k_next_uqmpc] = Car.EVModeling(State_EV_UQMPC(:, k), u_EV_k_uqmpc);
        Control_EV_UQMPC(:, k) = u_EV_k_uqmpc;
        State_EV_UQMPC(:, k + 1) = x_EV_k_next_uqmpc;
        State_RM_UQMPC(:, k + 1) = x_EV_k_next_uqmpc - x_LV_k_next - x_des;
    end
    Control_EV{i} = Control_EV_UQMPC;
    State_RM{i} = State_RM_UQMPC;
    X_Error = State_RM_UQMPC(1, :);
    exceed_up_control = Control_EV_UQMPC(Control_EV_UQMPC > acc_bound);
    exceed_low_control = Control_EV_UQMPC(Control_EV_UQMPC < -acc_bound);
    exceed_up_x = X_Error(X_Error > distance_error);
    exceed_low_x = X_Error(X_Error < -distance_error);
    index = length(exceed_up_control) + length(exceed_low_control) + length(exceed_up_x) + length(exceed_low_x);
    if index~= 0
       Infeasible_Index(i) = 1;
    end
end
sum(Infeasible_Index)/N_MC % success rate
%
if opts_ini_set.N_pre_sam == 10
    Results_5_10.x_RM_0 = x_RM_0;
    Results_5_10.N_pre_sam = opts_ini_set.N_pre_sam;
    Results_5_10.Control_EV =  Control_EV;
    Results_5_10.State_RM = State_RM;
    Results_5_10.Infeasible_Index = Infeasible_Index;
    save('Results_5_10.mat', 'Results_5_10');
end

if opts_ini_set.N_pre_sam == 100
    Results_5_100.x_RM_0 = x_RM_0;
    Results_5_100.N_pre_sam = opts_ini_set.N_pre_sam;
    Results_5_100.Control_EV =  Control_EV;
    Results_5_100.State_RM = State_RM;
    Results_5_100.Infeasible_Index = Infeasible_Index;
    save('Results_5_100.mat', 'Results_5_100');
end

if opts_ini_set.N_pre_sam == 500
    Results_5_500.x_RM_0 = x_RM_0;
    Results_5_500.N_pre_sam = opts_ini_set.N_pre_sam;
    Results_5_500.Control_EV =  Control_EV;
    Results_5_500.State_RM = State_RM;
    Results_5_500.Infeasible_Index = Infeasible_Index;
    save('Results_5_500.mat', 'Results_5_500');
end

if opts_ini_set.N_pre_sam == 2000
    Results_5_2000.x_RM_0 = x_RM_0;
    Results_5_2000.N_pre_sam = opts_ini_set.N_pre_sam;
    Results_5_2000.Control_EV =  Control_EV;
    Results_5_2000.State_RM = State_RM;
    Results_5_2000.Infeasible_Index = Infeasible_Index;
    save('Results_5_2000.mat', 'Results_5_2000');
end

if opts_ini_set.N_pre_sam == 5000
    Results_5_5000.x_RM_0 = x_RM_0;
    Results_5_5000.N_pre_sam = opts_ini_set.N_pre_sam;
    Results_5_5000.Control_EV =  Control_EV;
    Results_5_5000.State_RM = State_RM;
    Results_5_5000.Infeasible_Index = Infeasible_Index;
    save('Results_5_5000.mat', 'Results_5_5000');
end




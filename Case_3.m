% In this scenario, we try to compare the feasible region of RMPC and
% UQ-RMPC
clc
clear
close all
load('parameters.mat');
opts_ini_set = parameters.opts_ini_set;
opts_Car = parameters.opts_Car;
opts_UQMPC = parameters.opts_UQMPC;
F_N_True = parameters.F_N_True;
x_des = parameters.x_des;
nx = opts_UQMPC.nx;
nu = opts_UQMPC.nu;
nc = opts_UQMPC.nc;
A = opts_UQMPC.A;
B = opts_UQMPC.B;
T = parameters.T;
%%
import casadi.*
% call UQMPC and Car modeling
UQRobustMPC = UQRobustMPC(opts_UQMPC); % call the function of UQ-Robust MPC
Car = ModelingCar(opts_Car); % call the function for modeling the car

% change the number of samples for initial quantification
%opts_ini_set.N_pre_sam = 100; % |I_0^w| = 100
opts_ini_set.N_pre_sam = 20000; % |I_0^w| = 20000
IniSet = InitialSetComputation(opts_ini_set);
[alpha_ini, v_ini, samples] = IniSet.solve();

% Online iteraction
K_N = 50; % the MPC iterates 50 steps online
x_LV_0 = [100; 10]; % initial state of LV
x_RM_0 = [-12; 5];  % initial relative state
x_EV_0 = x_LV_0 + x_des + x_RM_0; % initial state of EV

State_LV = ones(nx, K_N + 1); % state of the LV
State_LV(:, 1) = x_LV_0;   
Input_LV = ones(nx, K_N);     % input of the LV, the input is control + disturbance

Control_EV_UQMPC = ones(nu, K_N); % control input of EV
Noise_EV_UQMPC = ones(nx, K_N);   % disturbance of the EV
State_EV_UQMPC = ones(nx, K_N + 1); % state of the EV
State_RM_UQMPC = ones(nx, K_N + 1); % relative state
Nominal_RM_State_UQMPC = ones(nx, K_N + 1); % nominal relative state
Hs_UQMPC       = ones(nc, K_N); % value of hs vector
State_EV_UQMPC(:, 1) = x_EV_0;
State_RM_UQMPC(:, 1) = x_EV_0 - x_LV_0 - x_des;
Nominal_RM_State_UQMPC(:, 1) = x_RM_0;

Alpha = ones(1, K_N + 1);
V = ones(nx, K_N + 1);
S_Hat = cell(K_N, 1);
W_Hat = cell(K_N, 1);
Alpha(:, 1) = alpha_ini;
V(:, 1) = v_ini;

Time_UQMPC = ones(K_N, 1);

tic
for k = 1:K_N
    [xi_LV_k, x_LV_k_next] = Car.LVModeling(State_LV(:, k));
    State_LV(:, k + 1) = x_LV_k_next;
    Input_LV(:, k) = xi_LV_k;

    if k == 1
        w_new = [0; 0];
        alpha_before = alpha_ini;
        v_before = v_ini;
    else
        w_new = State_RM_UQMPC(:, k) - A*State_RM_UQMPC(:, k-1) - B*Control_EV_UQMPC(:, k-1);
        alpha_before = alpha_k;
        v_before = v_k;
    end
    tic
    [s_k_uqmpc, u_EV_k_uqmpc, alpha_k, v_k, hs_uqmpc, S_hat, W_hat] = UQRobustMPC.solve(State_RM_UQMPC(:, k), alpha_before, v_before, w_new);
    toc
    Time_UQMPC(k) = toc;
    [xi_EV_k_uqmpc, x_EV_k_next_uqmpc] = Car.EVModeling(State_EV_UQMPC(:, k), u_EV_k_uqmpc);
    Control_EV_UQMPC(:, k) = u_EV_k_uqmpc;
    Noise_EV_UQMPC(:, k)   = xi_EV_k_uqmpc;
    State_EV_UQMPC(:, k + 1) = x_EV_k_next_uqmpc;
    State_RM_UQMPC(:, k + 1) = x_EV_k_next_uqmpc - x_LV_k_next - x_des;
    Nominal_RM_State_UQMPC(:, k) = s_k_uqmpc;
    Hs_UQMPC(:, k) = hs_uqmpc;
    S_Hat{k} = S_hat;
    W_Hat{k} = W_hat;
    
    Alpha(:, k + 1) = alpha_k;
    V(:, k + 1) = v_k;
    samples = [samples w_new];
end
toc
fprintf('Every Excecution Time is %d.\n', toc/K_N);
t = 0:T:T*K_N;
Time_UQMPC = Time_UQMPC(2:end);
%%
if opts_ini_set.N_pre_sam == 100
    Results_3_small.Control_EV_UQMPC_small = Control_EV_UQMPC;
    Results_3_small.State_RM_UQMPC_small = State_RM_UQMPC;
    Results_3_small.Nominal_RM_State_UQMPC_small = Nominal_RM_State_UQMPC;
    Results_3_small.Hs_UQMPC_small       = Hs_UQMPC;
    Results_3_small.S_Hat_small = S_Hat;
    Results_3_small.W_Hat_Small = W_Hat;
    Results_3_small.t_small = t;
    Results_3_small.Time_UQMPC_small = Time_UQMPC;
    save('Results_3_small.mat', 'Results_3_small');
end

if opts_ini_set.N_pre_sam == 20000
    Results_3_large.Control_EV_UQMPC_large = Control_EV_UQMPC;
    Results_3_large.State_RM_UQMPC_large = State_RM_UQMPC;
    Results_3_large.Nominal_RM_State_UQMPC_large = Nominal_RM_State_UQMPC;
    Results_3_large.Hs_UQMPC_large = Hs_UQMPC;
    Results_3_large.S_Hat_large = S_Hat;
    Results_3_large.W_Hat_large = W_Hat;
    Results_3_large.t_large = t;
    Results_3_large.Time_UQMPC_large = Time_UQMPC;
    save('Results_3_large.mat', 'Results_3_large');
end








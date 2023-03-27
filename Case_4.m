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
nx = opts_UQMPC.nx;
nu = opts_UQMPC.nu;
nc = opts_UQMPC.nc;
A = opts_UQMPC.A;
B = opts_UQMPC.B;
T = parameters.T;
%
import casadi.*
% call UQMPC and Car modeling
UQRobustMPC = UQRobustMPC(opts_UQMPC);
Car = ModelingCar(opts_Car);

opts_ini_set.N_pre_sam = 10; % value of |I_0^w|
IniSet = InitialSetComputation(opts_ini_set);

% Online iteraction
K_N = 10000; % the MPC iterates 10000 steps on line
x_LV_0 = [100; 10]; % initial state of LV
x_RM_0 = [-12; 5];  % relative initial state
x_EV_0 = x_LV_0 + x_des + x_RM_0;
N_MC = 30;
W_Hat_Vol = ones(N_MC, K_N);
Samples = cell(N_MC, 1);
for i = 1:1:N_MC
    % change the number of samples for initial quantification
    [alpha_ini, v_ini, samples] = IniSet.solve();

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
        W_Hat_Vol_loop(k) = W_hat.volume( );
        samples = [samples w_new];
    end
    W_Hat_Vol(i, :) = W_Hat_Vol_loop;
    Samples{i} = samples;
    
    figure(i)
    plot(samples(1, :),samples(2, :), 'k.', 'markersize', 8);
    hold on
    plot(parameters.W_true, 'wire', 1, 'edgecolor', 'm', 'linewidth', 2.5);
    
end
%
Results_4.t = 0:T:T*K_N;
Results_4.K_N = K_N;
Results_4.W_Hat_Vol =  W_Hat_Vol;
Results_4.W_hat_opt = parameters.W_hat_opt;
Results_4.Samples = Samples;
save('Results_4.mat', 'Results_4');







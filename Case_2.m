clc
clear
close all
load('parameters.mat');
opts_ini_set = parameters.opts_ini_set;
W = opts_ini_set.W;
%%
N_sam_ini = [ ];
for i = 2:0.5:14
    N_sam_ini = [N_sam_ini round(2^i)]; % change |I_0^w|
end
N_MC = 30; % number of Monte-Carlo simulations for each |I_0^w|
W_Hat = cell(length(N_sam_ini), 1);
W_Hat_Vol = ones(length(N_sam_ini), N_MC);
for i = 1:1:length(N_sam_ini)
    W_hat_mc = cell(N_MC, 1);
    W_hat_mc_vol = ones(N_MC, 1);
    for j = 1:1:N_MC
        opts_ini_set.N_pre_sam = N_sam_ini(i);
        IniSet = InitialSetComputation(opts_ini_set);
        [alpha_ini, v_ini, ~] = IniSet.solve();
        W_hat_j = (1 - alpha_ini)*v_ini + alpha_ini*W;
        W_hat_mc{j} = W_hat_j;
        W_hat_mc_vol(j) = W_hat_j.volume( );
    end
    W_Hat_Vol(i, :) = W_hat_mc_vol';
    W_Hat{i} = W_hat_mc;
end

Results_2.N_sam_ini = N_sam_ini;
Results_2.W_Hat_Vol = W_Hat_Vol;
Results_2.W_Hat = W_Hat;
Results_2.W = W;
Results_2.W_hat_opt = parameters.W_hat_opt;
Results_2.W_true = parameters.W_true;
Results_2.F_N_True = parameters.F_N_True;
Results_2.F_N_Hat_Opt = parameters.F_N_Hat_Opt;
Results_2.F_N_RMPC = parameters.F_N_RMPC;
save('Results_2.mat', 'Results_2')

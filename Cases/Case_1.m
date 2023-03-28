clc
clear
close all
load('parameters.mat');
opts_ini_set = parameters.opts_ini_set;
W = opts_ini_set.W;
%%
N_sam_ini = [5; 50; 500; 2000; 20000]; % change |I_0^w|
Alpha_ini = ones(length(N_sam_ini), 1);
V_ini = ones(2, length(N_sam_ini));
Samples_ini = cell(length(N_sam_ini), 1);
W_Hat_ini = cell(length(N_sam_ini), 1);
for k = 1:1:length(N_sam_ini)
    opts_ini_set.N_pre_sam = N_sam_ini(k);
    IniSet = InitialSetComputation(opts_ini_set); 
    [alpha_ini, v_ini, samples] = IniSet.solve(); % return \alpha_0, v_0, and samples
    Alpha_ini(k) = alpha_ini;
    V_ini(:, k) = v_ini;
    Samples_ini{k} = samples;
    W_Hat_ini{k} = (1 - alpha_ini)*v_ini + alpha_ini*W;
end

Results_1.N_sam_ini = N_sam_ini;
Results_1.Alpha_ini = Alpha_ini;
Results_1.V_ini = V_ini;
Results_1.Samples_ini = Samples_ini;
Results_1.W_Hat_ini = W_Hat_ini;
Results_1.W = W;
Results_1.W_hat_opt = parameters.W_hat_opt;
Results_1.W_true = parameters.W_true;
Results_1.F_N_True = parameters.F_N_True;
Results_1.F_N_Hat_Opt = parameters.F_N_Hat_Opt;
Results_1.F_N_RMPC = parameters.F_N_RMPC;
save('Results_1.mat', 'Results_1')

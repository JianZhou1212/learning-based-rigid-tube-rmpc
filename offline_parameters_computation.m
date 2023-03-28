% General Parameters
import casadi.*
clc
clear
close all
T = 0.5; % sampling interval
N = 8;   % prediction horizon
A = [1 T; 0 1]; % A matrix
B = [0; T];     % B matrix
Q = [1 0; 0 1]; % Q matrix
R = 0.1;        % R matrix
[K, Px] = dlqr(A, B, Q, R); % solving the discrete-time algebraic Riccati function
K = -K;                     
Phi = A + B*K;              % \Phi matrix, strictly stable
theta_EV = pi/40;           % define the rotation angle of \Xi_true_EV (rotation just makes the set more complex)
theta_LV = pi/40;           % define the rotation angle of \Xi_true_LV
Rotation_EV = [cos(theta_EV) -sin(theta_EV); sin(theta_EV) cos(theta_EV)]; % rotation matrix of \Xi_true_EV
Rotation_LV = [cos(theta_LV) -sin(theta_LV); sin(theta_LV) cos(theta_LV)]; % rotation matrix of \Xi_true_LV
Xi_true_EV = Rotation_EV*Polyhedron([-0.06 -0.015;0.06 -0.015; 0.01 0.025; -0.01 0.025]); % \Xi_true_EV
Xi_true_LV = Rotation_LV*Polyhedron([-0.06 0.015;0.06 0.015; 0.01 -0.025; -0.01 -0.025]); % \Xi_true_LV
min_u_LV = -1/20;  % lower bound of acceleration of LV
max_u_LV = 1/16;   % upper bound of acceleration of LV
U_true_LV  = Polyhedron([1/max_u_LV; 1/min_u_LV], [1; 1]); % U_true_LV
W_true = Xi_true_EV + (-1*Xi_true_LV) + (-B*U_true_LV);    % define the true disturbance set of the relative model, W_true
W = Polyhedron([-0.5 -0.2; 0.5 -0.2; 0.5 0.2; -0.5 0.2]);  % the admissible disturbance set of the relative model, a superset of W_true, overly estimated
Xi_true_EV = minHRep(Xi_true_EV); % reduce redundant vertices of the set
Xi_true_LV = minHRep(Xi_true_LV); % reduce redundant vertices of the set
W_true = minHRep(W_true);         % reduce redundant vertices of the set
W = minHRep(W);                   % reduce redundant vertices of the set
samples_opt = W_true.V;           % the optimal samples, i.e., the vertices of W_true, to quantify W_true based on W
[N_sam_opt, ~] = size(samples_opt);

% impose mix constraints
distance_error = 15; % constraint on the relative distance between EV and LV
acc_bound = 2;       % constraint on acceleration of EV
F = [1/distance_error 0; -1/distance_error 0; 0 0;0 0]; % F matrix of the mixed constraint
G = [0; 0; 1/acc_bound; -1/acc_bound];  % G matrix of the mixed constrain
x_des = [-35; 0]; % desired state: relative long. distance and relative long. speed
nx = 2; % number of state
nu = 1; % number of input
nc = 4; % column of G matrix
% following to define the matrices in the OPT problems
E = zeros(1, N);
E(1) = 1; % the E matrix
diag_ele = ones(N, 1);
for i = 1:1:N
    diag_ele(i) = B'*Px*B + R;
end
Pc = diag(diag_ele); % the Pc matrix
M = diag(ones(N-1, 1), 1); % the M matrix in Psi
Psi = cell(2, 2); % the Psi matrix
Psi{1, 1} = Phi;
Psi{1, 2} = B*E;
Psi{2, 1} = zeros(nu*N, nx);
Psi{2, 2} = M;
Psi = cell2mat(Psi);
F_bar = [F + G*K G*E]; % the \bar{F} matrix

% define parameters for geting an initial quantification, i.e., computed
% \hat{W}_0^* using \mathcal{I}_0^w
opts_ini_set.B = B;
opts_ini_set.N_pre_sam = 5; % the is |\mathcal{I}_0^w|, size of initial samples, will be changed in the implementatation
opts_ini_set.N_sam_opt = N_sam_opt;
opts_ini_set.W = W;
opts_ini_set.min_u_LV = min_u_LV;
opts_ini_set.max_u_LV = max_u_LV;
opts_ini_set.Xi_true_EV = Xi_true_EV;
opts_ini_set.Xi_true_LV = Xi_true_LV;
opts_ini_set.nx = nx;
IniSet = InitialSetComputation(opts_ini_set);
[alpha_opt, v_opt] = IniSet.solve_opt(samples_opt');
W_hat_opt = (1 - alpha_opt)*v_opt + alpha_opt*W; % here returns the optimal quantified set, \hat{W}_{opt}, based on vertices of W_true

S_hat_opt = compute_mrpi_set(Phi, W_hat_opt, 1e-2); % RMPI set corresponding to \hat{W}_{opt} (outer approximation)
S_true = compute_mrpi_set(Phi, W_true, 1e-2); % RMPI set corresponding to W_true (outer approximation)
S = compute_mrpi_set(Phi, W, 1e-2); % RMPI set corresponding to W (outer approximation)
S_hat_opt = minHRep(S_hat_opt); % reduce redundant vertices of the set
S_true = minHRep(S_true);       % reduce redundant vertices of the set
S = minHRep(S);                 % reduce redundant vertices of the set
num_half_space_S = length(S.b); % number of half spaces of the set S, will be used in the implementation

% define parameters for computing the feasible region of state s_0|k
opts_feasible_region.S = S;
opts_feasible_region.S_true = S_true;
opts_feasible_region.S_hat_opt = S_hat_opt;
opts_feasible_region.F = F;
opts_feasible_region.G = G;
opts_feasible_region.K = K;
opts_feasible_region.F_bar = F_bar;
opts_feasible_region.Psi = Psi;
opts_feasible_region.nc = nc;
opts_feasible_region.nx = nx;
opts_feasible_region.nu = nu;
opts_feasible_region.N = N;
opts_feasible_region.Ps = ones(nx + N*nu, nx + N*nu);
opts_feasible_region.num_half_space_S = num_half_space_S;
IniFeaReg = ComputeFeasibleRegion(opts_feasible_region);
[F_N_RMPC, hs_RMPC, Nu_RMPC] = IniFeaReg.ComFeasibleRegion_RMPC( ); % the feasible region of RMPC, F_N_RMPC is the feasible region, hs_RMPC is the value of h_s, Nu_RMPC is the value of \nu 
[F_N_True, hs_True] = IniFeaReg.ComFeasibleRegion_True( ); % the true feasible region and hs value with respect to W_true, just for analysis, unknown to the EV
[F_N_Hat_Opt, hs_Hat_Opt] = IniFeaReg.ComFeasibleRegion_UQOPT( ); % the feasible region and hs value with respect to \hat{W}_{opt}, just for analysis, unknown to the EV
opts_feasible_region.hs = hs_RMPC;

% define parameters for modeling the EV and LV
opts_Car.A = A;
opts_Car.B = B;
opts_Car.min_u_LV = min_u_LV;
opts_Car.max_u_LV = max_u_LV;
opts_Car.Xi_true_EV = Xi_true_EV;
opts_Car.Xi_true_LV = Xi_true_LV;

% define parameters of RMPC, i.e., conventional robust MPC
opts_RMPC.N = N;
opts_RMPC.nx = nx;
opts_RMPC.nu = nu;
opts_RMPC.nc = nc;
opts_RMPC.Nu = Nu_RMPC;
opts_RMPC.A = A;
opts_RMPC.B = B;
opts_RMPC.F = F;
opts_RMPC.G = G;
opts_RMPC.K = K;
opts_RMPC.Px = Px;
opts_RMPC.Pc = Pc;
opts_RMPC.F_bar = F_bar;
opts_RMPC.Psi = Psi;
opts_RMPC.Phi = Phi;
opts_RMPC.E = E;
opts_RMPC.W = W;
opts_RMPC.S = S;
opts_RMPC.num_half_space_S = num_half_space_S;
opts_RMPC.hs = hs_RMPC;

% define parameters for UQ-RMPC, i.e., uncertainty-quantification robust
% MPC
opts_UQMPC = opts_RMPC;
opts_UQMPC.Nu = Nu_RMPC; % here the UQ-RMPC does not update \nu at every time step,...

% Computed Ps in (22)
addpath('SDPT3-4.0');
H = [ ];
h = [ ];
for i = 0:1:Nu_RMPC
    H = [H; F_bar*(Psi^i)];
    h = [h; 1 - hs_RMPC];
end
Poly = Polyhedron(H, h);
Poly = minHRep(Poly);
V = Poly.V;
[m, ~] = size(V); % number of vertices
Ps = sdpvar(nx + N*nu, nx + N*nu); % Ps >= 0

Constraints = [ ];
for i = 1:1:m
    Constraints = [Constraints, norm(Ps*V(i, :)', 2) <= 1]; % Boyd, convex optimization, Eq.(8.11)
end
Objective = -logdet(Ps); % log det[A^(-1)]
sol = optimize(Constraints, Objective, sdpsettings('solver','sdpt3'));
Ps = value(Ps);
opts_UQMPC.Ps = Ps;
opts_feasible_region.Ps = Ps;

% save the parameters
parameters.opts_ini_set = opts_ini_set;
parameters.opts_feasible_region = opts_feasible_region;
parameters.opts_Car = opts_Car;
parameters.opts_RMPC = opts_RMPC;
parameters.opts_UQMPC = opts_UQMPC;
parameters.W_hat_opt = W_hat_opt;
parameters.F_N_RMPC = F_N_RMPC;
parameters.F_N_RMPC = F_N_RMPC;
parameters.F_N_True = F_N_True;
parameters.F_N_Hat_Opt = F_N_Hat_Opt;
parameters.x_des = x_des;
parameters.distance_error = distance_error;
parameters.acc_bound = acc_bound;
parameters.T = T;
parameters.W_true = W_true;
save('parameters.mat', 'parameters')



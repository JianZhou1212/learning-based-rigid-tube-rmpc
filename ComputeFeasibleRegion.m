classdef ComputeFeasibleRegion < handle

    properties (SetAccess = public)
        S;
        SA;
        Sb;
        num_half_space_S
        S_true;
        SA_true;
        Sb_true;
        S_hat_opt;
        SA_hat_opt;
        Sb_hat_opt;
        F;
        G;
        K;
        F_bar;
        Psi;
        nc;
        nx;
        nu;
        N;
        Com_hs;
        Com_hs_true;
        Com_hs_hat_opt;
        Com_hs_ini;
        Com_z;
    end

    methods (Access = public)
        function obj = ComputeFeasibleRegion(parameters)
            obj.S  = parameters.S;   % S setA*x <= b
            obj.SA = parameters.S.A; % A matrix of S
            obj.Sb = parameters.S.b; % b vector of S
            obj.num_half_space_S = parameters.num_half_space_S;
            obj.S_true  = parameters.S_true;   % S setA*x <= b
            obj.SA_true = parameters.S_true.A; % A matrix of S
            obj.Sb_true = parameters.S_true.b; % b vector of S
            obj.S_hat_opt  = parameters.S_hat_opt;   % S setA*x <= b
            obj.SA_hat_opt = parameters.S_hat_opt.A; % A matrix of S
            obj.Sb_hat_opt = parameters.S_hat_opt.b; % b vector of S
            obj.F  = parameters.F;
            obj.G  = parameters.G;
            obj.K  = parameters.K;
            obj.F_bar = parameters.F_bar;
            obj.Psi = parameters.Psi;
            obj.nc = parameters.nc;
            obj.nx = parameters.nx;
            obj.nu = parameters.nu;
            obj.N = parameters.N;
            obj.Com_hs = obj.LP( );
            obj.Com_hs_true = obj.LP_true( );
            obj.Com_hs_hat_opt = obj.LP_hat_opt( );
            obj.Com_hs_ini = obj.LP_hs_ini( );
            obj.Com_z = obj.LP_Com_z( );
        end

        function Nu = ComNu_RMPC(obj, hs)
            % compute the value of \nu according to hs
            Nu = obj.N; % Nu must be longer than N
            while Nu
                goal = zeros(obj.nc, 1); % save vector value max{F_bar*(Pis)^(Nu + 1)*z}
                M = obj.F_bar*(obj.Psi^(Nu + 1));
                M_before = obj.F_bar*(obj.Psi^(Nu));
                for i = 1:1:obj.nc
                    goal(i) = full(obj.Com_z(M(i, :), hs, M_before));
                end
                if goal <= (1 - hs) % max{F_bar*(Psi)^(Nu + 1)*z} <= 1 - hs
                    break
                else
                    Nu = Nu + 1;
                    continue
                end
            end
        end

        function [F_N, hs, Nu] = ComFeasibleRegion_RMPC(obj)
            % compute the feasible region of s_0|k with S
            % outputs the feasible region F_N, hs value, and \nu value Nu
            F_Com = obj.F + obj.G*obj.K; % [F + GK]
            hs = zeros(obj.nc, 1);
            for i = 1:1:obj.nc
                hs(i) = full(obj.Com_hs(F_Com(i, :))); % elementwisely get hs vector
            end
            Nu = ComNu_RMPC(obj, hs);

            Sam_num = 5000;
            Sample = 50*rand(2,Sam_num) - 25;
            yalmip('clear')
            z = sdpvar(size(obj.Psi,2), 1); 
            proj_sample = sdpvar(obj.nx, 1);
            Input_sample = sdpvar(obj.nx, 1);
            cns = [];
            for i = 0:1:Nu
                cns = [cns, obj.F_bar*(obj.Psi^i)*z <= ones(obj.nc, 1) - hs];
            end
            cns = [cns,z(1:obj.nx,1) == proj_sample];
            J = (proj_sample-Input_sample)'*(proj_sample-Input_sample);
            ops = sdpsettings('relax', 0);
            Fea_Set = optimizer(cns,J,ops,Input_sample,proj_sample);
            
            % evaluate the function and form a convexhull
            parfor k = 1:Sam_num
                [sample_proj(:,k),~] = Fea_Set(Sample(:,k));
            end

            [convhull_index, ~] = convhull(sample_proj');
            F_Ns = sample_proj(:,convhull_index);
            F_Ns = Polyhedron(F_Ns');
            F_N = F_Ns + obj.S;
        end

        function [F_N, hs_true] = ComFeasibleRegion_True(obj)
            % compute the feasible region of s_0|k with S_true
            % outputs the feasible region F_N, hs value hs_true
            F_Com = obj.F + obj.G*obj.K; % [F + GK]
            hs_true = zeros(obj.nc, 1);
            for i = 1:1:obj.nc
                hs_true(i) = full(obj.Com_hs_true(F_Com(i, :))); % elementwisely get hs vector
            end

            Nu_true = ComNu_RMPC(obj, hs_true);
            
            Sam_num = 5000;
            Sample = 50*rand(2, Sam_num) - 25;
            yalmip('clear')
            z = sdpvar(size(obj.Psi, 2), 1);
            proj_sample =  sdpvar(obj.nx, 1);
            Input_sample = sdpvar(obj.nx, 1);
            cns = [];
            for i = 0:1:Nu_true
                cns = [cns,obj.F_bar*(obj.Psi^i)*z <= ones(obj.nc, 1) - hs_true];
            end
            cns = [cns,z(1:obj.nx) == proj_sample];
            J = (proj_sample - Input_sample)'*(proj_sample - Input_sample);
            ops = sdpsettings('relax', 0);
            Fea_Set = optimizer(cns, J, ops, Input_sample, proj_sample);
            
            % evaluate the function and form a convexhull
            parfor k = 1:Sam_num
                [sample_proj(:,k),~] = Fea_Set(Sample(:,k));
            end
            [convhull_index, ~] = convhull(sample_proj');
            F_Ns = sample_proj(:,convhull_index);
            F_Ns = Polyhedron(F_Ns');
            F_N = F_Ns + obj.S_true;
        end
        
        function [F_N_hat_opt, hs_hat_opt] = ComFeasibleRegion_UQOPT(obj)
            % compute the feasible region of s_0|k with \hat{S}_{opt}
            % outputs the feasible region F_N_hat_opt, the hs value
            % hs_hat_opt
            F_Com = obj.F + obj.G*obj.K; % [F + GK]
            hs_hat_opt = zeros(obj.nc, 1);
            for i = 1:1:obj.nc
                hs_hat_opt(i) = full(obj.Com_hs_hat_opt(F_Com(i, :))); % elementwisely get hs vector
            end

            Nu_true = ComNu_RMPC(obj, hs_hat_opt);
            
            Sam_num = 5000;
            Sample = 50*rand(2, Sam_num) - 25;
            yalmip('clear')
            z = sdpvar(size(obj.Psi, 2), 1);
            proj_sample =  sdpvar(obj.nx, 1);
            Input_sample = sdpvar(obj.nx, 1);
            cns = [];
            for i = 0:1:Nu_true
                cns = [cns,obj.F_bar*(obj.Psi^i)*z <= ones(obj.nc, 1) - hs_hat_opt];
            end
            cns = [cns,z(1:obj.nx) == proj_sample];
            J = (proj_sample - Input_sample)'*(proj_sample - Input_sample);
            ops = sdpsettings('relax', 0);
            Fea_Set = optimizer(cns, J, ops, Input_sample, proj_sample);
            
            % evaluate the function and form a convexhull
            parfor k = 1:Sam_num
                [sample_proj(:,k),~] = Fea_Set(Sample(:,k));
            end
            [convhull_index, ~] = convhull(sample_proj');
            F_Ns_hat_opt = sample_proj(:,convhull_index);
            F_Ns_hat_opt = Polyhedron(F_Ns_hat_opt');
            F_N_hat_opt = F_Ns_hat_opt + obj.S_hat_opt;
        end

        function [F_Ns, F_N] = ComFeasibleRegion_UQMPC(obj, S_hat, h_k, Nu)
            % compute the feasible region of s_0|k of UQMPC with S_hat
            % the region should be time-varying
            % it takes input S_hat and h_k, which are computed online
            Sam_num = 5000;
            Sample = 50*rand(2, Sam_num) - 25;
            yalmip('clear')
            z = sdpvar(size(obj.Psi, 2), 1);
            proj_sample =  sdpvar(obj.nx, 1);
            Input_sample = sdpvar(obj.nx, 1);
            cns = [];
            for i = 0:1:Nu
                cns = [cns,obj.F_bar*(obj.Psi^i)*z <= ones(obj.nc, 1) - h_k];
            end
            cns = [cns,z(1:obj.nx) == proj_sample];
            J = (proj_sample - Input_sample)'*(proj_sample - Input_sample);
            ops = sdpsettings('relax', 0);
            Fea_Set = optimizer(cns, J, ops, Input_sample, proj_sample);
            
            % evaluate the function and form a convexhull
            parfor k = 1:Sam_num
                [sample_proj(:,k),~] = Fea_Set(Sample(:,k));
            end
            [convhull_index, ~] = convhull(sample_proj');
            F_Ns = sample_proj(:,convhull_index);
            F_Ns = Polyhedron(F_Ns');
            F_N = F_Ns + S_hat;
        end
        
        function [hs_ini] = ComFeasibleRegion_UQMPC_Hs_ini(obj, S_hat_ini)
            % given the initial samples, and the estimation of disturbance
            % set and invariant set, compute hs_ini
            F_Com = obj.F + obj.G*obj.K; % [F + GK]
            hs_ini = zeros(obj.nc, 1);
            for i = 1:1:obj.nc
                hs_ini(i) = full(obj.Com_hs_ini(S_hat_ini.A, S_hat_ini.b, F_Com(i, :))); % elementwisely get hs vector
            end
        end

        function Com_hs = LP(obj)
            % LP for computing the hs with S
            % the problem is LP, while the LP solver clp does not work well
            % the NLP solver ipopt is accurate enough, takes longer time
            opti = casadi.Opti( );
            F_com = opti.parameter(1, obj.nx);
            e = opti.variable(obj.nx, 1);

            opti.minimize(-F_com*e);
            opti.subject_to(obj.SA*e <= obj.Sb); % e \in S

            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            Com_hs = opti.to_function('f', {F_com}, {F_com*e});
        end

        function Com_hs_true = LP_true(obj)
            % LP for computing the hs_true with S_true
            % the problem is LP, while the LP solver clp does not work well
            % the NLP solver ipopt is accurate enough, takes longer time
            opti = casadi.Opti( );
            F_com = opti.parameter(1, obj.nx);
            e = opti.variable(obj.nx, 1);

            opti.minimize(-F_com*e);
            opti.subject_to(obj.SA_true*e <= obj.Sb_true); % e \in S

            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            Com_hs_true = opti.to_function('f', {F_com}, {F_com*e});
        end
        
        function Com_hs_hat_opt = LP_hat_opt(obj)
            % LP for computing the \hat{hs}_opt with S_true
            % the problem is LP, while the LP solver clp does not work well
            % the NLP solver ipopt is accurate enough, takes longer time
            opti = casadi.Opti( );
            F_com = opti.parameter(1, obj.nx);
            e = opti.variable(obj.nx, 1);

            opti.minimize(-F_com*e);
            opti.subject_to(obj.SA_hat_opt*e <= obj.Sb_hat_opt); % e \in S

            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            Com_hs_hat_opt = opti.to_function('f', {F_com}, {F_com*e});
        end
        
        function Com_hs_ini = LP_hs_ini(obj)
            % LP for computing the hs with S
            % the problem is LP, while the LP solver clp does not work well
            % the NLP solver ipopt is accurate enough, takes longer time
            opti = casadi.Opti( );
            F_com = opti.parameter(1, obj.nx);
            S_A  = opti.parameter(obj.num_half_space_S, obj.nx);
            S_b  = opti.parameter(obj.num_half_space_S, 1);
            e = opti.variable(obj.nx, 1);

            opti.minimize(-F_com*e);
            opti.subject_to(S_A*e <= S_b); % e \in S

            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            Com_hs_ini = opti.to_function('f', {S_A, S_b, F_com}, {F_com*e});
        end

        function Com_z = LP_Com_z(obj)
            % LP for computing the element of max{F_bar*(Pis)^(Nu + 1)*z}
            % the problem is LP, while the LP solver clp does not work well
            % the NLP solver ipopt is accurate enough, takes longer time
            opti = casadi.Opti( );
            M_i = opti.parameter(1, obj.nx + obj.N*obj.nu); % the i-th row of the matrix F_bar*(Psi)^(Nu + 1)
            hs = opti.parameter(obj.nc, 1);
            M_before = opti.parameter(obj.nc, obj.nx + obj.N*obj.nu);
            z = opti.variable(obj.nx + obj.N*obj.nu, 1);

            opti.minimize(-M_i*z);
            opti.subject_to(M_before*z <= ones(obj.nc, 1) - hs); % z \in \Omega(1-hs, Nu_s)

            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            Com_z = opti.to_function('f', {M_i, hs, M_before}, {M_i*z});
        end
        
    end
end

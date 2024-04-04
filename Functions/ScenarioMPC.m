classdef ScenarioMPC < handle
    
    properties (SetAccess = public)
        A;
        B;
        N_SC;
        Q;
        R;
        N;
        F;
        G;
        nu;
        nx;
        min_u_LV;
        max_u_LV;
        Xi_true_EV;
        Xi_true_EV_A;
        Xi_true_EV_b;
        Xi_true_LV;
        Xi_true_LV_A;
        Xi_true_LV_b;
        optimizer;
    end
    
    methods (Access = public)
        function obj = ScenarioMPC(parameters)
            obj.A = parameters.A;
            obj.B = parameters.B;
            obj.N_SC = parameters.N_SC;
            obj.Q = parameters.Q;
            obj.R = parameters.R;
            obj.N = parameters.N;
            obj.F = parameters.F;
            obj.G = parameters.G;
            obj.nu = parameters.nu;
            obj.nx = parameters.nx;
            obj.min_u_LV = parameters.min_u_LV;
            obj.max_u_LV = parameters.max_u_LV;
            obj.Xi_true_EV = parameters.Xi_true_EV;
            obj.Xi_true_EV_A = parameters.Xi_true_EV.A;
            obj.Xi_true_EV_b = parameters.Xi_true_EV.b;
            obj.Xi_true_LV = parameters.Xi_true_LV;
            obj.Xi_true_LV_A = parameters.Xi_true_LV.A;
            obj.Xi_true_LV_b = parameters.Xi_true_LV.b;
            obj.optimizer = obj.MPCFormulation( );
        end

        function [u_k] = solve(obj, x_k)
            scenarios = obj.sample_scenarios( );
            x_k_full = repmat(x_k, obj.N_SC, 1);
            u_opt = obj.optimizer(x_k_full, scenarios);
            u_opt = full(u_opt);
            u_k   = u_opt(:, 1);
        end

        function scenarios = sample_scenarios(obj)
            scenarios = ones(obj.N_SC*obj.nx, obj.N);
            for k = 1:1:obj.N_SC
                for i = 1:1:obj.N
                    u_LV_random = (obj.max_u_LV - obj.min_u_LV)*rand(1) + obj.min_u_LV;
                    scenarios((k-1)*obj.nx + 1:k*obj.nx, i) = polytope_sample_EV(obj, 1) - polytope_sample_LV(obj, 1) - obj.B*u_LV_random;
                end
            end
        end

        function optimizer = MPCFormulation(obj)
            opti = casadi.Opti( );
            U = opti.variable(1, obj.N);
            X = opti.variable(obj.N_SC*obj.nx, obj.N + 1);

            x_k_full = opti.parameter(obj.N_SC*obj.nx, 1);
            scenarios = opti.parameter(obj.N_SC*obj.nx, obj.N);

            X(:, 1) = x_k_full;

            J = 0;

            for k = 1:1:obj.N_SC
                for i = 1:1:obj.N
                    x_current = X((k-1)*obj.nx + 1:k*obj.nx, i);
                    u_current = U(:, i);
                    w_current = scenarios((k-1)*obj.nx + 1:k*obj.nx, i);
                    x_next = obj.A*x_current + obj.B*u_current + w_current;
                    
                    X((k-1)*obj.nx + 1:k*obj.nx, i + 1) = x_next;
                    opti.subject_to(obj.F*x_next + obj.G*u_current <= 1);

                    J = J + x_next'*obj.Q*x_next + u_current'*obj.R*u_current;
                end
            end
            
            opti.minimize(J);

            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            optimizer = opti.to_function('f', {x_k_full, scenarios}, {U});
        end
        
        function EV_samples = polytope_sample_EV(obj, N_sam)
            V = obj.Xi_true_EV.V;
            ver_x = V(:, 1);
            ver_y = V(:, 2);
            min_x = min(ver_x);
            max_x = max(ver_x);
            min_y = min(ver_y);
            max_y = max(ver_y);
            EV_samples = zeros(2, N_sam);
            i = 1;
            while i <= N_sam
                x = (max_x-min_x).*rand(1) + min_x;
                y = (max_y-min_y).*rand(1) + min_y;
                if obj.Xi_true_EV_A*[x; y] <= obj.Xi_true_EV_b
                    EV_samples(:, i) = [x; y];
                    i = i + 1;
                else
                    continue;
                end
            end
        end
        
        function LV_samples = polytope_sample_LV(obj, N_sam)
            V = obj.Xi_true_LV.V;
            ver_x = V(:, 1);
            ver_y = V(:, 2);
            min_x = min(ver_x);
            max_x = max(ver_x);
            min_y = min(ver_y);
            max_y = max(ver_y);
            LV_samples = zeros(2, N_sam);
            i = 1;
            while i <= N_sam
                x = (max_x-min_x).*rand(1) + min_x;
                y = (max_y-min_y).*rand(1) + min_y;
                if obj.Xi_true_LV_A*[x; y] <= obj.Xi_true_LV_b
                    LV_samples(:, i) = [x; y];
                    i = i + 1;
                else
                    continue;
                end
            end
        end
        
    end
end

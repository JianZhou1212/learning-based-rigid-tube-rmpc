classdef InitialSetComputation < handle
    
    properties (SetAccess = public)
        B;
        N_pre_sam;
        N_sam_opt;
        nx;
        W;
        W_A;
        W_b;
        min_u_LV;
        max_u_LV;
        Xi_true_EV;
        Xi_true_EV_A;
        Xi_true_EV_b;
        Xi_true_LV;
        Xi_true_LV_A;
        Xi_true_LV_b;
        com_ini_set;
        com_ini_set_opt;
    end
    
    methods (Access = public)
        function obj = InitialSetComputation(parameters)
            obj.B = parameters.B;
            obj.N_pre_sam = parameters.N_pre_sam;
            obj.N_sam_opt = parameters.N_sam_opt;
            obj.nx = parameters.nx;
            obj.W = parameters.W;
            obj.W_A = parameters.W.A;
            obj.W_b = parameters.W.b;
            obj.min_u_LV = parameters.min_u_LV;
            obj.max_u_LV = parameters.max_u_LV;
            obj.Xi_true_EV = parameters.Xi_true_EV;
            obj.Xi_true_EV_A = parameters.Xi_true_EV.A;
            obj.Xi_true_EV_b = parameters.Xi_true_EV.b;
            obj.Xi_true_LV = parameters.Xi_true_LV;
            obj.Xi_true_LV_A = parameters.Xi_true_LV.A;
            obj.Xi_true_LV_b = parameters.Xi_true_LV.b;
            obj.com_ini_set = obj.LP( );
            obj.com_ini_set_opt = obj.LP_opt( );
        end

        function [alpha_ini, v_ini, samplesinput] = solve(obj)
            samplesinput = ones(obj.nx, obj.N_pre_sam);
            for i = 1:1:obj.N_pre_sam
                u_LV_random = (obj.max_u_LV - obj.min_u_LV)*rand(1) + obj.min_u_LV;
                samplesinput(:, i) = polytope_sample_EV(obj, 1) - polytope_sample_LV(obj, 1) - obj.B*u_LV_random;
            end
            out = full(obj.com_ini_set(samplesinput));
            beta_ini = out(1);
            y_ini = out(2:end);
            alpha_ini  = 1 - beta_ini;
            v_ini      = y_ini/beta_ini;
        end
        
        function [alpha_opt, v_opt] = solve_opt(obj,samples_opt)
            out = full(obj.com_ini_set_opt(samples_opt));
            beta_opt = out(1);
            y_opt = out(2:end);
            alpha_opt  = 1 - beta_opt;
            v_opt      = y_opt/beta_opt;
        end
        
        function com_ini_set = LP(obj)
            opti = casadi.Opti( );
            var = opti.variable(3, 1);
            beta = var(1);
            y = var(2:end);
            samplesinput = opti.parameter(obj.nx, obj.N_pre_sam);
            opti.minimize(-beta);
            opti.subject_to(0 <= beta <= 1);
            for i = 1:1:obj.N_pre_sam
                opti.subject_to(obj.W_A*(samplesinput(:, i) - y) <= (1 - beta)*obj.W_b);
            end
            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            com_ini_set = opti.to_function('f', {samplesinput}, {var});
        end
        
        function com_ini_set_opt = LP_opt(obj)
            opti = casadi.Opti( );
            var = opti.variable(3, 1);
            beta = var(1);
            y = var(2:end);
            samplesinput = opti.parameter(obj.nx, obj.N_sam_opt);
            opti.minimize(-beta);
            opti.subject_to(0 <= beta <= 1);
            for i = 1:1:obj.N_sam_opt
                opti.subject_to(obj.W_A*(samplesinput(:, i) - y) <= (1 - beta)*obj.W_b);
            end
            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            com_ini_set_opt = opti.to_function('f', {samplesinput}, {var});
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

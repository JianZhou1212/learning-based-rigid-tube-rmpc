classdef UQRobustMPC < handle
    
    properties (SetAccess = public)
        N;
        nx;
        nu;
        nc;
        F;
        G;
        K;
        Px;
        Pc;
        F_bar;
        Psi;
        Phi;
        E;
        W;
        W_A;
        W_b;
        S;
        SA;
        Sb;
        num_half_space_S;
        Nu;
        optimizer;
        Com_hs;
        UQ;
        
    end
    
    methods (Access = public)
        function obj = UQRobustMPC(parameters)
            obj.N = parameters.N;
            obj.nx = parameters.nx;
            obj.nu = parameters.nu;
            obj.nc = parameters.nc;
            obj.F = parameters.F;
            obj.G = parameters.G;
            obj.K = parameters.K;
            obj.Px = parameters.Px;
            obj.Pc = parameters.Pc;
            obj.F_bar = parameters.F_bar;
            obj.Psi = parameters.Psi;
            obj.Phi = parameters.Phi;
            obj.E = parameters.E;
            obj.W = parameters.W;
            obj.W_A = parameters.W.A;
            obj.W_b = parameters.W.b;
            obj.S = parameters.S;
            obj.SA = parameters.S.A;
            obj.Sb = parameters.S.b;
            obj.num_half_space_S = parameters.num_half_space_S;
            obj.Nu = parameters.Nu;
            obj.optimizer = obj.MPCFormulation( );
            obj.Com_hs = obj.LP( );
            obj.UQ = obj.LPforUQ( );
        end

        function [s_k_opt, u_k, alpha_k, v_k, hs, S_hat_star, W_hat] = solve(obj, x_k, alpha_before, v_before, w_new)
            out = full(obj.UQ(alpha_before, v_before, w_new));
            beta_k = out(1);
            y_k = out(2:end);
            alpha_k  = 1 - beta_k;
            v_k      = y_k/beta_k;
            
            W_hat = (1 - alpha_k)*v_k + alpha_k*obj.W;
            S_hat_star = alpha_k*obj.S + (1 - alpha_k)*(inv(1 - obj.Phi)*v_k);
            S_A = S_hat_star.A;
            S_b = S_hat_star.b;
        
            F_Com = obj.F + obj.G*obj.K;
            hs = zeros(obj.nc, 1);
            for i = 1:1:obj.nc
                hs(i) = full(obj.Com_hs(F_Com(i, :), S_A, S_b));
            end
            var_opt = full(obj.optimizer(x_k, hs, S_A, S_b));
            s_k_opt = var_opt(1:obj.nx);
            c_k_opt = var_opt(obj.nx+1:end);
            u_k = obj.K*x_k + c_k_opt(1:obj.nu);
        end
        
        function Com_hs = LP(obj)
            opti = casadi.Opti( );
            F_com = opti.parameter(1, obj.nx);
            S_A  = opti.parameter(obj.num_half_space_S, obj.nx);
            S_b  = opti.parameter(obj.num_half_space_S, 1);
            e = opti.variable(obj.nx, 1);
            
            opti.minimize(-F_com*e);
            opti.subject_to(S_A*e <= S_b);
            
            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            Com_hs = opti.to_function('f', {F_com, S_A, S_b}, {F_com*e});
        end
        
        function UQ = LPforUQ(obj)
            opti = casadi.Opti( );
            var = opti.variable(3, 1);
            beta = var(1);
            y = var(2:end);
            alpha_before = opti.parameter( );
            v_before = opti.parameter(obj.nx, 1);
            w_new = opti.parameter(obj.nx, 1);
            opti.minimize(-beta);
            opti.subject_to(0 <= beta <= 1);
            opti.subject_to(obj.W_A*y - beta*obj.W_b <= 0);
            opti.subject_to(alpha_before*obj.W_b + (1-alpha_before)*(obj.W_A*v_before) <= (1 - beta)*obj.W_b + obj.W_A*y);
            opti.subject_to(obj.W_A*w_new <= (1 - beta)*obj.W_b + obj.W_A*y);
            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            UQ = opti.to_function('f', {alpha_before, v_before, w_new}, {var});
        end
        
        function optimizer = MPCFormulation(obj)
            opti = casadi.Opti( );
            var = opti.variable(obj.nx + obj.N*obj.nu, 1);
            s_0k = var(1:obj.nx);
            c_k = var(obj.nx+1:end);
            x_k = opti.parameter(obj.nx, 1);
            hs = opti.parameter(obj.nc, 1);
            S_A  = opti.parameter(obj.num_half_space_S, obj.nx);
            S_b  = opti.parameter(obj.num_half_space_S, 1);
            J = s_0k'*obj.Px*s_0k + c_k'*obj.Pc*c_k;
            
            opti.minimize(J);
            opti.subject_to(S_A*(x_k - s_0k) <= S_b);
            for i = 0:1:obj.Nu
                opti.subject_to(obj.F_bar*(obj.Psi^i)*var <= ones(obj.nc, 1) - hs);
            end
            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            optimizer = opti.to_function('f', {x_k, hs, S_A, S_b}, {var});
        end
    end
end

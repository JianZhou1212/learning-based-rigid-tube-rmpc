classdef NominalRobustMPC < handle
    
    properties (SetAccess = public)
        N;
        nx;
        nu;
        nc;
        Nu;
        Px;
        Pc;
        K;
        F_bar;
        Psi;
        E;
        S;
        SA;
        Sb;
        hs;
        optimizer;
        
    end
    
    methods (Access = public)
        function obj = NominalRobustMPC(parameters)
            obj.N = parameters.N;
            obj.nx = parameters.nx;
            obj.nu = parameters.nu;
            obj.nc = parameters.nc;
            obj.Nu = parameters.Nu;
            obj.Px = parameters.Px;
            obj.Pc = parameters.Pc;
            obj.K = parameters.K;
            obj.F_bar = parameters.F_bar;
            obj.Psi = parameters.Psi;
            obj.E = parameters.E;
            obj.S = parameters.S;
            obj.SA = parameters.S.A;
            obj.Sb = parameters.S.b;
            obj.hs = parameters.hs;
            obj.optimizer = obj.MPCFormulation( );
        end

        function [s_k_opt, u_k] = solve(obj, x_k)
            var_opt = full(obj.optimizer(x_k));
            s_k_opt = var_opt(1:obj.nx);
            c_k_opt = var_opt(obj.nx+1:end);
            u_k = obj.K*x_k + c_k_opt(1:obj.nu);
        end
        
        function optimizer = MPCFormulation(obj)
            opti = casadi.Opti( );
            var = opti.variable(obj.nx + obj.N*obj.nu, 1);
            s_0k = var(1:obj.nx);
            c_k = var(obj.nx+1:end);
            x_k = opti.parameter(obj.nx, 1);
            J = s_0k'*obj.Px*s_0k + c_k'*obj.Pc*c_k;
            
            opti.minimize(J);
            opti.subject_to(obj.SA*(x_k - s_0k) <= obj.Sb);
            for i = 0:1:obj.Nu
                opti.subject_to(obj.F_bar*(obj.Psi^i)*var <= ones(obj.nc, 1) - obj.hs);
            end
            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            optimizer = opti.to_function('f', {x_k}, {var});
        end
    end
end

classdef ModelingCar < handle
    
    properties (SetAccess = public)
        A;
        B;
        Xi_true_EV;
        Xi_true_EV_A;
        Xi_true_EV_b;
        Xi_true_LV;
        Xi_true_LV_A;
        Xi_true_LV_b;
        max_u_LV;
        min_u_LV;
    end
    
    methods (Access = public)
        function obj = ModelingCar(parameters)
            obj.A = parameters.A;
            obj.B = parameters.B;
            obj.Xi_true_EV = parameters.Xi_true_EV;
            obj.Xi_true_EV_A = parameters.Xi_true_EV.A;
            obj.Xi_true_EV_b = parameters.Xi_true_EV.b;
            obj.Xi_true_LV = parameters.Xi_true_LV;
            obj.Xi_true_LV_A = parameters.Xi_true_LV.A;
            obj.Xi_true_LV_b = parameters.Xi_true_LV.b;
            obj.max_u_LV = parameters.max_u_LV;
            obj.min_u_LV = parameters.min_u_LV;
        end

        function [xi_EV_k, x_EV_k_next] = EVModeling(obj, x_EV_k, u_EV_k)
            % system dynamic of the EV
            xi_EV_k = polytope_sample_EV(obj, 1); % random sample from \Xi_EV_true
            x_EV_k_next = obj.A*x_EV_k + obj.B*u_EV_k + xi_EV_k;
        end
        
        function [xi_LV_k, x_LV_k_next] = LVModeling(obj, x_LV_k)
            % system dynamic of the LV
            xi_LV_k = polytope_sample_LV(obj, 1); % random sample from \Xi_LV_true
            u_LV_random = (obj.max_u_LV - obj.min_u_LV)*rand(1) + obj.min_u_LV; % random sample from U_LV_true
            x_LV_k_next = obj.A*x_LV_k + xi_LV_k + obj.B*u_LV_random;
        end
        
        function EV_samples = polytope_sample_EV(obj, N_sam)
            % generate samples from \Xi_EV_true, N_sam is the number of
            % samples
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
            % generate samples from \Xi_LV_true, N_sam is the number of
            % samples
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

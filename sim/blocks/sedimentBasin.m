classdef sedimentBasin
    %SEDIMENTATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        maxdp       % Maximum particle diameter which will settle (m)
        Q           % Volumetric flow rate (m3/s)
    
        model       % Model reference
    end
    
    properties (Dependent)
        V           % Sedimentation basin volume (m3)
    end
    
    properties (Constant)
        RHOp=1050   % Influent particle density (kg/m3)
    end
    
    methods
        function obj = sedimentBasin(maxdp, Q, model)
            %SEDIMENTATION Construct an instance of this class
            %   Detailed explanation goes here
            obj.maxdp = maxdp;
            obj.Q = Q;
            obj.model = model;
        end
        
        function [L, W, D] = dimensions(obj)
            % Using L:W of 5:1 and D=4m
            W = sqrt(obj.A/5);
            L = 5*W;
            D = 4;
        end
        
        function A = A(obj)
            RHOw = refEQ.RHOw(obj.model.Ta_C);
            mu = refEQ.mu(obj.model.Ta_C);
            
            vs = 9.81*(obj.RHOp-RHOw)*obj.maxdp^2/(18*mu);
            % Check if flow conditions are correct (i.e. Re < 2)
            Re = RHOw*vs*obj.maxdp/mu;
            % Recalculate vs if Re > 2
            if Re > 2
                vs = (9.81*(obj.RHOp-RHOw)*obj.maxdp^1.6/(13.9*RHOw^0.4*mu^0.6))^(1/1.4);
            end
            
            A = obj.Q/vs;
        end
        
        function V = get.V(obj)
            [L, W, D] = obj.dimensions;
            V = L*W*D;
        end
    end
end


classdef sedimentBasin
    %SEDIMENTATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        maxdp       % Maximum particle diameter which will settle (m)
        Q           % Volumetric flow rate (m3/d)
    
        model       % Model reference
    end
    
    properties (Dependent)
        A           % Sedimentation basin area (m2)
        V           % Sedimentation basin volume (m3)
        TSS         % Effluent TSS distribution (mg/L)
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
        
        function A = get.A(obj)
            RHOp = obj.model.WW_RHOp;
            RHOw = refEQ.RHOw(obj.model.Ta_C);
            mu = refEQ.mu(obj.model.Ta_C);
            
            vs = 9.81*(RHOp-RHOw)*obj.maxdp^2/(18*mu);
            % Check if flow conditions are correct (i.e. Re < 2)
            Re = RHOw*vs*obj.maxdp/mu;
            % Recalculate vs if Re > 2
            if Re > 2
                vs = (9.81*(RHOp-RHOw)*obj.maxdp^1.6/(13.9*RHOw^0.4*mu^0.6))^(1/1.4);
            end
            
            A = obj.Q/3600/vs;
        end
        
        function V = get.V(obj)
            [L, W, D] = obj.dimensions;
            V = L*W*D;
        end
        
        function TSS = get.TSS(obj)
            % Calculate the effluent TSS based on influent particles
            % diameters
            dp = obj.model.WW_dp;
            TSS = zeros(1, length(dp));
            
            % Calculate stoke's settling velocity for maxdp
            RHOp = obj.model.WW_RHOp;
            RHOw = refEQ.RHOw(obj.model.Ta_C);
            mu = refEQ.mu(obj.model.Ta_C);

            vd = 9.81*(RHOp-RHOw)*obj.maxdp^2/(18*mu);
            % Check if flow conditions are correct (i.e. Re < 2)
            Re = RHOw*vd*obj.maxdp/mu;
            % Recalculate vs if Re > 2
            if Re > 2
                vd = (9.81*(RHOp-RHOw)*obj.maxdp^1.6/(13.9*RHOw^0.4*mu^0.6))^(1/1.4);
            end
            
            i=1;
            for d=dp
                if d > obj.maxdp
                    i = i + 1;
                    continue
                end
                % Calculate the Stoke's settling velocity for each d
                RHOw = refEQ.RHOw(obj.model.Ta_C);
                mu = refEQ.mu(obj.model.Ta_C);

                vs = 9.81*(RHOp-RHOw)*d^2/(18*mu);
                % Check if flow conditions are correct (i.e. Re < 2)
                Re = RHOw*vs*d/mu;
                % Recalculate vs if Re > 2
                if Re > 2
                    vs = (9.81*(RHOp-RHOw)*d^1.6/(13.9*RHOw^0.4*mu^0.6))^(1/1.4);
                end
                
                R = vs/vd;
                TSS(i) = R;
                
                i = i + 1;
            end
            
            TSS_0 = obj.model.WW_TSS_0;
            dist  = obj.model.WW_dist;
            TSS = TSS_0*dist.*TSS;
        end
    end
end


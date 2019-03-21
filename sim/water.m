classdef water < handle
    %WATER Tracks the flow of water throughout the system
    %   Water is compartmentalized in this class into
    %   different areas throughout the system.
    
    properties
        Vcrop       % Volume of water in soil (m3)
        Vwetland    % Volume of water in vertical wetlands (m3)
        Vvapor      % Volume of water in air as vapor (m3)
        Vstore      % Volume of water in water bodies (m3)
        Vsed        % Volume of water in the sediment basin (m3)
        Vreject     % Volume of water rejected from RO unit (m3)
        
        % Temporary variables (valid for one time step only)
        ET          % Crop evapotranspiration (m3)
        WW          % Wastewater flow (m3)
        evap        % Evaporated water (m3)
        cond        % Condensated water (m3)
        VRO         % RO flow (m3)
        ROreject    % RO reject (m3)
        dVwetland   % Exiting wetland (m3)
        
        model       % Handle for main model
    end
    
    properties (Dependent)
        V           % Total volume of water in system (m3)
        M           % Total mass of water in the system (kg)
        Mcrop       % Mass of water in soil (kg)
        Mwetland    % Mass of water in vertical wetlands (kg)
        Mvapor      % Mass of water in air as vapor (kg)
        Mstore      % Mass of water in water bodies (kg)
        Msed        % Mass of water in the sediment basin (kg)
        Mreject     % Mass of water rejected from RO unit (kg)
    end
    
    methods
        function obj = water(model)
            %WATER Construct an instance of this class
            %   Mtank - initial water amount in storage
            %   Mbody - initial water amount in open water bodies
            obj.model = model;
            Vstore = model.waterBlock.Vw;

            % Make some initial assumptions about the distribution of water
            % in the system
            obj.Vstore = Vstore;
            obj.Vvapor = 0;
            obj.Vcrop = 0;
            obj.Vwetland = 0;
            obj.Vsed = 0;
            obj.Vreject = 0;
        end
        
        function step(obj)
            % Reset all temporary variables
            obj.WW = 0; obj.ET = 0; obj.evap = 0; obj.cond = 0;
            obj.ROreject = 0; obj.dVwetland = 0;
            
            % Run once every hour
            m = obj.model;
            dt = 3600/m.dt;
            
            % wastewater flux
            obj.WW = m.WW_Q/(24*dt);
            obj.Vstore = obj.Vstore - obj.WW;
            
            % waterBody heat exchange
            if ~isempty(m.dMv) && m.dMv ~=0
                if m.dMv > 0
                    obj.evap = obj.evap + m.dMv / 1000;
                else
                    obj.cond = obj.cond + m.dMv / 1000;
                end
                obj.Vvapor = obj.Vvapor + m.dMv/1000;
                obj.Vstore = obj.Vstore - m.dMv/1000;
                m.dMv = 0;
            end
            
            % crop irrigation/ET
            if ~isempty(m.Vi) && m.Vi ~= 0
                obj.ET = m.Vi;
                obj.Vstore = obj.Vstore - m.Vi;
                obj.Vvapor = obj.Vvapor + m.Vi;
                m.Vi = 0;
            end
            
            % sedimentation basin
            obj.dVwetland = (m.wetlandBlock.qd*m.wetlandBlock.A)/(24*dt);
            if obj.Vsed == 0
                obj.Vsed = m.sedBlock.V;
                obj.Vstore = obj.Vstore - m.sedBlock.V;
            else
                obj.Vsed = obj.Vsed + obj.WW;
                obj.Vsed = obj.Vsed - obj.dVwetland;
                obj.Vwetland = obj.Vwetland + obj.dVwetland;
            end
            
            % wetland
            if obj.Vwetland == 0
                obj.Vwetland = m.wetlandBlock.Vw;
                obj.Vstore = obj.Vstore - m.wetlandBlock.Vw;
            else
                obj.Vwetland = obj.Vwetland - obj.dVwetland;
            end
            
            % filtration treatment
            no_RO = ceil(mod(m.WW_Q/24, filterUnit.RO_min_flux));
            obj.VRO = max(m.WW_Q/24,filterUnit.RO_min_flux*no_RO)/dt;
            obj.Vstore = obj.Vstore - max(obj.VRO - obj.dVwetland, 0);
            recoveryEfficiency = 0.75;
            obj.ROreject = obj.VRO * (1-recoveryEfficiency);
            
            % Track amount of reject to track salt.
            % Not included in total amount of water in system
            obj.Vreject = obj.ROreject;
            
            % Evaporate/distill the reject water and transfer clean to storage
            obj.Vstore = obj.Vstore + obj.VRO;
            
            % Update storage height
            m.waterBlock.hw = obj.Vstore / m.waterBlock.Aw;
            
            % The amount actually evaporated cannot exceed
            % the desired RH of the air
            RHOvs = obj.model.RHOvs * obj.model.RHmax;
            Vsat = RHOvs * obj.model.Volume / 1000;
            if obj.Vvapor > Vsat
                Vdew = obj.Vvapor - Vsat;
                obj.Vstore = obj.Vstore + Vdew;
                obj.Vvapor = obj.Vvapor - Vdew;
                obj.cond = obj.cond + Vdew;
                
                % Remove heat from air due to condensation
%                 Ld = refEQ.L(obj.model.Ta_C);
%                 qd = Ld * Vdew * 1000;
%                 obj.model.Ta_C = obj.model.Ta_C - qd/...
%                     (obj.Mvapor*obj.model.cv + obj.model.Ma*obj.model.ca);
            end
        end
        
        function V = get.V(obj)
            V = obj.Vstore + obj.Vvapor + obj.Vcrop + obj.Vwetland + ...
                obj.Vsed;
        end
        
        function M = get.M(obj)
            M = obj.V * 1000;
        end
        
        function Mstore = get.Mstore(obj)
            Mstore = obj.Vstore * 1000;
        end
        
        function Mvapor = get.Mvapor(obj)
            Mvapor = obj.Vvapor * 1000;
        end
        
        function Mcrop = get.Mcrop(obj)
            Mcrop = obj.Vcrop * 1000;
        end
        
        function Mwetland = get.Mwetland(obj)
            Mwetland = obj.Vwetland *1000;
        end
        
        function Mreject = get.Mreject(obj)
            Mreject = obj.Vreject * 1000;
        end
    end
end


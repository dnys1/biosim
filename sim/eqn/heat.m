classdef heat < handle
    %HEAT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Ts_C        % Water body surface temperature (C)
        
        model       % Main model handle
    end
    
    properties (Constant)
        sigma = 5.670367e-8;    % Stefan-Boltzmann constant (W m^2 K^-4)
        A = 0.6;                % Atmospheric attenuation constant (between 0.5-0.7)
    end
    
    properties (Dependent)
        Qan         % Atmospheric long-wave heat
        Qsn         % Short-wave heat
        Qbr         % Blackbody radiation
        Qc          % Conduction/convection heat
        Qe          % Evaporative heat
        dMv         % Mass of water evaporated(+)/condensed(-)
        
        Aw          % Area of water where heat transfer occurs (m2)
        wv          % Wind velocity (m/s)
        Pv          % Water vapor pressure (kPa)
        dt          % Time delta (s)
        Ta_C        % Air temperature (C)
        Ts_K        % Surface temperature (K)
        Ta_K        % Air temperature (K)
    end
    
    methods
        function obj = heat(model)
            %HEAT Construct an instance of this class
            %   Detailed explanation goes here
            obj.model = model;
        end
        
        function Q = water(obj, Ts_C, Qsn)
            obj.Ts_C = Ts_C;
            Q = obj.Qan + Qsn - (obj.Qbr + obj.Qe + obj.Qc);
        end
        
        function Qan = get.Qan(obj)
            Qan = obj.sigma * obj.Aw * obj.model.Ta_K^4 * (obj.A + 0.031 * sqrt(obj.Pv/0.133322)) * (1 - 0.03) * obj.dt;
        end
        
        function Qbr = get.Qbr(obj)
            Qbr = 0.97 * obj.sigma * obj.Aw * obj.Ts_K^4 * obj.dt;
        end
        
        function Qc = get.Qc(obj)
            if obj.Ts_C > obj.Ta_C
                hc = 10.45 - obj.wv + 10*sqrt(obj.wv);
                Qc = hc * obj.Aw * (obj.Ts_C - obj.Ta_C) * obj.dt;
            else
                L = sqrt(obj.Aw);
                kWater = 600e-3;
                muWater=refEQ.mu(obj.Ts_C);
                Cp = 4150;
                Pr = muWater*Cp/kWater;
                Re = refEQ.RHOw(obj.Ts_C)*obj.wv*L/muWater;
                h = (0.037*Re^(4/5)-871)*Pr^(1/3)/L*kWater;
                Qc = h * obj.Aw * (obj.Ts_C - obj.Ta_C) * obj.dt;
            end
        end
        
        function Qe = get.Qe(obj)
            L = refEQ.L(obj.Ts_C);
            Qe = L * obj.dMv * 1000;
        end
        
        function dMv = get.dMv(obj)
            if obj.Ta_C < obj.Ts_C
                RHOvs = refEQ.Pvs(obj.Ta_C) / (refEQ.Rv * obj.Ta_K);
            else
                RHOvs = refEQ.Pvs(obj.Ts_C) / (refEQ.Rv * obj.Ts_K);
            end
            
            x_s = RHOvs * obj.model.Volume / obj.model.M;
            x = obj.model.W.Mvapor / obj.model.M;
            ae = (25 + 19*obj.wv) / 3600;
            dMv = obj.Aw * (x_s - x) * ae * obj.dt;
        end
        
        function Aw = get.Aw(obj)
            Aw = obj.model.waterBlock.Aw;
        end
        
        function wv = get.wv(obj)
            wv = obj.model.wv;
        end
        
        function Pv = get.Pv(obj)
            Pv = obj.model.Pv;
        end
        
        function dt = get.dt(obj)
            dt = obj.model.dt;
        end
        
        function Ta_C = get.Ta_C(obj)
            Ta_C = obj.model.Ta_C;
        end
        
        function Ts_K = get.Ts_K(obj)
            Ts_K = obj.Ts_C + 273.15;
        end
        
        function Ta_K = get.Ta_K(obj)
            Ta_K = obj.Ta_C + 273.15;
        end
    end
end
classdef heat < handle
    %HEAT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Ta_C        % Air temperature (C)
        Ts_C        % Water body surface temperature (C)
        Pv          % Water vapor pressure (kPa)
        Aw          % Water body area (m2)
        qm          % Quantity of water evaporated (kg)
        wv          % Wind velocity (m/s)
        dt          % Time delta (s)
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
    end
    
    methods
        function obj = heat(dt)
            %HEAT Construct an instance of this class
            %   Detailed explanation goes here
            obj.dt = dt;
        end
        
        function Q = water(obj, Ts_C, Ta_C, Pv, qm, Qsn, wv, Aw)
            obj.Ta_C = Ta_C;
            obj.Ts_C = Ts_C;
            obj.Pv = Pv;
            obj.Aw = Aw;
            obj.qm = qm;
            obj.wv = wv;
            Q = obj.Qan + Qsn - (obj.Qbr + obj.Qe + obj.Qc);
        end
        
        function Qan = get.Qan(obj)
            Qan = obj.sigma * obj.Aw * (obj.Ta_C + 273.15).^4 * (obj.A + 0.031 * sqrt(obj.Pv/0.133322)) * (1 - 0.03) * obj.dt;
        end
        
        function Qbr = get.Qbr(obj)
            Qbr = 0.97 * obj.sigma * obj.Aw * (obj.Ts_C + 273.15).^4 * obj.dt;
        end
        
        function Qc = get.Qc(obj)
            hc = 10.45 - obj.wv + 10*sqrt(obj.wv);
            Qc = hc * obj.Aw * (obj.Ts_C - obj.Ta_C) * obj.dt;
        end
        
        function Qe = get.Qe(obj)
            L = refEQ.L(obj.Ts_C);
            Qe = L * obj.qm * 1000;
        end
    end
end
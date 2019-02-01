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
        Qan
        Qsn
        Qbr
        Qc
        Qe
    end
    
    methods
        function obj = heat(dt)
            %HEAT Construct an instance of this class
            %   Detailed explanation goes here
            obj.dt = dt;
        end
        
        function Q = water(obj, Ts_C, Ta_C, Pv, qm, Qsn, wv, Aw)
            if ~isempty(Ts_C) && ~isempty(Ta_C) && ~isempty(Pv) && ...
                    ~isempty(qm) && ~isempty(Qsn) && ~isempty(wv) && ...
                    ~isempty(Aw)
                obj.Ta_C = Ta_C;
                obj.Ts_C = Ts_C;
                obj.Pv = Pv;
                obj.Aw = Aw;
                obj.qm = qm;
                obj.wv = wv;
            end
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
    
%     methods (Static)
%         function Qan = Qan(Area, Ta, Pv, dt)
%             Qan = sigma * Area * (Ta + 273.15).^4 * (A +  0.031 * sqrt(Pv/0.133322)) * (1 - 0.03) * dt;
%         end
%         
%         function Qc = Qc(Area, Ts, Ta, wv, dt)
%             hc = 10.45 - wv + 10.45 * sqrt(wv);
%             Qc = hc * Area * (Ts - Ta) * dt;
%         end
%         
%         function Qe = Qe(Ts_C, qm)
%             L = refEQ.L(Ts_C);
%             Qe = L * qm * 1000;
%         end
%     end
end
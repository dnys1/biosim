classdef refEQ
    %REFEQ Reference Equations for computation
    %   All static equations
    
    properties (Constant)
        Rv = 0.463; % gas constant for water vapor kJ/(kg-K)
    end
    
    methods (Static)
        function Pvs = Pvs(Ta_C)
            %Pvs Return the saturation vapor pressure estimate
            %   for a certain air temperature
            %
            %   Ta_C = Air temperature (C)
            Pvs = 0.611 * exp(17.3*Ta_C./(Ta_C + 237.3));
        end
        
        function RHOw = RHOw(Tw_C)
            %RHOw Return the water density estimate
            %   for a certain water temperature
            %
            %   Tw_C = Water temperature (C)
            RHOw = 999.5308+6.32693e-2*Tw_C-8.523829e-3*Tw_C.^2+6.943248e-5*Tw_C.^3-3.821216e-7*Tw_C.^4;
        end
        
        function L = L(Tw_C)
            %L Returns the latent heat of evaporation for water
            %   of a certain temperature (kJ/kg)
            %
            %   Tw_C = Water temperature (C)
            L = (-0.000000000000006178*Tw_C.^6+0.000000000008001366*Tw_C.^5-0.000000004930049019*Tw_C.^4+0.000000544072948524*Tw_C.^3-0.000044511189869516*Tw_C.^2-0.56229260095658*Tw_C+1093.20841244749)*2.326;
        end
        
        function Td_K = Td_K(Pv)
            %Td_K Returns the dew point temperature given the vapor pressure
            %
            %   Pv = Water vapor pressure (kPa)
            Td_K = (1/273.15 - refEQ.Rv/2.495e3 * log(Pv/0.611))^-1;
        end
        
        function Td_C = Td_C(Pv)
            %Td_C Returns the dew point temperature given the vapor pressure
            %
            %   Pv = Water vapor pressure (kPa)
            Td_C = refEQ.Td_K(Pv) - 273.15;
        end
        
        function gdd = GDD(trackTa)
           %GDD Calculate growing degree days(GDD)
           % https://en.wikipedia.org/wiki/Growing_degree-day
           Tbase = 10;
           T = trackTa - Tbase;
           T = T(T>0);
           gdd = max(trapz(T) / 24, 0);
        end
    end
end


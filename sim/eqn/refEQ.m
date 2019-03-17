classdef refEQ
    %REFEQ Reference Equations for computation
    %   All static equations
    
    properties (Constant)
        Ra = 0.288;     % gas constant for air kJ/(kg-K)
        Rv = 0.463;     % gas constant for water vapor kJ/(kg-K)
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
        
        function mu = mu(Tw_C)
            %mu Returns the dynamic viscosity of water over the range
            %   0 < T < 100 degC
            Tw_K = Tw_C + 273.15;
            mu = 10^-3 * exp(-3.7188+578.919/(Tw_K-137.546));
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
        
        function GDD = GDD(crop, trackTa)
           %GDD Calculate growing degree days(GDD)
           % https://en.wikipedia.org/wiki/Growing_degree-day
           T = min(trackTa, crop.Tcutoff) - crop.Tbase;
           T = T(T>0);
           GDD = max(trapz(T) / 24, 0);
        end
        
        function kAir = kAir(Ta_C)
            Ta = Ta_C + 273.15;
            kAir = 1.5207e-11*Ta^3-4.8574e-8*Ta^2+1.0184e-4*Ta-3.9333e-4;
        end
        
        function ETo = ETo(model)
            %ETo Returns the reference evapotranspiration per hour
            %
            %   Ta_C = air temperature (C)
            %   Rs = "solar" radiation (kJ)
            %   wv = wind velocity (m/s)
            %   P = Total air pressure (kPa)
            %   Pv = Vapor pressure (kPa)
            %
            %   Returns:
            %   ETo = reference evapotranspiration (mm/hr)
            
            % Steps from:
            % http://edis.ifas.ufl.edu/pdffiles/AE/AE45900.pdf

            sigma = 2.0413e-10; % Stefan-Boltzmann constant (MJ m^-2 K^-4 h^-1)
            Ta_C = model.Ta_C;
            Rs = model.Rs;
            wv = model.wv;
            P = model.P;
            Pv = model.Pv;

            % Calculate slop of saturation vapor pressure curve
            Delta = 4098 * (0.6108 * exp(17.27 * Ta_C / (Ta_C + 237.3)))/(Ta_C + 237.3).^2;

            % Calculate the psychrometric constant
            lamda = 2.45; % latent heat of vaporization (MJ/kg)
            Cp = 1.013e-3; % specific heat at constant pressure (MJ/kg)
            epsilon = 0.622; % ratio of molecular weight water vapor / dry air
            gamma = Cp * P / (epsilon * lamda);

            % Calculate Delta Term (DT)
            DT = Delta / (Delta + gamma * (1 + 0.34 * wv));

            % Calculate Psi Term (PT)
            PT = gamma / (Delta + gamma * (1 + 0.34 * wv));

            % Calcualte Temperature Term (TT)
            % 24 is correction so it's per hour instead of per day
            TT = ((900/24) / (Ta_C + 273.15)) * wv;

            % Calculate saturation vapor pressure
            Pvs = refEQ.Pvs(Ta_C);

            % Calculate short-wave & long-wave radiation terms
            Rns = Rs;
            Rnl = sigma * (Ta_C + 273.15)^4 * (0.34 - 0.14 * sqrt(Pv));
            Rn = Rns - Rnl;
            Rng = 0.408 * Rn;

            % Radiation Term (ET_rad)
            ET_rad = DT * Rng;

            % Wind Term (ET_wind)
            ET_wind = PT * TT * (Pvs - Pv);

            ETo = max(ET_rad + ET_wind, 0); % mm/(m^2-h)
        end
        
        function randn = rand(x1,x2)
            randn = x1 + rand * (x2-x1);
        end
    end
end


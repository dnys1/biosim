clc; clearvars;

hours = 1:1:24;

% Space dimensions
A = 10*10; % 100 m^2 room
V = 10*10*10; % 10 x 10 x 10 m room (m3)

% Lighting requirements for living space
% https://www.engineeringtoolbox.com/light-level-rooms-d_708.html
lumensperm2 = 1000;
% Bulbs used: http://www.alconlighting.com/specsheets/alcon/13124%20-%20Sinch.pdf
lumensperbulb = 1000;
heatperbulb = 12; % (W)
% No. of bulbs needed
no_bulbs = ceil(lumensperm2 * A / lumensperbulb);
heatperm2 = no_bulbs * heatperbulb / A;

Ta_C = 20; % Air temp = 20 deg C
Tw_C = Ta_C; % Water temp = air temp
To_C = 24; % Temp outside = temp underground
P = 101.325; % Pressure held constant
Pa = P; % kPa
Pv = 0; % kPa
Ra = 0.288; % gas constant for air kJ/(kg-K)
Rv = 0.463; % gas constant for water vapor kJ/(kg-K)
RHOa = Pa / ((Ta_C+273.15) * Ra);
RHOv = Pv / ((Ta_C+273.15) * Rv);
RHO = RHOa + RHOv;
Mv = RHOv * V; % No water vapor at time = 0
Ma = RHOa * V;
M = Mv + Ma;
Aw = 30; % Area of pond = 30 m^2
hw = 1; % depth of pond = 1 m
Vw = Aw * hw; % Volume of pond
epsilon = Ra/Rv; % epsilon in gas law equations
RHmax = 0.5; % maximum humidity allowed
ca = 1.000; % air specific heat (kJ/kg-deg C)
cw = 4.186; % water specific heat (kJ/kg-deg C)
cv = 1.996; % water vapor specific heat (kJ/kg-deg C)

kWall = 1.0; % thermal conductivity of building wall material (W/m-K)
sWall = 0.7; % wall thickness (m)

trackTa = zeros(size(hours));
trackTw = zeros(size(hours));
trackTd = zeros(size(hours));
trackGDD = zeros(size(hours));

for i = hours
   % Calculate water density based off temp
   % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909168/
   % Valid for 5 < T_C < 40
   if Tw_C < 5 || Tw_C > 40
       disp("ERROR: density formula not valid")
   else
       RHOw = 999.5308+6.32693e-2*Tw_C-8.523829e-3*Tw_C^2+6.943248e-5*Tw_C^3-3.821216e-7*Tw_C^4;
   end
   Mw = RHOw * Vw;
   
   %% Calculate the amount of heat exchange
   % Calculate the heat absorbed by water due to lighting
   % Valid for "daylight" hours (5AM - 8PM, for example)
   if mod(i,24) > 5 && mod(i,24) < 20
       Qsn = heatperm2 * Aw; % (W/m^2) * m^2
   else
       Qsn = 0;
   end
   
   % Update temperature of water
   % Calculate heat of evaporation
   % https://www.engineeringtoolbox.com/water-properties-d_1573.html
   % Valid for 0 < T_C < 350
   % H = round((-0.000000000000006178*Tw_C*Tw_C*Tw_C*Tw_C*Tw_C*Tw_C+0.000000000008001366*Tw_C*Tw_C*Tw_C*Tw_C*Tw_C-0.000000004930049019*Tw_C*Tw_C*Tw_C*Tw_C+0.000000544072948524*Tw_C*Tw_C*Tw_C-0.000044511189869516*Tw_C*Tw_C-0.56229260095658*Tw_C+1093.20841244749)*2.3260*10)/10;
   % Tw_C = -H*qm/(Mw*cw) + Tw_C; % negative because cooling

   dMv = 0; % Quantity of water evaporated
   for k = 1:1:3600 % seconds in an hour
       %% Calculate the amount of water released into air
       % Estimate saturation vapor pressure + density
       Pvs = 611 * exp(17.3*Ta_C/(Ta_C + 237.3)) / 1000; % Pa -> kPa
       RHOvs = Pvs / (Rv * (Tw_C + 273.15));

       % Wind velocity
       wv = 0;
       % Evaporation coefficient (kg/m2-s)
       ae = (25 + 19*wv) / 3600;
       % Hypothetical saturated vapor mass per total mass (kg/kg)
       % i.e. saturation humidity ratio
       x_s = RHOvs * V / Ma;
       % Actual (current) vapor mass per total mass (kg/kg)
       % i.e. humidity ratio
       x = Mv / Ma;
       % Amount of water vapor evaporated (kg/s)
       dMv = Aw * (x_s - x) * ae;
       
       % Solve nonlinear equation for equilibrium temp. of water
       fun = @(T) heat(T, Ta_C, Aw, dMv, Pv, wv, Qsn);
       [Tw_C, Q] = fzero(fun, Tw_C);
       
       % Calculate heat exchange due to walls of BioSim
       qWall = (kWall / sWall) * 10*10*4 * (To_C - Ta_C) / 1000;
       
       % Update temperature of air
       Ta_C = (qWall-(Q-Qsn)*dMv)/(Mv*cv + Ma*ca) + Ta_C;
   
       % Update mass of air
       % Assumption: Complete mixing of dry air and water vapor
       Mv = Mv + dMv;
       M = Ma + Mv;

       % Update density of air; V = constant
       RHOv = Mv / V;
       RHO = RHOa + RHOv;

       % Update partial pressures
       Pa = RHOa * Ra * (Ta_C + 273.15);
       Pv = RHOv * Rv * (Ta_C + 273.15);
       P = Pa + Pv;
       
       % Update mass of water
       Mw = Mw - dMv;

       % Update water body specs
       RHOw = 999.5308+6.32693e-2*Tw_C-8.523829e-3*Tw_C^2+6.943248e-5*Tw_C^3-3.821216e-7*Tw_C^4;
       Vw = Mw / RHOw;
       hw = Vw / Aw;
       
       % Dehumidify to achieve RHmax
       % https://www.sylvane.com/quest-hi-e-dry-195-dehumidifier.html
       RH = Pv / Pvs;
       if RH > RHmax
           disp('Dehumidifying');
           disp(RH);
           Pv0 = Pvs * RHmax;
           RHOv0 = Pv0 / (Rv * (Ta_C + 273.15));
           Mv0 = RHOv0 * V;
           
           % Calculate quantity of water removed
           Md = Mv - Mv0;
           Mw = Mw + Md;
           Mv = Mv0;
           Pv = Pv0;
           RHOv = RHOv0;
           
           % Calculate heat removed & update air temp
           H = round((-0.000000000000006178*Ta_C*Ta_C*Ta_C*Ta_C*Ta_C*Ta_C+0.000000000008001366*Ta_C*Ta_C*Ta_C*Ta_C*Ta_C-0.000000004930049019*Ta_C*Ta_C*Ta_C*Ta_C+0.000000544072948524*Ta_C*Ta_C*Ta_C-0.000044511189869516*Ta_C*Ta_C-0.56229260095658*Ta_C+1093.20841244749)*2.3260*10)/10;
           Ta_C = -H*Md/(Mv*cv + Ma*ca) + Ta_C;
           Ta_K = Ta_C + 273.15;
           
           % Update global values
           P = Pa + Pv;
           RHO = RHOa + RHOv;
           M = Ma + Mv;
       end
   end

   % Calculate dewpoint
   % Estimate saturation vapor pressure + density
   Pvs = 611 * exp(17.3*Ta_C/(Ta_C + 237.3)) / 1000; % Pa -> kPa
   Td_K = (1/273.15 - Rv/2.495e3 * log(Pv/0.611))^-1;
   Td_C = Td_K - 273.15;

   if Ta_C < Td_C
       disp("temp below dewpoint")
   end
   
   trackTa(i) = Ta_C;
   trackTw(i) = Tw_C;
   trackTd(i) = Td_C;
   trackGDD(i) = GDD(trackTa);
end

plot(hours, trackTa, hours, trackTw, hours, trackTd);
hold on
legend('T_{air}', 'T_{water}', 'T_{dewpoint}')
hold off
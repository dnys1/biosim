function Q = heat(Tw_C, Ta_C, Aw, qm, Pv, wv, Qsn)

sigma = 5.670367e-8; % Stefan-Boltzmann constant (W m^2 K^-4)

% Calculate heat from lamps (Jsn)
% Jsn = lamps;

% Calculate atmospheric long-wave radiation (Jan)
A = 0.6; % Atmospheric attenuation constant (between 0.5-0.7)
Qan = sigma * Aw * (Ta_C + 273.15).^4 * (A + 0.031 * sqrt(Pv/0.133322)) * (1 - 0.03);

% Calculate outgoing water long-wave radiation
Qbr = 0.97 * sigma * Aw * (Tw_C + 273.15).^4;

% Calculate outgoing heat due to conduction/convection
hc = 10.45 - wv + 10*sqrt(wv); % Convective heat transfer coefficient (W m^-2 K^-1)
Qc = hc * Aw * (Tw_C - Ta_C);

% Calculate outgoing heat due to evaporation
% L = (kJ/kg)
L = (-0.000000000000006178*Tw_C.^6+0.000000000008001366*Tw_C.^5-0.000000004930049019*Tw_C.^4+0.000000544072948524*Tw_C.^3-0.000044511189869516*Tw_C.^2-0.56229260095658*Tw_C+1093.20841244749)*2326;
Qe = L * qm;

Q = Qan + Qsn - (Qbr + Qc + Qe);
Q = Q / 1000;
end
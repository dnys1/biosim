function ETo = penman(Ta_C, Rs, wv, P, Pv)

% Steps from:
% http://edis.ifas.ufl.edu/pdffiles/AE/AE45900.pdf

sigma = 2.0413e-10; % Stefan-Boltzmann constant (MJ m^-2 K^-4 h^-1)

% Convert Rs in kJ to megajoules (MJ)
Rs = Rs / 1000;

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
Pvs = 0.6108 * exp(17.27*Ta_C / (Ta_C + 237.3));

% Calculate short-wave & long-wave radiation terms
Rns = Rs;
Rnl = sigma * (Ta_C + 273.15)^4 * (0.34 - 0.14 * sqrt(Pv));
Rn = Rns - Rnl;
Rng = 0.408 * Rn;

% Radiation Term (ET_rad)
ET_rad = DT * Rng;

% Wind Term (ET_wind)
ET_wind = PT * TT * (Pvs - Pv);

ETo = ET_rad + ET_wind; % mm/(m^2-h)
end
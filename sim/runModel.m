addpath blocks;
addpath data;
addpath eqn;

clc; clearvars;

% BioSim dimensions
L = 10;
W = 10;
H = 10;

% Initial pressure and temperatures
P = 101.325;
Ta = 20;
To = 20;
RH = 0.5;

% Time step
dt = time.MINUTE;

% Length of simulation
hours = 1:1:168*2;

% Build and run the model
m = model(L, W, H, P, Ta, To, RH, dt);
out = m.run(hours);

% Plot the output of air and water temperatures
plot(hours, out.Ta, hours, out.Tw);
hold on
legend('Ta', 'Tw');
hold off

% Plot the output of crop simulation
figure;
c=m.cropBlocks(1);
plot(hours,c.trackLeafB/1000,hours,c.trackStorB/1000,hours,c.trackStemB/1000,hours,c.trackRootB/1000);
hold on
legend('LeafB','StorB','StemB','RootB');
hold off
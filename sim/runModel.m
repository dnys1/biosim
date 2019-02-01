addpath blocks;
addpath data;
addpath eqn;

clc; clearvars;

L = 10;
W = 10;
H = 10;
P = 101.325;
Ta = 20;
To = 20;
RH = 0.5;

hours = 1:1:168;
m = model(L, W, H, P, Ta, To, RH);
out = m.run(hours);

plot(hours, out.Ta, hours, out.Tw);
hold on
legend('Ta', 'Tw');
hold off
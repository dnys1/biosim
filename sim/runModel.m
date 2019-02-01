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
dt = time.MINUTE;

hours = 1:1:24;
m = model(L, W, H, P, Ta, To, RH, dt);
out = m.run(hours);

plot(hours, out.Ta, hours, out.Tw);
hold on
legend('Ta', 'Tw');
hold off
figure;
plot(hours, out.ETo);
hold on
legend('ETo');
hold off
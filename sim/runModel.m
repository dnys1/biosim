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

hours = 1:1:168*15;
m = model(L, W, H, P, Ta, To, RH, dt);
out = m.run(hours);

plot(hours, out.Ta, hours, out.Tw);
hold on
legend('Ta', 'Tw');
hold off
figure;
c=m.cropBlocks(1);
plot(hours,c.trackLeafB/1000,hours,c.trackStorB/1000,hours,c.trackStemB/1000,hours,c.trackRootB/1000);
hold on
legend('LeafB','StorB','StemB','RootB');
hold off
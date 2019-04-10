addpath blocks;
addpath data;
addpath eqn;

clc; clearvars;

% Number of people in BioSim
no_people = 8;

% Initial pressure and temperatures
P = 101.325;
Ta = 20;
To = 20;
RH = 0.5;

% Time step
dt = time.HOUR;

% Length of simulation
hours = 1:1:173;

% Build and run the model
m = model(no_people, P, Ta, To, RH, dt);
out = m.run(hours);

% Plot the output of air and water temperatures
figure;
plot(hours, out.Ta, hours, out.Tw);
hold on
xlabel('Hour');
ylabel('Temperature (C)');
legend('Ta', 'Tw');
hold off

% Plot the output of crop simulation
figure;
trackLeafB=m.trackLeafB/1000;
trackStorB=m.trackStorB/1000;
trackStemB=m.trackStemB/1000;
trackRootB=m.trackRootB/1000;
plot(hours,trackLeafB,hours,trackStorB,hours,trackStemB,hours,trackRootB);
hold on
xlabel('Hour');
ylabel('Biomass (kg)');
legend('LeafB','StorB','StemB','RootB');
hold off
% figure;
% % DVS=c.trackDVS(c.trackDVS >= 0);
% % index_DVS0=length(hours)-length(DVS)+1;
% % trackLeafB=c.trackLeafB(index_DVS0:end)/1000;
% % trackStorB=c.trackStorB(index_DVS0:end)/1000;
% % trackStemB=c.trackStemB(index_DVS0:end)/1000;
% % trackRootB=c.trackRootB(index_DVS0:end)/1000;
% DVS=c.trackDVS;
% plot(DVS,trackLeafB,DVS,trackStorB,DVS,trackStemB,DVS,trackRootB);
% hold on
% xlabel('Development Stage');
% ylabel('Biomass (kg)');
% legend('LeafB','StorB','StemB','RootB');
% hold off

figure;
trackVcrop=m.trackVcrop;
trackVwetland=m.trackVwetland;
trackVvapor=m.trackVvapor;
trackVstore=m.trackVstore;
trackVsed=m.trackVsed;
trackVreject=m.trackVreject;
plot(hours, trackVstore, hours, trackVvapor, hours, trackVcrop, hours, trackVsed, ...
    hours, trackVwetland, hours, trackVreject);
hold on
legend('Storage','Air','Crop','Sed. Basin','Wetland','Reject');
xlabel('Hour');
ylabel('Volume of water (m^3)');
hold off
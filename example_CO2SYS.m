%% Example CO2SYSv3 function run
%
% This example function uses CO2SYSv3 to calculate carbonate system
% parameters using an input of total alkalinity and dissolved inorganic
% carbon, and compares aragonite saturation states calculated at laboratory
% conditions to aragonite saturation states calculated at in situ
% conditions.

% This is data from Station 156 of the P16N cruise in 2015 (44' N, 155' W)
PRS    = [3.4 40.6 69.5 109.1 174.6 250.0 332.7 432.4 550.1 701.0 850.9 1000.2 1167.2 1366.8 1599.9 1933.4 2366.6 2866.0 3366.6 3865.3 4367.8 4866.4 5282.4 5727.4];
TMP    = [11.805 10.0835 9.0949 8.7028 8.1385 7.2186 6.0406 4.8193 4.2995 3.9305 3.4850 3.1329 2.8238 2.5597 2.3071 2.0166 1.7961 1.6101 1.5137 1.4757 1.4998 1.5522 1.6036 1.6633];
SAL    = [33.095 33.0915 33.1315 33.7722 33.9339 33.942 33.9295 33.9477 34.053 34.1861 34.2706 34.3405 34.4022 34.4614 34.5131 34.5741 34.6161 34.6497 34.6675 34.6809 34.6850 34.6858 34.6859 34.6873];
SIL    = [4.21 7.44 8.51 15.06 25.43 38.92 53.79 72.66 89.96 106.57 121.23 132.96 142.73 151.81 159.73 165.97 167.73 164.21 159.80 157.56 160.39 162.83 163.42 163.52];
PO4    = [0.61 0.74 0.84 1.05 1.38 1.76 2.16 2.57 2.82 2.97 3.06 3.12 3.11 3.11 3.08 3.00 2.87 2.73 2.62 2.55 2.52 2.52 2.52 2.51];
DIC    = [2021.9 2028.9 2047.0 2096.7 2135.5 2179.2 2223.2 2269.0 2306.4 2338.4 2360.0 2375.0 2384.2 2393.9 2401.2 2387.7 2382.0 2362.0 2348.2 2338.6 2337.5 2337.1 2337.2 2338.8];
TA     = [2216.6 2215.2 2216.7 2250.9 2261.1 2278.2 2280.7 2301.5 2320.7 2336.8 2355.3 2369.7 2381.0 2393.9 2405.2 2413.0 2420.8 2422.2 2423.4 2425.4 2429.8 2431.6 2431.0 2432.4];

% Define constants for carbonate system calculations
SCALE  = 1; % Total pH scale
K1K2   = 10; % Leuker et al (2000) K1K2
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB

% Performs carbonate system calculations
C = CO2SYS(TA,DIC,1,2,SAL,20,TMP,0,PRS,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);

% Plot Omega(Ar)
OmA_in  = C(:,18);
OmA_out = C(:,35);
figure; clf; hold on
scatter(OmA_in, PRS,70,'Filled','Marker','o');
scatter(OmA_out,PRS,70,'Filled','Marker','^');
set(gcf,'Position',[200,200,600,600]);
set(gca,'Fontsize',15);
set(gca,'XAxisLocation', 'top')
xlabel('\Omega_{Ar}','Fontsize',18);
ylabel('Depth (dbar)');
legend('\Omega_{Ar,in}','\Omega_{Ar,out}','Location','southeast')
set(gca,'Ydir','reverse'); % reverse y-axis direction
hold off

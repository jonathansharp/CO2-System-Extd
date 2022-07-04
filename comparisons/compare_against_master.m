% Compares your latest commit with CO2SYS v3.2.0

%% Run CO2SYS v3.2.0
disp('Running CO2SYS v3.2.0')
system('git checkout v3.2.0'); % Checkout the v3.2.0 tag

%% Set up input conditions
PARvalues = [2250 2100 8.1 400 405];
PARTYPEs = 1:5;
pHSCALEIN_opts = 1:4;
K1K2CONSTANTS_opts = 1:15;
KSO4CONSTANTS_opts = 1:4;
KFCONSTANT_opts = 1;
SALvalue = 33.1;
[P1, P2, P1type, P2type, sal, pHscales, K1K2, ~, KSO4, KF, ...
    BSal, ~, ~] = CO2SYSigen(PARvalues, PARTYPEs, SALvalue, pHSCALEIN_opts, ...
    K1K2CONSTANTS_opts, KSO4CONSTANTS_opts, KFCONSTANT_opts);
tempin = 24;
tempout = 12;
presin = 1;
presout = 1647;
si = 10;
phos = 1;

%% Determine whether to calculate each input row or not
% xrow = 1 + 210; % just do one row, or...
xrow = 1:numel(P1); % ... do all rows (do this for saving output file)
P1 = P1(xrow);
P2 = P2(xrow);
P1type = P1type(xrow);
P2type = P2type(xrow);
sal = sal(xrow);
pHscales = pHscales(xrow);
K1K2 = K1K2(xrow);

tic
[DATA_v3, HEADERS_v3] = ...
    CO2SYS(P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
toc



% Circumvent caching of CO2SYS function
rehash path



%% Run current commit CO2SYS
disp('Running CO2SYS current commit/head')
system('git switch -'); % Switch back to your commit

%% Set up input conditions
PARvalues = [2250 2100 8.1 400 405];
PARTYPEs = 1:5;
pHSCALEIN_opts = 1:4;
K1K2CONSTANTS_opts = 1:15;
KSO4CONSTANTS_opts = 1:4;
KFCONSTANT_opts = 1;
SALvalue = 33.1;
[P1, P2, P1type, P2type, sal, pHscales, K1K2, ~, KSO4, KF, ...
    BSal, ~, ~] = CO2SYSigen(PARvalues, PARTYPEs, SALvalue, pHSCALEIN_opts, ...
    K1K2CONSTANTS_opts, KSO4CONSTANTS_opts, KFCONSTANT_opts);
tempin = 24;
tempout = 12;
presin = 1;
presout = 1647;
si = 10;
phos = 1;

%% Determine whether to calculate each input row or not
% xrow = 1 + 210; % just do one row, or...
xrow = 1:numel(P1); % ... do all rows (do this for saving output file)
P1 = P1(xrow);
P2 = P2(xrow);
P1type = P1type(xrow);
P2type = P2type(xrow);
sal = sal(xrow);
pHscales = pHscales(xrow);
K1K2 = K1K2(xrow);

tic
[DATA, HEADERS] = ...
    CO2SYS(P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
toc

fprintf("\n\nRelative change vs v3.2.0:\n"); ...
fprintf("%20s %20s %20s %20s %20s\n", "Variable", "Mean rel. change", "Min rel. change", "Max rel. change", "# of samples"); ...
for V = 1:length(HEADERS_v3)
    x = DATA_v3(:,V);
    y = DATA(:,V);
    relerr = abs(y - x) ./ abs(x);
    ix = x ~= -999; % Only compare non fill-in values
    maxrelerr = max(relerr(idx));
    minrelerr = min(relerr(idx));
    meanrelerr = mean(relerr(idx));
    fprintf("%20s %20.2g %20.2g %20.2g %20i\n", HEADERS_v3{V}, meanrelerr, minrelerr, maxrelerr, length(ix))
end
% Compares CO2SYS v3.2.0 with CO2SYS v2.0.5.
%
% CO2SYS v2.0.5 comes from https://github.com/jamesorr/CO2SYS-MATLAB,
%  but you must first rename the function to CO2SYSv2_0_5 (both inside the
%  file and in the file name).
%
% CO2SYS v3.2.0 comes from
% https://github.com/jonathansharp/CO2-System-Extd/v2_0_5_compatible
%
% CO2SYSigen comes from
% https://github.com/jonathansharp/CO2-System-Extd/comparisons/CO2SYSigen.m,
%
% compare_versions from
% https://github.com/jonathansharp/CO2-System-Extd/comparisons/compare_versions.m
%
% Corrctions for KSO4, KF, and BSal inputs and column
% headers from J. Sharp, 10 June 2020
%
% Differences are expected in outputs between the two versions due to minor
% differences in the way [CO2] is determined and a slight error in
% CO2SYSv2.0.5 in calculating Revelle factor

%% Set up input conditions
PARvalues = [2250 2100 8.1 400 405];
PARTYPEs = 1:5;
pHSCALEIN_opts = 1:4;
K1K2CONSTANTS_opts = 1:15;
KSO4CONSTANTS_opts = 1:4;
KFCONSTANT_opts = 1;
SALvalue = 33.1;
[P1, P2, P1type, P2type, sal, pHscales, K1K2, KSO4_only, KSO4, KF, ...
    BSal, U1, U2] = CO2SYSigen(PARvalues, PARTYPEs, SALvalue, pHSCALEIN_opts, ...
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
KSO4_only = KSO4_only(xrow);

%% Run CO2SYSv2.0.5
disp('Running CO2SYS v2.0.5...')
tic
[DATA_v2, HEADERS_v2] = ...
    CO2SYSv2_0_5(P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, pHscales, K1K2, KSO4_only);
toc

%% Run CO2SYSv3
disp('Running CO2SYS v3...')
tic
[DATA_v3, HEADERS_v3] = ...
    CO2SYS_adjusted_to_v2_0_5(P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
toc

%% Put results in tables
clear co2s_v2
for V = 1:numel(HEADERS_v2)
    co2s_v2.(HEADERS_v2{V}) = DATA_v2(:, V);
end
co2s_v2 = struct2table(co2s_v2);
clear co2s_v3
for V = 1:numel(HEADERS_v3)
    co2s_v3.(HEADERS_v3{V}) = DATA_v3(:, V);
end
co2s_v3 = struct2table(co2s_v3);

%% Calculate differences
clear co2s_diff
H=1;
for V = 1:numel(HEADERS_v3)
    if H < numel(HEADERS_v2)
    if isequal(HEADERS_v3{V},HEADERS_v2{H})
        co2s_diff.(HEADERS_v2{H}) = ...
           abs((co2s_v3.(HEADERS_v3{V}) - co2s_v2.(HEADERS_v2{H})) ./ co2s_v2.(HEADERS_v2{H})).*100;
    elseif isequal(HEADERS_v2{H},'KSO4CONSTANTS') && isequal(HEADERS_v3{V},'KSO4CONSTANT')
        co2s_diff.(HEADERS_v2{H}) = ...
           abs((co2s_v3.(HEADERS_v3{V}) - co2s_v2.(HEADERS_v2{H})) ./ co2s_v2.(HEADERS_v2{H})).*100;
    else
        H = H-1;
    end
    end
    H = H+1;
end % for V
co2s_diff = struct2table(co2s_diff);

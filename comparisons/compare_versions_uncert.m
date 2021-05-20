% Compares CO2SYS v3 with CO2SYS v2.0.5.
%
% errors v2.0.5 comes from https://github.com/jamesorr/CO2SYS-MATLAB,
%  but you must first rename the function to errorsv2_0_5 (both inside the
%  file and in the file name), and you must make a change in line 485 of
%  CO2SYSv2.0.5: case {'KW'} should be changed to case {'BOR'}
%
% errors v3 comes from
% https://github.com/jonathansharp/CO2-System-Extd/v2_0_5_compatible
%
% CO2SYSigen comes from
% https://github.com/jonathansharp/CO2-System-Extd/comparisons/CO2SYSigen.m,
%
% compare_versions_uncert from
% https://github.com/jonathansharp/CO2-System-Extd/comparisons/compare_versions.m

%% Add tools to path (if you need to!)
% addpath('/home/matthew/github/PyCO2SYS/validate')

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

%% Run uncertainties using CO2SYSv2.0.5 **SLOWWWWWW**
% Define dissociation constant uncertainties
epK = [0.002, 0.0075, 0.015, 0.01, 0.01, 0.02, 0.02];
disp('Running errors.m (CO2SYS v2.0.5)...')
tic
ERR_v2 = nan(size(P1,1),20);
ERR_HEADERS_v2 = cell(size(P1,1),20);
UNITS_v2 = cell(size(P1,1),20);
for n = 1:size(P1,1)
    [err, head, units] = ...
        errorsv2_0_5(P1(n), P2(n), P1type(n), P2type(n), sal(n), tempin, tempout, presin, ...
        presout, si, phos, U1(n), U2(n), 0.01, 0.02, 0.1, 0.01, epK, 0.02, 0.1, ...
        pHscales(n), K1K2(n), KSO4_only(n));
    ERR_v2(n,:) = err;
    ERR_HEADERS_v2 = head;
    UNITS_v2 = units;
end
toc

%% Run uncertainties using CO2SYSv3
epK = [0.002, 0.0075, 0.015, 0.01, 0.01, 0.02, 0.02];
disp('Running errors.m (CO2SYS v3)...')
tic
[ERR_v3, ERR_HEADERS_v3, UNITS_v3] = ...
    errors_adjusted_to_v2_0_5(P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, U1, U2, 0.01, 0.02, 0.1, 0.01, 0, 0, epK, ...
    0.02, 0.1, pHscales, K1K2, KSO4, KF, BSal, 0);
toc

%% Put results in tables
clear errs_v2
ERR_HEADERS_v2 = strrep(ERR_HEADERS_v2,'(','_');
ERR_HEADERS_v2 = strrep(ERR_HEADERS_v2,')','_');
for V = 1:numel(ERR_HEADERS_v2)
    errs_v2.(ERR_HEADERS_v2{V}) = ERR_v2(:, V);
end
errs_v2 = struct2table(errs_v2);
clear errs_v3
ERR_HEADERS_v3 = strrep(ERR_HEADERS_v3,'(','_');
ERR_HEADERS_v3 = strrep(ERR_HEADERS_v3,')','_');
for V = 1:numel(ERR_HEADERS_v3)
    errs_v3.(ERR_HEADERS_v3{V}) = ERR_v3(:, V);
end
errs_v3 = struct2table(errs_v3);

%% Calculate differences in errors
clear errs_diff
H=1;
for V = 1:8
    errs_diff.(ERR_HEADERS_v3{V}) = abs((errs_v3.(ERR_HEADERS_v3{V}) - ...
        errs_v2.(ERR_HEADERS_v2{V})) ./ errs_v2.(ERR_HEADERS_v2{V})).*100;
end
for V = 9:17
    errs_diff.(ERR_HEADERS_v3{V+1}) = abs((errs_v3.(ERR_HEADERS_v3{V+1}) - ...
        errs_v2.(ERR_HEADERS_v2{V})) ./ errs_v2.(ERR_HEADERS_v2{V})).*100;
end
for V = 18:20
    errs_diff.(ERR_HEADERS_v3{V+2}) = abs((errs_v3.(ERR_HEADERS_v3{V+2}) - ...
        errs_v2.(ERR_HEADERS_v2{V})) ./ errs_v2.(ERR_HEADERS_v2{V})).*100;
end
errs_diff = struct2table(errs_diff);

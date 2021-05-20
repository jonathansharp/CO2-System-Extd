% Generate derivatives across set of input conditions using CO2SYSv3.2.0

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

%% Run derivatives using CO2SYSv3
disp('Running derivnum v3...')
tic
[DERIV_PAR1_v3, HEADERS_PAR1_v3] = ...
    derivnum('PAR1', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_PAR2_v3, HEADERS_PAR2_v3] = ...
    derivnum('PAR2', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_TEMP_v3, HEADERS_TEMP_v3] = ...
    derivnum('t', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_SAL_v3, HEADERS_SAL_v3] = ...
    derivnum('s', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_SIL_v3, HEADERS_SIL_v3] = ...
    derivnum('sil', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_PHOS_v3, HEADERS_PHOS_v3] = ...
    derivnum('phos', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_K0_v3, HEADERS_K0_v3] = ...
    derivnum('k0', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_K1_v3, HEADERS_K1_v3] = ...
    derivnum('k1', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_K2_v3, HEADERS_K2_v3] = ...
    derivnum('k2', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_KB_v3, HEADERS_KB_v3] = ...
    derivnum('kb', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_KW_v3, HEADERS_KW_v3] = ...
    derivnum('kw', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
[DERIV_BOR_v3, HEADERS_BOR_v3] = ...
    derivnum('bor', P1, P2, P1type, P2type, sal, tempin, tempout, presin, ...
    presout, si, phos, 0, 0, pHscales, K1K2, KSO4, KF, BSal);
toc

%% Put results in tables
clear par1_v3
HEADERS_PAR1_v3 = strrep(HEADERS_PAR1_v3,'<','');
HEADERS_PAR1_v3 = strrep(HEADERS_PAR1_v3,'>','');
HEADERS_PAR1_v3 = strrep(HEADERS_PAR1_v3,'/','_');
for V = 1:numel(HEADERS_PAR1_v3)
    par1_v3.(HEADERS_PAR1_v3{V}) = DERIV_PAR1_v3(:, V);
end
par1_v3 = struct2table(par1_v3);
clear par2_v3
HEADERS_PAR2_v3 = strrep(HEADERS_PAR2_v3,'<','');
HEADERS_PAR2_v3 = strrep(HEADERS_PAR2_v3,'>','');
HEADERS_PAR2_v3 = strrep(HEADERS_PAR2_v3,'/','_');
for V = 1:numel(HEADERS_PAR2_v3)
    par2_v3.(HEADERS_PAR2_v3{V}) = DERIV_PAR2_v3(:, V);
end
par2_v3 = struct2table(par2_v3);
clear temp_v3
HEADERS_TEMP_v3 = strrep(HEADERS_TEMP_v3,'<','');
HEADERS_TEMP_v3 = strrep(HEADERS_TEMP_v3,'>','');
HEADERS_TEMP_v3 = strrep(HEADERS_TEMP_v3,'/','_');
for V = 1:numel(HEADERS_TEMP_v3)
    temp_v3.(HEADERS_TEMP_v3{V}) = DERIV_TEMP_v3(:, V);
end
temp_v3 = struct2table(temp_v3);
clear sal_v3
HEADERS_SAL_v3 = strrep(HEADERS_SAL_v3,'<','');
HEADERS_SAL_v3 = strrep(HEADERS_SAL_v3,'>','');
HEADERS_SAL_v3 = strrep(HEADERS_SAL_v3,'/','_');
for V = 1:numel(HEADERS_SAL_v3)
    sal_v3.(HEADERS_SAL_v3{V}) = DERIV_SAL_v3(:, V);
end
sal_v3 = struct2table(sal_v3);
clear sil_v3
HEADERS_SIL_v3 = strrep(HEADERS_SIL_v3,'<','');
HEADERS_SIL_v3 = strrep(HEADERS_SIL_v3,'>','');
HEADERS_SIL_v3 = strrep(HEADERS_SIL_v3,'/','_');
for V = 1:numel(HEADERS_SIL_v3)
    sil_v3.(HEADERS_SIL_v3{V}) = DERIV_SIL_v3(:, V);
end
sil_v3 = struct2table(sil_v3);
clear phos_v3
HEADERS_PHOS_v3 = strrep(HEADERS_PHOS_v3,'<','');
HEADERS_PHOS_v3 = strrep(HEADERS_PHOS_v3,'>','');
HEADERS_PHOS_v3 = strrep(HEADERS_PHOS_v3,'/','_');
for V = 1:numel(HEADERS_PHOS_v3)
    phos_v3.(HEADERS_PHOS_v3{V}) = DERIV_PHOS_v3(:, V);
end
phos_v3 = struct2table(phos_v3);
clear k0_v3
HEADERS_K0_v3 = strrep(HEADERS_K0_v3,'<','');
HEADERS_K0_v3 = strrep(HEADERS_K0_v3,'>','');
HEADERS_K0_v3 = strrep(HEADERS_K0_v3,'/','_');
for V = 1:numel(HEADERS_K0_v3)
    k0_v3.(HEADERS_K0_v3{V}) = DERIV_K0_v3(:, V);
end
k0_v3 = struct2table(k0_v3);
clear k1_v3
HEADERS_K1_v3 = strrep(HEADERS_K1_v3,'<','');
HEADERS_K1_v3 = strrep(HEADERS_K1_v3,'>','');
HEADERS_K1_v3 = strrep(HEADERS_K1_v3,'/','_');
for V = 1:numel(HEADERS_K1_v3)
    k1_v3.(HEADERS_K1_v3{V}) = DERIV_K1_v3(:, V);
end
k1_v3 = struct2table(k1_v3);
clear k2_v3
HEADERS_K2_v3 = strrep(HEADERS_K2_v3,'<','');
HEADERS_K2_v3 = strrep(HEADERS_K2_v3,'>','');
HEADERS_K2_v3 = strrep(HEADERS_K2_v3,'/','_');
for V = 1:numel(HEADERS_K2_v3)
    k2_v3.(HEADERS_K2_v3{V}) = DERIV_K2_v3(:, V);
end
k2_v3 = struct2table(k2_v3);
clear kb_v3
HEADERS_KB_v3 = strrep(HEADERS_KB_v3,'<','');
HEADERS_KB_v3 = strrep(HEADERS_KB_v3,'>','');
HEADERS_KB_v3 = strrep(HEADERS_KB_v3,'/','_');
for V = 1:numel(HEADERS_KB_v3)
    kb_v3.(HEADERS_KB_v3{V}) = DERIV_KB_v3(:, V);
end
kb_v3 = struct2table(kb_v3);
clear kw_v3
HEADERS_KW_v3 = strrep(HEADERS_KW_v3,'<','');
HEADERS_KW_v3 = strrep(HEADERS_KW_v3,'>','');
HEADERS_KW_v3 = strrep(HEADERS_KW_v3,'/','_');
for V = 1:numel(HEADERS_KW_v3)
    kw_v3.(HEADERS_KW_v3{V}) = DERIV_KW_v3(:, V);
end
kw_v3 = struct2table(kw_v3);
clear bor_v3
HEADERS_BOR_v3 = strrep(HEADERS_BOR_v3,'<','');
HEADERS_BOR_v3 = strrep(HEADERS_BOR_v3,'>','');
HEADERS_BOR_v3 = strrep(HEADERS_BOR_v3,'/','_');
for V = 1:numel(HEADERS_BOR_v3)
    bor_v3.(HEADERS_BOR_v3{V}) = DERIV_BOR_v3(:, V);
end
bor_v3 = struct2table(bor_v3);

writetable(par1_v3,'derivatives/par1_v3.csv')
writetable(par2_v3,'derivatives/par2_v3.csv')
writetable(temp_v3,'derivatives/temp_v3.csv')
writetable(sal_v3,'derivatives/sal_v3.csv')
writetable(sil_v3,'derivatives/sil_v3.csv')
writetable(phos_v3,'derivatives/phos_v3.csv')
writetable(k0_v3,'derivatives/k0_v3.csv')
writetable(k1_v3,'derivatives/k1_v3.csv')
writetable(k2_v3,'derivatives/k2_v3.csv')
writetable(kb_v3,'derivatives/kb_v3.csv')
writetable(kw_v3,'derivatives/kw_v3.csv')
writetable(bor_v3,'derivatives/bor_v3.csv')

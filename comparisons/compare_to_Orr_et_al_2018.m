%% Determine path
fpath = mfilename('fullpath');
path = fileparts(fpath) ;
 
%% Reproduce values in Table 2 of Orr et al., 2018
% Calculate derivatives
ders = derivnum_adjusted_to_v2_0_5('PAR1',2300,2000,1,2,35,18,NaN,0,NaN,0,0,0,0,1,10,1,1,2);
dv2(1,1:6) = ders([3 4 8 6 7 11]);
ders = derivnum_adjusted_to_v2_0_5('PAR2',2300,2000,1,2,35,18,NaN,0,NaN,0,0,0,0,1,10,1,1,2);
dv2(2,1:6) = ders([3 4 8 6 7 11]);
ders = derivnum_adjusted_to_v2_0_5('t',2300,2000,1,2,35,18,NaN,0,NaN,0,0,0,0,1,10,1,1,2);
dv2(3,1:6) = ders([3 4 8 6 7 11]);
ders = derivnum_adjusted_to_v2_0_5('s',2300,2000,1,2,35,18,NaN,0,NaN,0,0,0,0,1,10,1,1,2);
dv2(4,1:6) = ders([3 4 8 6 7 11]);
% Load Table 2
tab2 = readtable(strcat(path,'/data/orr2018-table2.csv'));
tab2 = table2array(tab2(:,3:end));
tab2 = tab2([1 7 13 17],:);
% Calculate differences
table2_differences = tab2-dv2

%% Reproduce values in Table 3 of Orr et al., 2018
% Calculate derivatives
ders = derivnum_adjusted_to_v2_0_5('phos',2300,2000,1,2,35,18,NaN,0,NaN,60,2,0,0,1,10,1,1,2);
dv3(1,1:6) = ders([3 4 8 6 7 11]);
ders = derivnum_adjusted_to_v2_0_5('sil',2300,2000,1,2,35,18,NaN,0,NaN,60,2,0,0,1,10,1,1,2);
dv3(2,1:6) = ders([3 4 8 6 7 11]);
% Load Table 3
tab3 = readtable(strcat(path,'/data/orr2018-table3.csv'));
tab3 = table2array(tab3(:,3:end));
tab3 = tab3([1 5],:);
% Calculate differences
table3_differences = tab3-dv3

%% Reproduce values in Table 4 of Orr et al., 2018
errs = errors_adjusted_to_v2_0_5(2300,2000,1,2,35,18,NaN,0,NaN,60,2,0,0,2,2,...
    0,0,4,0.1,0,0,[0 0 0 0 0 0 0],0,0,1,10,1,1,1);
err4(1,1:8) = errs([3 8 5 4 6 7 11 10]);
errs = errors_adjusted_to_v2_0_5(2300,2000,1,2,35,18,NaN,0,NaN,60,2,0,0,2,2,...
    0,0,4,0.1,0,0,[0.002 0.0075 0.015 0.01 0.01 0.02 0.02],0.02,0,1,10,1,1,1);
err4(2,1:8) = errs([3 8 5 4 6 7 11 10]);
% Load Table 3
tab4 = readtable(strcat(path,'/data/orr2018-table4.csv'));
tab4 = table2array(tab4(:,4:end));
tab4 = tab4([1 5],:);
% Calculate differences
table4_differences = tab4-err4

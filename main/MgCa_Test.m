% Mg Calcite Test

% read check values
TA=2300;
DIC=2150;
dt = readtable('check-values.xlsx');
x = table2array(dt(:,3));
T = repmat(25,36,1);%table2array(dt(:,4));
S = repmat(35,36,1);%table2array(dt(:,5));
P = table2array(dt(:,6));

index_2 = (1:4)';
index_1 = (5:8)';
index_3 = (9:12)';

% 
carb2 = CO2SYS(TA,DIC,1,2,S(index_2),T(index_2),NaN,P(index_2),NaN,...
    0,0,0,0,1,10,1,2,2,'MgContent',x(index_2),'MgTopt',2);
carb1 = CO2SYS(TA,DIC,1,2,S(index_1),T(index_1),NaN,P(index_1),NaN,...
    0,0,0,0,1,10,1,2,2,'MgContent',x(index_1),'MgTopt',1);
carb3 = CO2SYS(TA,DIC,1,2,S(index_3),T(index_3),NaN,P(index_3),NaN,...
    0,0,0,0,1,10,1,2,2,'MgContent',x(index_3),'MgTopt',3);

tt2 = readtable('testtable2.csv');
tt1 = readtable('testtable1.csv');
tt3 = readtable('testtable3.csv');

new_table = [tt2(1:4,:);tt1(1:4,:);tt3(1:4,:);tt2(5:8,:);tt1(5:8,:);...
    tt2(5:8,:);tt2(9:12,:);tt1(9:12,:);tt3(9:12,:)];
writetable(new_table,'testtable_all.csv');

% mgca_table = carb.OmegaMgCain1
% carb2 = CO2SYS(TA,DIC,1,2,35,25,NaN,1,NaN,...
%     0,0,0,0,1,10,1,2,2,'MgContent',x(index_2),'MgTopt',2);

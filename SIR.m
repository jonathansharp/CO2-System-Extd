function Ratio = SIR(CO2SYSDATA,CO2SYSHEADERS)
%**************************************************************************
%
% This function calculates the substrate inhibitor ratio (SIR) of free
% proton concentration to bicarbonate concentration from CO2SYS output.
%
% If HEADERS are provided, this function will work with output from any
% version of CO2SYS. If HEADERS are not provided, this function will work
% only with output from CO2SYS v3.0, found at:
% https://github.com/jonathansharp/CO2-System-Extd
%
% SIR Citation: Bach, L. T.: Reconsidering the role of carbonate ion
% concentration in calcification by marine organisms, Biogeosciences, 12,
% 4939?4951, https://doi.org/10.5194/bg-12-4939-2015, 2015.
%
% Function written by Jonathan Sharp
% University of South Florida
% June 11, 2020
%
%**************************************************************************

if exist('CO2SYSHEADERS','var')
    
   Hfreeidx = find(strcmp(CO2SYSHEADERS,'Hfreeout'));
   Hfree    = CO2SYSDATA(:,Hfreeidx);
   HCO3idx  = find(strcmp(CO2SYSHEADERS,'HCO3out'));
   HCO3     = CO2SYSDATA(:,HCO3idx);
   Ratio    = HCO3./(Hfree.*1e6);

else
    
   Hfree = CO2SYSDATA(:,32);
   HCO3  = CO2SYSDATA(:,23);
   Ratio = HCO3./(Hfree.*1e6);

end

end
%CO Oxidation Surface Reaction. Programmed by E. M.

%Grid Size
clear all;
clc;
dbstop if error

N=20; %number of sites on the grid= NxN
max_trials=100000;
CO2_time_step=200;
res=10; %plot resolution

yCO_V=linspace(0.3,0.6,res);
Coverage_fraction_O=zeros(1,res);
Coverage_fraction_CO=zeros(1,res);
CO2_Production_rate=zeros(1,res);
for i=1:res
    yCO=yCO_V(i);
    [Coverage_fraction_O(i), Coverage_fraction_CO(i), CO2_Production_rate(i)]=SteadyStateRun(N,max_trials,yCO,CO2_time_step);
end

plot(yCO_V,Coverage_fraction_O,'-r', yCO_V,Coverage_fraction_CO,'-.b', yCO_V,CO2_Production_rate,'-.or');
legend('O coverage fraction','CO coverage','CO2 coverage')


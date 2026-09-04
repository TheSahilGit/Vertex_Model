clear; clc; close all; 

it = 100; 
nrun = 2; 

[~, ~, ~, ~, ~, ~, ~, ~, all_end_data] = LoadData(it, nrun);

% Was reading these via readtable()+magic row numbers, which silently
% broke when if_PBC was inserted into para_Simulation.dat (every row
% after it shifted) -- see ReadPara1Params.m's header for why that
% pattern is fragile in the first place. Switched to the name-based
% reader so this can't silently desync again.
p1 = ReadPara1Params("../para_Simulation.dat");
dt = p1.dt;
eps0 = p1.Oscl_shearStrength;
w0 = p1.Oscl_freq_wo;


shearStress = all_end_data(:,2);
time = (1:length(shearStress))*dt; 

inputStrain = eps0 * sin(w0 * time);

figure('Position',[100 100 800 800])
plot(time, inputStrain, 'LineWidth',4, 'DisplayName','Strain');
hold on;
plot(time, shearStress, 'LineWidth',4, 'DisplayName','Stress');

axis square

xlabel("Time")

legend()
set(gca, 'FontSize', 36)

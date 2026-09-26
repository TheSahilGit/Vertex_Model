clear; clc; close all;
% Thin wrapper around compute_Circularity.m -- standalone circularity-vs-
% time plot, kept for interactive/ad-hoc use. For a multi-panel figure
% alongside other diagnostics, use PlotAnalysis('doCircularity', true, ...)
% instead, which calls the same compute_Circularity.m function.

nrun = 1;

p1 = ReadPara1Params("../para_Simulation.dat");
it_dump = p1.it_dump;
totT    = p1.totT;

itEnd = FindLatestAvailableIt(nrun, it_dump, totT);
itList = GetSnapshotItList(it_dump, 1, itEnd, []);

[time, circularity] = compute_Circularity(itList, nrun);

writematrix([time' circularity'], "circularity.dat");

figure("Position",[100 100 800 800])
plot(time, circularity, 'o', "LineWidth", 3, 'MarkerSize', 20);
hold on;
yline(0.95)

xlabel("Time");
ylabel("Circularity")
set(gca, "FontSize", 32)
axis square
legend()

clear; clc; close all;
% Thin wrapper around compute_MSD_cellID.m -- standalone MSD plot, kept
% for interactive/ad-hoc use. For a multi-panel figure alongside other
% diagnostics, use PlotAnalysis('doMSD', true, ...) instead, which calls
% the same compute_MSD_cellID.m function (tracked by persistent
% cell_identity, robust to T2 cell removal -- see that file's header).

nrun = 1;

p1 = ReadPara1Params("../para_Simulation.dat");
it_dump = p1.it_dump;
totT    = p1.totT;

itEnd = FindLatestAvailableIt(nrun, it_dump, totT);
itList = GetSnapshotItList(it_dump, it_dump, itEnd, 100000);

[time, MSD] = compute_MSD_cellID(itList, nrun);

writematrix([time' MSD'], "msd.dat");

figure("Position",[100 100 800 800])
loglog(time, MSD, 'o', "LineWidth", 3, 'MarkerSize', 20);
hold on;
loglog(time, 5e-4*time.^2, "LineWidth", 3, 'DisplayName', "t^2");
hold on;
loglog(time, 5e-4*time, "LineWidth", 3, "DisplayName", 't')

xlabel("Time");
ylabel("MSD")
set(gca, "FontSize", 32)
set(gcf, "Renderer", "Painter")
axis square
legend()

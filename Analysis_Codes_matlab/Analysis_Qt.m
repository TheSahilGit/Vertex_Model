clear; clc; close all;
% Thin wrapper around compute_Qt.m -- standalone Q(t) plot, kept for
% interactive/ad-hoc use. For a multi-panel figure alongside other
% diagnostics, use PlotAnalysis('doQt', true, ...) instead, which calls
% the same compute_Qt.m function.

nrun = 1;
ac = 1;

p1 = ReadPara1Params("../para1_in.dat");
it_dump = p1.it_dump;
totT    = p1.totT;

itEnd = FindLatestAvailableIt(nrun, it_dump, totT);
itList = GetSnapshotItList(it_dump, 1000, itEnd, 1000);

[time, Qt] = compute_Qt(itList, nrun, ac);

writematrix([time' Qt'], 'Qt.dat');

figure("Position",[100 100 800 800])
semilogx(time, Qt, 'o', "LineWidth", 3, 'MarkerSize', 20);
xlabel("Time");
ylabel("Q(t)")
set(gca, "FontSize", 32)
axis square
legend()

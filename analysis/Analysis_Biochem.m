clear; clc; close all;
% Biochemistry (Rho/ROCK/Myosin) diagnostics for runs with if_RhoROCK (or
% if_active_contractility, which also writes into the Myosin array) on
% -- see PARAMETERS.md. Reworked from an older ad hoc script
% (RhoROCK_Area_data.m, which had a stale/hardcoded it-range and no
% PBC-aware area calculation) into a general-purpose, standalone tool in
% the same style as Analysis_MSD_cellID.m/Analysis_COM.m: mean-field
% timeseries (Rho/ROCK/Myosin/area), a Myosin-vs-Rho phase portrait, and
% a Myosin/Rho cross-correlation vs lag.
%
% Thin wrapper around compute_BiochemMeans.m -- the same function
% PlotAnalysis.m's Biochem panel uses -- plus a standalone mean-area
% timeseries and two exploratory panels (phase portrait, cross-
% correlation) that aren't part of the multi-panel PlotAnalysis.m figure.

nrun = 1;

p1 = ReadPara1Params("../para_Simulation.dat");
it_dump = p1.it_dump;
totT    = p1.totT;

itEnd = FindLatestAvailableIt(nrun, it_dump, totT);
itList = GetSnapshotItList(it_dump, it_dump, itEnd, it_dump);

[time, meanRho, meanROCK, meanMyosin] = compute_BiochemMeans(itList, nrun);

% ---- Mean cell area vs time (PBC-aware unwrap -- same convention as
% ComputeCellColorData.m/Calculate_Total_Stress.m) ----
meanArea = zeros(size(itList));
for k = 1:numel(itList)
    [Lx, Ly, v, inn, num] = LoadData(itList(k), nrun);
    Nc = find(num ~= 0, 1, 'last');
    areaE = zeros(Nc,1);
    for ic = 1:Nc
        vx = v(inn(ic,1:num(ic)), 1);
        vy = v(inn(ic,1:num(ic)), 2);
        if numel(vx) > 1
            dx = vx(2:end) - vx(1); dx = dx - Lx.*round(dx./Lx); vx(2:end) = vx(1) + dx;
            dy = vy(2:end) - vy(1); dy = dy - Ly.*round(dy./Ly); vy(2:end) = vy(1) + dy;
        end
        areaE(ic) = polyarea(vx, vy);
    end
    meanArea(k) = mean(areaE);
end

%% ---- Timeseries + phase-portrait panel ----
figure('Position', [100 100 1600 500])

subplot(1,3,1)
plot(time, meanRho, 'LineWidth', 3, 'DisplayName', 'Rho'); hold on
plot(time, meanROCK, 'LineWidth', 3, 'DisplayName', 'ROCK');
plot(time, meanMyosin, '-o', 'LineWidth', 3, 'DisplayName', 'Myosin');
xlabel('Time'); ylabel('Mean field value'); legend('Location', 'best'); axis square
set(gca, 'FontSize', 20, 'LineWidth', 2)

subplot(1,3,2)
plot(time, meanArea, '--', 'LineWidth', 3);
xlabel('Time'); ylabel('Mean cell area'); axis square
set(gca, 'FontSize', 20, 'LineWidth', 2)

subplot(1,3,3)
plot(meanRho, meanMyosin, 'LineWidth', 3); hold on
scatter(meanRho(1), meanMyosin(1), 100, 'g', 'filled', 'DisplayName', 'start');
scatter(meanRho(end), meanMyosin(end), 100, 'r', 'filled', 'DisplayName', 'end');
xlabel('Mean Rho'); ylabel('Mean Myosin'); legend('Location', 'best'); axis square
title('Myosin vs Rho phase portrait')
set(gca, 'FontSize', 20, 'LineWidth', 2)

%% ---- Cross-correlation panel (Myosin against Rho, vs lag) ----
figure()
[xc, lags] = xcorr(meanMyosin - mean(meanMyosin), meanRho - mean(meanRho), 'coeff');
dt_frame = time(2) - time(1);
plot(lags*dt_frame, xc, 'LineWidth', 3)
xlabel('Lag (time)'); ylabel('Normalized cross-correlation (Myosin, Rho)')
axis square
set(gca, 'FontSize', 20, 'LineWidth', 2)

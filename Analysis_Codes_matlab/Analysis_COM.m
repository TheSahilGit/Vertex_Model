clear; clc;
% close all;

%% Load parameters
para2 = load("../para2_in.dat"); 
para1 = readtable("../para1_in.dat");

deltat = table2array(para1(8,1));

fname_etas = sprintf('../data/motility_store.dat');
fid = fopen(fname_etas);
dum4 = fread(fid,1,'float32');
etas = fread(fid,100000000,'float64');
fclose(fid);

%% Read domain info
Lx = para2(1);
Ly = para2(2);
numdim  = para2(3);
vdim1   = para2(4);
vdim2   = para2(5);
inndim1 = para2(6);
inndim2 = para2(7);

nrun = 2;

%% ---------------------------------------------------------
%% Compute COM trajectory + physical time
ct = 1;
for it = 5000:5000:10000000

    [Lx, Ly, v, inn, num, forces] = LoadData(it, nrun);    

    [cmX, cmY] = calculate_cellCentre(Lx, Ly, v, inn, num);

    % Global COM of tissue
    COMx(ct) = mean(cmX);
    COMy(ct) = mean(cmY);

    % Physical time
    time(ct) = deltat * it;

    ct = ct + 1
end

%% ---------------------------------------------------------
%% Plot COM magnitude vs physical time
r_com = sqrt(COMx.^2 + COMy.^2);

figure()
plot(time, r_com-r_com(1), 'o', 'MarkerSize', 32, 'LineWidth', 3);
hold on
yline(0, ':','LineWidth',3)
%ylim([11.5 12.5])
xlabel('Time');
ylabel('COM magnitude');
set(gca, 'FontSize', 32, 'LineWidth', 2)
set(gcf, 'Renderer', 'Painter')
axis square


%% ---------------------------------------------------------
%% Compute MSD of COM
N = length(COMx);
maxLag = floor(N/2);

MSD = zeros(maxLag,1);

for dt = 1:maxLag
    dx = COMx(1+dt:N) - COMx(1:N-dt);
    dy = COMy(1+dt:N) - COMy(1:N-dt);
    MSD(dt) = mean(dx.^2 + dy.^2);
end

%% Physical time-lag array
dt_phys = time(2) - time(1);
tau = (1:maxLag) * dt_phys;

%% ---------------------------------------------------------
%% Log-log MSD vs physical time lag
figure()
loglog(tau, MSD, 'o', 'MarkerSize', 24, 'LineWidth', 3);
hold on;
loglog(tau, 1e-4 * tau, '-', 'LineWidth', 3, 'DisplayName',"t");   % diffusive guide
hold on;
loglog(tau, 7e-7 * tau.^2 , '-', 'LineWidth', 3, 'DisplayName',"t^2");   % diffusive guide

legend('Location','northwest')
xlabel('Time lag');
ylabel('MSD of COM');
set(gca, 'FontSize', 32, 'LineWidth', 2)
axis square

%% ---------------------------------------------------------
%% Function: Cell-wise COM
function [cmX, cmY] = calculate_cellCentre(Lx, Ly, v, inn, num)

idx = find(inn(:,1) == 0, 1, 'first');

cmX = zeros(idx-1,1);
cmY = zeros(idx-1,1);

for i = 1:idx-1
    vx = v(inn(i,1:num(i)),1);
    vy = v(inn(i,1:num(i)),2);

    cmX(i) = mean(vx);
    cmY(i) = mean(vy);
end

end

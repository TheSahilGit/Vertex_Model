clear; clc; close all;

nrun = 2;
para1 = readtable("../para1_in.dat");
dt = table2array(para1(8,1));

% Initial data
it = 1000; 
[Lx, Ly, v, inn, num, ~, ~, cell_identity_init] = LoadData(it, nrun);
[cmX_init, cmY_init] = calculate_cellCentre(v,inn,num);

% Map identity → initial position
initPos = containers.Map;

for i = 1:Lx*Ly
    key = cell_identity_init(i);
    initPos(key) = [cmX_init(i), cmY_init(i)];
end

it_times = (1000:1000:10000000);
MSD = zeros(length(it_times),1);
times = it_times * dt; 
tidx = 0;

for it = it_times
    tidx = tidx + 1

    [Lx, Ly, v, inn, num, ~, ~, cell_identity] = LoadData(it, nrun);
    [cmX, cmY] = calculate_cellCentre(v,inn,num);

    msd_sum = 0;
    count   = 0;

    for i = 1:length(cell_identity)
        key = cell_identity(i);

        if isKey(initPos, key)
            r0 = initPos(key);
            dx = cmX(i) - r0(1);
            dy = cmY(i) - r0(2);

            msd_sum = msd_sum + (dx^2 + dy^2);
            count = count + 1;
        end
    end

    MSD(tidx) = msd_sum / count;
end

writematrix([times' MSD],"msd.dat");

%%

figure("Position",[100 100 800 800])

loglog(times, MSD,'o', "LineWidth",3, 'MarkerSize',20);
hold on; 
loglog(times, 5e-4*times.^2, "LineWidth",3, 'DisplayName',"t^2");
hold on; 
loglog(times, 5e-4*times, "LineWidth",3, "DisplayName",'t')


xlabel("Time");
ylabel("MSD")
set(gca, "FontSize", 32)
set(gcf, "Renderer", "Painter")

axis square
legend()


function [cmX, cmY] = calculate_cellCentre(v,inn,num)
idx = find(inn(:,1) == 0, 1, 'first');
for i = 1:idx-1
    vx=v(inn(i,1:num(i)),1);
    vy=v(inn(i,1:num(i)),2);

    cmX(i) = mean(vx); cmY(i) = mean(vy);
end
end
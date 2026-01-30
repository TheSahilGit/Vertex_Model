clear; clc ;
%close all;


para2 = load("../para2_in.dat");
%etas_in = load("../motility_in.dat");

para1 = readtable("../para1_in.dat");

deltat = table2array(para1(8,1));


fname_etas = sprintf('../data/motility_store.dat');
fid = fopen(fname_etas);
dum4 = fread(fid,1,'float32');
etas = fread(fid,100000000,'float64');


%%

Lx = para2(1);
Ly = para2(2);
numdim  = para2(3);
vdim1 = para2(4);
vdim2 = para2(5);
inndim1 = para2(6);
inndim2 = para2(7);


nrun = 2;

ct = 1;
for it = 10000:10000:10000000
    [Lx, Ly, v,inn,num, forces] = LoadData(it, nrun);    

    [globalMeanX, globalMeanY] = calculate_globalCentre(v, inn, num) ;

    r_com(ct) = sqrt(globalMeanX.^2 + globalMeanY.^2);
    time(ct) = deltat * it;


    ct = ct + 1
end

%%
figure()
plot(time, r_com-r_com(1), 'o', 'MarkerSize',32, 'LineWidth',3);
xlabel("Time");
ylabel("Global COM")
set(gca, "FontSize", 32)
set(gcf, "Renderer", "Painter")
axis square
legend()

%%



function [globalMeanX, globalMeanY] = calculate_globalCentre(v, inn, num)

    % Find number of cells (same logic you used)
    idx = find(inn(:,1) == 0, 1, 'first') - 1;

    % Preallocate (upper bound)
    totalVerts = sum(num(1:idx));
    all_vx = zeros(totalVerts,1);
    all_vy = zeros(totalVerts,1);

    ct = 1;
    for i = 1:idx
        ids = inn(i,1:num(i));
        nv  = num(i);

        all_vx(ct:ct+nv-1) = v(ids,1);
        all_vy(ct:ct+nv-1) = v(ids,2);

        ct = ct + nv;
    end

    % Global mean over ALL vertices
    globalMeanX = mean(all_vx);
    globalMeanY = mean(all_vy);

end

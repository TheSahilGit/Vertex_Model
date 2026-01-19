clear; clc ;
%close all;


para2 = load("../para2_in.dat");
%etas_in = load("../motility_in.dat");



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
for it = 1000:1000:620000
    [Lx, Ly, v,inn,num, forces] = LoadData(it, nrun);    


    [cmX, cmY] = calculate_cellCentre(Lx,Ly,v,inn,num);

    globalMeanX = mean(cmX);
    globalmeanY = mean(cmY);

    r_com(ct) = sqrt(globalMeanX.^2 + globalmeanY.^2);


    ct = ct + 1
end

%%
figure()
plot(r_com, 'o', 'MarkerSize',32, 'LineWidth',3);
xlabel("Time");
ylabel("COM")
set(gca, "FontSize", 32)
set(gcf, "Renderer", "Painter")
axis square
legend()


function [cmX, cmY] = calculate_cellCentre(Lx,Ly,v,inn,num)
idx = find(inn(:,1) == 0, 1, 'first');
for i = 1:idx-1

    vx=v(inn(i,1:num(i)),1);
    vy=v(inn(i,1:num(i)),2);

    cmX(i) = mean(vx); cmY(i) = mean(vy);
end
end


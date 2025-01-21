clear; clc; close all;


nrun = 1;


it = 1;
[Lx, Ly, v,inn,num, ~] = LoadData(it, nrun);

for i = 1:Lx*Ly

    vx = v(inn(i, 1:num(i)), 1); % x-coordinates of neighbors
    vy = v(inn(i, 1:num(i)), 2); % y-coordinates of neighbors

    plot(polyshape(vx,vy), FaceAlpha=0.001, EdgeColor='k'); hold on;


    cmx(i) = mean(vx);
    cmy(i) = mean(vy);
    
end

gmeanX1 = mean(cmx);
gmeanY1 = mean(cmy);



it = 2;
[Lx, Ly, v, inn, num, ~] = LoadData(it, nrun);

for i = 1:Lx*Ly

    vx = v(inn(i, 1:num(i)), 1); % x-coordinates of neighbors
    vy = v(inn(i, 1:num(i)), 2); % y-coordinates of neighbors
    
    plot(polyshape(vx,vy), FaceAlpha=0.001, EdgeColor='r'); hold on;

    cmx(i) = mean(vx);
    cmy(i) = mean(vy);
    
end


gmeanX2 = mean(cmx);
gmeanY2 = mean(cmy);




it = 500000;
[Lx, Ly, v, inn, num, ~] = LoadData(it, nrun);

for i = 1:Lx*Ly

    vx = v(inn(i, 1:num(i)), 1); % x-coordinates of neighbors
    vy = v(inn(i, 1:num(i)), 2); % y-coordinates of neighbors
    
    plot(polyshape(vx,vy), FaceAlpha=0.001, EdgeColor='b'); hold on;

    cmx(i) = mean(vx);
    cmy(i) = mean(vy);
    
end

axis square


axis square




gmeanX3 = mean(cmx);
gmeanY3 = mean(cmy);



shift1 = (gmeanX1 - gmeanX2)^2 + (gmeanY1 - gmeanY2)^2
shift2 = (gmeanX1 - gmeanX3)^2 + (gmeanY1 - gmeanY3)^2



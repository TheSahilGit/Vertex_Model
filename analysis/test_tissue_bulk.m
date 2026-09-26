clear; clc; close all;



nrun = 2;
it = 100;

[Lx, Ly, v, inn, num, ~, biochem, ~] = LoadData(it, nrun);

R = 4;
[inside_cells, outside_cells, COM, cell_centers] = Cells_Within_Radius(Lx, Ly, v, inn, num, R);

Nc = find(num ~= 0, 1, 'last');

colordata = rand(1,Nc);
colorbar_string = "\beta";

norm_flag = 'data';
norm_range = [];


TisuePlot(Lx,Ly,v,inn,num, colordata, colorbar_string, norm_flag, norm_range)
hold on
scatter(cell_centers(:,1), cell_centers(:,2))

for ii = 1:length(inside_cells)
    i  = inside_cells(ii);
    vx = v(inn(i,1:num(i)),1);
    vy = v(inn(i,1:num(i)),2);

    pl = polyshape(vx,vy);



    plot(pl, ...
        FaceColor = 'k', ...
        FaceAlpha = 0.5, ...
        LineWidth = 1.5);
end


for ii = 1:length(outside_cells)
    i  = outside_cells(ii);
    vx = v(inn(i,1:num(i)),1);
    vy = v(inn(i,1:num(i)),2);

    pl = polyshape(vx,vy);



    plot(pl, ...
        FaceColor = 'w', ...
        FaceAlpha = 0.5, ...
        LineWidth = 1.5);
end


%%

nrun = 2;
it = 100;

[Lx, Ly, v, inn, num, ~, biochem, ~] = LoadData(it, nrun);


figure()
TisuePlot(Lx,Ly,v,inn,num, colordata, colorbar_string, norm_flag, norm_range)
hold on


%%----

inside = load('../inside.dat');

for ic = 1:inside(1)
    i = inside(ic+1);
    vx = v(inn(i,1:num(i)),1);
    vy = v(inn(i,1:num(i)),2);

    pl = polyshape(vx,vy);

    plot(pl, ...
        FaceColor = 'k', ...
        FaceAlpha = 0.5, ...
        LineWidth = 1.5);
    hold on
end







%%




function [inside_cells, outside_cells, COM, cell_centers] = Cells_Within_Radius(Lx, Ly, v, inn, num, R)

Nc = find(num ~= 0, 1, 'last');

cell_centers = zeros(Nc,2);

% ---- Step 1: Compute cell centroids ----
for i = 1:Nc
    vx = v(inn(i,1:num(i)),1);
    vy = v(inn(i,1:num(i)),2);

    pl = polyshape(vx, vy);
    [cx, cy] = centroid(pl);

    cell_centers(i,:) = [cx, cy];
end

% ---- Step 2: Global COM ----
COM = mean(cell_centers, 1);

% ---- Step 3: Distance from COM ----
dx = cell_centers(:,1) - COM(1);
dy = cell_centers(:,2) - COM(2);

dist = sqrt(dx.^2 + dy.^2);

% ---- Step 4: Classification ----
inside_cells  = find(dist <= R);
outside_cells = find(dist > R);

end

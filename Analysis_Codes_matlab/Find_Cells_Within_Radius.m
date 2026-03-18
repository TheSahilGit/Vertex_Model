function [inside_cells, outside_cells, COM, cell_centers] = Find_Cells_Within_Radius(Lx, Ly, v, inn, num, R)

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
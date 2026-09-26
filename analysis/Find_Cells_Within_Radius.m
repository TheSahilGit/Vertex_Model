function [inside_cells, outside_cells, COM, cell_centers] = Find_Cells_Within_Radius(Lx, Ly, v, inn, num, R)

Nc = find(num ~= 0, 1, 'last');

cell_centers = zeros(Nc,2);

% ---- Step 1: Compute cell centroids ----
for i = 1:Nc
    vx = v(inn(i,1:num(i)),1);
    vy = v(inn(i,1:num(i)),2);

    % PBC-aware (log.txt): unwrap a wrap-straddling cell relative to its
    % own first vertex before building the polyshape -- otherwise a cell
    % that straddles the periodic wrap draws a self-intersecting/huge
    % polygon and centroid() returns garbage (same root cause as
    % TisuePlot.m's "wired" rendering bug; see ComputeCellColorData.m).
    if numel(vx) > 1
        dx = vx(2:end) - vx(1); dx = dx - Lx .* round(dx ./ Lx);
        vx(2:end) = vx(1) + dx;
        dy = vy(2:end) - vy(1); dy = dy - Ly .* round(dy ./ Ly);
        vy(2:end) = vy(1) + dy;
    end

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
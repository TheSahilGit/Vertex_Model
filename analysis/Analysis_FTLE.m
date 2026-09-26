clear; clc; close all;

nrun = 2;

para1 = readtable("../para1_in.dat");
dt = table2array(para1(8,1));

it0 =   10000;
it1 = 1000000;

plot_ftle_map      = 1;
plot_shape_map     = 0;
plot_stress_map    = 1;
plot_scatter       = 0;
plot_ftle_pdf      = 0;
plot_ftle_hist     = 0;
plot_ftle_vs_pbar  = 0;
plot_ftle_pressure_overlay = 1;

boundary_layers = 3;   % set to 2 or 4 (or any integer)

[Lx,Ly,v0,inn0,num0,~,~,cell_identity0] = LoadData(it0,nrun);
[~,~,v1,inn1,num1,~,~,cell_identity1] = LoadData(it1,nrun);

[cmX0,cmY0] = calculate_cellCentre(v0,inn0,num0);
[cmX1,cmY1] = calculate_cellCentre(v1,inn1,num1);

pos0 = [cmX0(:) cmY0(:)];
pos1 = [cmX1(:) cmY1(:)];

ncell0 = size(pos0,1);
ncell1 = size(pos1,1);

cell_identity0 = cell_identity0(1:ncell0);
cell_identity1 = cell_identity1(1:ncell1);

neigh = BuildNeighbors(inn0,num0,ncell0);

T = (it1-it0)*dt;

FTLE = ComputeFTLE( ...
            pos0,...
            pos1,...
            cell_identity0,...
            cell_identity1,...
            neigh,...
            T);

[Area,Perim,ShapeIndex] = ComputeShapeIndex(v0,inn0,num0);

pbar = ComputeNeighborhoodShapeIndex(ShapeIndex,neigh);

dp = ShapeIndex - pbar;

% [~,PressureIndividual,ShearStress_Individual] = ...
%     calculate_total_stress(Lx,Ly,v1,inn1,num1);

 [~, Pressure_Individual, ShearStress_Individual] = ...
     calculate_total_stress(Lx, Ly, v1, inn1, num1, boundary_layers);


PrintFTLEStats(FTLE);


%% --- Bulk mask: remove n boundary layers using cell index structure ---

%boundary_layers = 4;   % set to 2 or 4

% number of cells along each direction
Nx = round(Lx);   % number of columns
Ny = round(Ly);   % number of rows  (cells 1..Ny are column 1)

% column and row index for each cell (1-based)
col_idx = ceil((1:ncell0)' / Ny);
row_idx = mod((1:ncell0)' - 1, Ny) + 1;

bulk_mask = (col_idx > boundary_layers) & ...
            (col_idx <= Nx - boundary_layers) & ...
            (row_idx > boundary_layers) & ...
            (row_idx <= Ny - boundary_layers);

FTLE_bulk = FTLE;
FTLE_bulk(~bulk_mask) = NaN;

%%
if plot_ftle_map || plot_shape_map || plot_stress_map

    figure("Position",[200 200 1700 800])

    nplot = plot_ftle_map + plot_shape_map + plot_stress_map;
    pidx = 1;

    if plot_ftle_map

        subplot(1,nplot,pidx)

        TisuePlot(Lx,Ly,v0,inn0,num0,...
                  FTLE_bulk,...
                  "FTLE",...
                  'data',[]);

        title("FTLE")

        pidx = pidx + 1;

    end

    if plot_shape_map

        subplot(1,nplot,pidx)

        TisuePlot(Lx,Ly,v0,inn0,num0,...
                  ShapeIndex,...
                  " Shape Index",...
                  'data',[]);

        title(" Shape Index")

        pidx = pidx + 1;

    end

    if plot_stress_map

        subplot(1,nplot,pidx)
        % 
        % TisuePlot(Lx,Ly,v0,inn0,num0,...
        %           ShearStress_Individual,...
        %           "Shear Stress",...
        %           'data',[]);

        TisuePlot(Lx,Ly,v0,inn0,num0,...
            Pressure_Individual,...
            "Pressure",...
            'data',[]);


        title(" Stress")

    end

end

if plot_scatter

    PlotScatter(Pressure_Individual,...
                FTLE,...
                "Pressure",...
                "FTLE")

end

if plot_ftle_vs_pbar

    PlotScatter(pbar,...
                FTLE,...
                "Neighborhood Shape Index",...
                "FTLE")

end

if plot_ftle_hist

    PlotFTLEHistogram(FTLE,60);

end

if plot_ftle_pdf

    ratio_min = 1;
    ratio_max = it1/it0;

    npoints = 30;

    ratio_list = logspace(log10(ratio_min),...
                          log10(ratio_max),...
                          npoints);

    ratio_list = round(ratio_list);

    ratio_list = unique(ratio_list);

    dt_list = it0 * ratio_list - it0;

    dt_list = dt_list(dt_list > 0);

    PlotFTLE_PDF_TimeIntervals( ...
            it0,...
            dt_list,...
            nrun,...
            60)

end

%% -- Pressure colour map + top/bottom 20% FTLE contours ---

if plot_ftle_pressure_overlay

    figure('Position', [200 200 900 900])

    % --- background: pressure colour map (bulk only) ---
    Pressure_bulk = Pressure_Individual(:);
    if length(Pressure_bulk) < ncell0
        Pressure_bulk(end+1:ncell0) = NaN;
    end
    Pressure_bulk(~bulk_mask) = NaN;

    TisuePlot(Lx, Ly, v0, inn0, num0, ...
              Pressure_bulk, ...
              'Pressure', ...
              'data', []);

    hold on

    % --- interpolate FTLE onto regular grid ---
    Ngrid = 512;
    xg = linspace(0, Lx, Ngrid);
    yg = linspace(0, Ly, Ngrid);
    [Xg, Yg] = meshgrid(xg, yg);

    valid = bulk_mask & ~isnan(FTLE_bulk);

    F_interp = scatteredInterpolant( ...
                    cmX0(valid), cmY0(valid), FTLE_bulk(valid), ...
                    'natural', 'none');

    FTLE_grid = F_interp(Xg, Yg);

    % --- compute top/bottom 20% thresholds ---
    ftle_vals  = FTLE_bulk(valid);
    thresh_top = prctile(ftle_vals, 80);   % top 20%
    thresh_bot = prctile(ftle_vals, 20);   % bottom 20%

    % --- top 20% FTLE contour: green dashed ---
    contour(Xg, Yg, FTLE_grid, [thresh_top thresh_top], ...
            'LineWidth', 6, ...
            'LineColor', [0 0 0], ...   % green
            'LineStyle', '-')

    % --- bottom 20% FTLE contour: red dashed ---
    contour(Xg, Yg, FTLE_grid, [thresh_bot thresh_bot], ...
            'LineWidth', 6, ...
            'LineColor', [1 1 1], ...  % red
            'LineStyle', '-')

    % --- legend ---
    h1 = plot(nan, nan, '-', 'Color', [0.1 0.7 0.1], 'LineWidth', 6);
    h2 = plot(nan, nan, '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 6);
    legend([h1 h2], ...
           {'Top 20% FTLE', 'Bottom 20% FTLE'}, ...
           'Location', 'northeast', ...
           'FontSize', 16, ...
           'Box', 'off')

    title(sprintf('Pressure + FTLE contours  (\\Deltat = %d)', it1-it0))

    set(gca, 'FontSize', 22, 'LineWidth', 2)
    axis equal tight
    hold off

    set(gcf, 'Renderer', 'Painters')

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bulk_mask = ComputeBulkMask(cmX, cmY, Lx, Ly, neigh, ncell, n_layers)
% Erodes n_layers of cells from each edge of the periodic box.
%
% Strategy:
%   Layer 0 (boundary cells): any cell whose centre is within a margin
%   from any box edge. The margin is estimated as the mean nearest-neighbour
%   distance, so it adapts to cell size automatically.
%   Subsequent layers are grown inward by flood-fill on the neighbour graph.

% --- estimate typical cell diameter from mean nn distance ---
all_dist = [];
for i = 1:min(ncell, 500)          % sample up to 500 cells for speed
    nb = neigh{i};
    for k = 1:length(nb)
        j = nb(k);
        d = sqrt((cmX(i)-cmX(j))^2 + (cmY(i)-cmY(j))^2);
        all_dist(end+1) = d;
    end
end
cell_diam = mean(all_dist);         % ~1 cell diameter

% --- mark layer-0: cells touching the box boundary ---
margin = 0.5 * cell_diam;          % half a cell width from each wall

is_boundary = (cmX < margin)       | (cmX > Lx - margin) | ...
              (cmY < margin)       | (cmY > Ly - margin);

% --- flood-fill inward for n_layers ---
layer = zeros(ncell, 1);           % 0 = interior (not yet condemned)
layer(is_boundary) = 1;            % layer 1 = box-edge cells

for L = 2 : n_layers
    newly_added = false(ncell, 1);
    for i = 1:ncell
        if layer(i) ~= 0
            continue               % already condemned
        end
        nb = neigh{i};
        for k = 1:length(nb)
            j = nb(k);
            if layer(j) == L-1     % neighbour was added in previous layer
                newly_added(i) = true;
                break
            end
        end
    end
    layer(newly_added) = L;
end

bulk_mask = (layer == 0);          % true = bulk (never condemned)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FTLE = ComputeFTLE( ...
                    pos0,...
                    pos1,...
                    id0,...
                    id1,...
                    neigh,...
                    T)

ncell0 = size(pos0,1);
ncell1 = size(pos1,1);

FTLE = nan(ncell0,1);

idMap1 = containers.Map('KeyType','char',...
                        'ValueType','double');

for k = 1:ncell1
    idMap1(char(id1{k})) = k;
end

for i = 1:ncell0

    key_i = char(id0{i});

    if ~isKey(idMap1,key_i)
        continue
    end

    nb = neigh{i};

    if length(nb) < 3
        continue
    end

    ii1 = idMap1(key_i);

    xi0 = pos0(i,:);
    xi1 = pos1(ii1,:);

    R0 = [];
    R1 = [];

    for kk = 1:length(nb)

        j = nb(kk);

        if j > ncell0
            continue
        end

        key_j = char(id0{j});

        if ~isKey(idMap1,key_j)
            continue
        end

        jj1 = idMap1(key_j);

        r0 = pos0(j,:) - xi0;
        r1 = pos1(jj1,:) - xi1;

        R0(:,end+1) = r0';
        R1(:,end+1) = r1';

    end

    if size(R0,2) < 3
        continue
    end

    F = R1 * pinv(R0);

    C = F' * F;

    eigvals = eig(C);

    lam = max(real(eigvals));

    if lam <= 0
        continue
    end

    FTLE(i) = (1/T)*log(sqrt(lam));

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function neigh = BuildNeighbors(inn,num,ncell)

neigh = cell(ncell,1);

for i = 1:ncell

    vi = inn(i,1:num(i));

    for j = i+1:ncell

        vj = inn(j,1:num(j));

        ncommon = length(intersect(vi,vj));

        if ncommon >= 2

            neigh{i}(end+1) = j;
            neigh{j}(end+1) = i;

        end

    end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cmX,cmY] = calculate_cellCentre(v,inn,num)

idx = find(inn(:,1)==0,1,'first');

cmX = zeros(idx-1,1);
cmY = zeros(idx-1,1);

for i = 1:idx-1

    vx = v(inn(i,1:num(i)),1);
    vy = v(inn(i,1:num(i)),2);

    cmX(i) = mean(vx);
    cmY(i) = mean(vy);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Area, Perim, ShapeIndex] = ComputeShapeIndex(v,inn,num)

Nc = find(num~=0,1,'last');

Area       = nan(Nc,1);
Perim      = nan(Nc,1);
ShapeIndex = nan(Nc,1);

for i = 1:Nc

    verts = inn(i,1:num(i));

    vx = v(verts,1);
    vy = v(verts,2);

    Area(i) = polyarea(vx,vy);

    P = 0;

    nv = length(vx);

    for k = 1:nv

        kp = mod(k,nv) + 1;

        dx = vx(kp) - vx(k);
        dy = vy(kp) - vy(k);

        P = P + sqrt(dx^2 + dy^2);

    end

    Perim(i) = P;

    ShapeIndex(i) = P/sqrt(Area(i));

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pbar = ComputeNeighborhoodShapeIndex(ShapeIndex,neigh)

ncell = length(ShapeIndex);

pbar = nan(ncell,1);

for i = 1:ncell

    nb = neigh{i};

    vals = ShapeIndex(i);

    for k = 1:length(nb)

        j = nb(k);

        if ~isnan(ShapeIndex(j))
            vals(end+1) = ShapeIndex(j);
        end

    end

    pbar(i) = mean(vals);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotScatter(x,y,xlab,ylab)

figure('Position',[200 200 600 600])

scatter(x,y,120,'LineWidth',4)

xlabel(xlab)
ylabel(ylab)

set(gca,...
    'FontSize',30,...
    'FontName','Sans',...
    'LineWidth',4)

box on

set(gcf,'Renderer','Painters')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotFTLEHistogram(FTLE,nbins)

vals = FTLE(~isnan(FTLE));

figure('Position',[200 200 700 600])

histogram(vals,...
          nbins,...
          'Normalization','pdf')

xlabel('FTLE')
ylabel('PDF')

set(gca,...
    'FontSize',30,...
    'LineWidth',3)

box on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PrintFTLEStats(FTLE)

vals = FTLE(~isnan(FTLE));

fprintf('\n');

fprintf('Valid FTLE cells = %d / %d\n', ...
        sum(~isnan(FTLE)), ...
        length(FTLE));

fprintf('Mean FTLE = %.6e\n',mean(vals));
fprintf('Std  FTLE = %.6e\n',std(vals));
fprintf('Max  FTLE = %.6e\n',max(vals));
fprintf('Min  FTLE = %.6e\n',min(vals));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotFTLE_PDF_TimeIntervals( ...
                it0,...
                dt_list,...
                nrun,...
                nbins)

[Lx,Ly,v0,inn0,num0,~,~,cell_identity0] = ...
    LoadData(it0,nrun);

[cmX0,cmY0] = calculate_cellCentre(v0,inn0,num0);

pos0 = [cmX0(:) cmY0(:)];

ncell = length(cmX0);

cell_identity0 = cell_identity0(1:ncell);

neigh = BuildNeighbors(inn0,num0,ncell);

figure('Position',[100 100 900 700])

legend_entries = cell(length(dt_list),1);

for kk = 1:length(dt_list)

    dtime = dt_list(kk);

    it1 = it0 + dtime;

    [~,~,v1,inn1,num1,~,~,cell_identity1] = ...
        LoadData(it1,nrun);

    [cmX1,cmY1] = calculate_cellCentre(v1,inn1,num1);

    pos1 = [cmX1(:) cmY1(:)];

    ncell1 = length(cmX1);

    cell_identity1 = cell_identity1(1:ncell1);

    FTLE = ComputeFTLE( ...
                pos0,...
                pos1,...
                cell_identity0,...
                cell_identity1,...
                neigh,...
                dtime);

    FTLE = FTLE(~isnan(FTLE));

    [counts,edges] = histcounts( ...
                        FTLE,...
                        nbins,...
                        'Normalization','pdf');

    mean_FTLE(kk) = mean(FTLE);
    centers = 0.5*(edges(1:end-1)+edges(2:end));

    plot(centers,...
         counts,...
         'LineWidth',3)

    hold on

    legend_entries{kk} = sprintf('\\Deltat = %d',dtime);

end

xlabel('FTLE')
ylabel('PDF')

legend(legend_entries,...
       'Location','best')

set(gca,...
    'FontSize',32,...
    'LineWidth',2)

box on

title(sprintf('FTLE PDF from t_0 = %d',it0))

figure('Position',[500 500 900 700])

loglog(dt_list, mean_FTLE, 'o','LineWidth', 4,'MarkerSize',40, 'DisplayName','Data')
hold on;
loglog(dt_list, 0.08*dt_list.^(-1), 'DisplayName','\Delta T ^{-1}')

xlabel("\Delta T")
ylabel("<FTLE>")

legend()
set(gca,...
    'FontSize',32,...
    'LineWidth',2)
set(gcf, "Renderer","painters")

end

function [time, Qt] = compute_Qt(itList, nrun, ac, progressFcn)
% COMPUTE_QT  Self-overlap order parameter Q(t) vs. time.
%
%   [time, Qt] = compute_Qt(itList, nrun, ac, progressFcn)
%
% itList(1) is treated as the reference ("t=0") configuration; Qt(k) is
% the fraction of cells whose center-of-mass has moved less than `ac`
% (a cage/persistence length scale) from its position at itList(1).
% time(k) = itList(k)*dt, dt read from para_Simulation.dat.
% progressFcn : optional function handle progressFcn(k, n), called after
%          each of the n requested frames is processed (used by
%          PlotAnalysis.m to print live progress); omit or pass [] for none.
%
% Extracted from Analysis_Qt.m's CalculateQt; the only behavior change is
% looping over Nc = find(num~=0,1,'last') instead of a hardcoded Lx*Ly, so
% this stays correct with if_cell_division on.

if nargin < 3 || isempty(ac)
    ac = 1.0;
end
if nargin < 4
    progressFcn = [];
end

p1 = ReadPara1Params("../para_Simulation.dat");
time = itList * p1.dt;

Qt = zeros(size(itList));
vcmXIn = []; vcmYIn = [];
n = numel(itList);

for k = 1:n
    it = itList(k);
    [Lx, Ly, v, inn, num] = LoadData(it, nrun);

    Nc = find(num ~= 0, 1, 'last');
    vcmX = zeros(Nc,1); vcmY = zeros(Nc,1);
    for ic = 1:Nc
        vx = v(inn(ic, 1:num(ic)), 1);
        vy = v(inn(ic, 1:num(ic)), 2);
        % PBC-aware (log.txt): unwrap a wrap-straddling cell's vertices
        % relative to its own first vertex before averaging -- otherwise
        % the raw mean lands in the middle of the box, not at the cell's
        % true centroid (same technique as ComputeCellColorData.m).
        if numel(vx) > 1
            dxu = vx(2:end) - vx(1); dxu = dxu - Lx .* round(dxu ./ Lx);
            vx(2:end) = vx(1) + dxu;
            dyu = vy(2:end) - vy(1); dyu = dyu - Ly .* round(dyu ./ Ly);
            vy(2:end) = vy(1) + dyu;
        end
        vcmX(ic) = mean(vx);
        vcmY(ic) = mean(vy);
    end

    if k == 1
        vcmXIn = vcmX;
        vcmYIn = vcmY;
    end

    % Only compare against cells that existed at the reference time too
    % (division can grow Nc beyond the reference frame's cell count).
    NcRef = min(Nc, numel(vcmXIn));
    ddx = vcmX(1:NcRef) - vcmXIn(1:NcRef);
    ddy = vcmY(1:NcRef) - vcmYIn(1:NcRef);
    % PBC-aware: minimum-image the displacement itself too, so a cell
    % whose centroid legitimately drifted across the periodic wrap since
    % the reference frame doesn't read as a ~Lx spurious jump.
    ddx = ddx - Lx .* round(ddx ./ Lx);
    ddy = ddy - Ly .* round(ddy ./ Ly);
    dis = sqrt(ddx.^2 + ddy.^2);
    Qt(k) = sum(ac - dis > 0) / NcRef;

    if ~isempty(progressFcn); progressFcn(k, n); end
end

end

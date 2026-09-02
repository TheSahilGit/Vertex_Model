function [time, circularity] = compute_Circularity(itList, nrun, progressFcn)
% COMPUTE_CIRCULARITY  Tissue outer-boundary circularity vs. time.
%
%   [time, circularity] = compute_Circularity(itList, nrun, progressFcn)
%
% itList : vector of Fortran timesteps to evaluate (see GetSnapshotItList
%          -- each must have an actual dump file on disk).
% nrun   : 1 or 2 (LoadData convention).
% progressFcn : optional function handle progressFcn(k, n), called after
%          each of the n requested frames is processed (used by
%          PlotAnalysis.m to print live progress); omit or pass [] for none.
%
% circularity(k) = isoperimetric ratio 4*pi*Area/Perimeter^2 of the convex
% hull of the tissue's outer-boundary vertices at itList(k); 1 for a
% perfect circle, less than 1 otherwise. time(k) = itList(k)*dt, with dt
% read from para1_in.dat (not hardcoded) so the time axis stays correct
% regardless of which run this is pointed at.
%
% Extracted from Analysis_Circularity.m's core loop; the only behavior
% change is that FindBorderVertices below now scans all live cells
% (Nc = find(num~=0,1,'last')) instead of a hardcoded 1:Lx*Ly, so this
% stays correct for runs with if_cell_division on (Nc can exceed Lx*Ly).

if nargin < 3
    progressFcn = [];
end

p1 = ReadPara1Params("../para1_in.dat");
time = itList * p1.dt;

circularity = zeros(size(itList));
n = numel(itList);

for k = 1:n
    it = itList(k);
    [~, ~, v, inn, num] = LoadData(it, nrun);

    borderver = FindBorderVertices(inn, num);
    vx = v(borderver, 1);
    vy = v(borderver, 2);
    circularity(k) = calculateCircularityWithQhull(vx, vy);

    if ~isempty(progressFcn); progressFcn(k, n); end
end

end


function borderver = FindBorderVertices(inn, num)
% A vertex is on the tissue's outer boundary iff it is shared by fewer
% than 3 cells (an interior vertex is always shared by exactly 3).
Nc = find(num ~= 0, 1, 'last');
bordercount = zeros(1, max(inn(:)));

for ic = 1:Nc
    nn = num(ic);
    for jc = 1:nn
        bordercount(inn(ic, jc)) = bordercount(inn(ic, jc)) + 1;
    end
end

borderver = find(bordercount > 1 & bordercount < 3);
end


function circularity = calculateCircularityWithQhull(x, y)
K = convhull(x, y);
hullX = x(K);
hullY = y(K);

area = polyarea(hullX, hullY);

dx = diff([hullX; hullX(1)]);
dy = diff([hullY; hullY(1)]);
perimeter = sum(sqrt(dx.^2 + dy.^2));

circularity = (4 * pi * area) / (perimeter^2);
end

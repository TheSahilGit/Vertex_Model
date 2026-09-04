function [time, meanArea, meanPerimeter, meanShapeFactor] = compute_ShapeMetrics(itList, nrun, progressFcn)
% COMPUTE_SHAPEMETRICS  Mean per-cell area, perimeter, and shape factor
% (P/sqrt(A)) vs. time.
%
%   [time, meanArea, meanPerimeter, meanShapeFactor] = compute_ShapeMetrics(itList, nrun, progressFcn)
%
% Reuses ComputeCellColorData.m's 'Area'/'Perimeter' computation directly
% (the same PBC-aware, minimum-image-unwrapped per-cell logic the movie
% coloring options use) -- this is the tissue-wide MEAN of those same
% per-cell values, not a separate calculation.
%
% Shape factor is computed PER CELL first (perimeter(i)/sqrt(area(i))),
% then averaged -- NOT perimeter-of-means/sqrt(area-of-means), a different
% and less meaningful quantity. This matches the standard convention in
% the confluent-tissue jamming literature (Bi & Manning et al.), where the
% mean shape index p0 = <P_i/sqrt(A_i)> is the usual solid-fluid jamming
% order parameter (p0* ~ 3.81 is the commonly-cited transition value).
%
% itList : vector of Fortran timesteps to evaluate (see GetSnapshotItList
%          -- each must have an actual dump file on disk).
% nrun   : 1 or 2 (LoadData convention).
% progressFcn : optional function handle progressFcn(k, n), called after
%          each of the n requested frames is processed (used by
%          PlotAnalysis.m to print live progress); omit or pass [] for none.
%
% time(k) = itList(k)*dt, with dt read from para_Simulation.dat.

if nargin < 3
    progressFcn = [];
end

p1 = ReadPara1Params("../para_Simulation.dat");
time = itList * p1.dt;

meanArea = zeros(size(itList));
meanPerimeter = zeros(size(itList));
meanShapeFactor = zeros(size(itList));
n = numel(itList);

for k = 1:n
    it = itList(k);
    [Lx, Ly, v, inn, num] = LoadData(it, nrun);

    cellArea = ComputeCellColorData('Area', v, inn, num, [], [], [], Lx, Ly);
    cellPerimeter = ComputeCellColorData('Perimeter', v, inn, num, [], [], [], Lx, Ly);
    cellShapeFactor = cellPerimeter ./ sqrt(cellArea);

    meanArea(k) = mean(cellArea);
    meanPerimeter(k) = mean(cellPerimeter);
    meanShapeFactor(k) = mean(cellShapeFactor);

    if ~isempty(progressFcn); progressFcn(k, n); end
end

end

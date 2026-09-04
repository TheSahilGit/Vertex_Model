function [time, meanRho, meanROCK, meanMyosin] = compute_BiochemMeans(itList, nrun, progressFcn)
% COMPUTE_BIOCHEMMEANS  Mean per-cell Rho, ROCK, and Myosin vs. time (the
% RhoROCK mechanochemical fields -- if_RhoROCK in para_Simulation.dat,
% see PARAMETERS.md).
%
%   [time, meanRho, meanROCK, meanMyosin] = compute_BiochemMeans(itList, nrun, progressFcn)
%
% biochemdata (LoadData.m) is numdim x 3 = [Rho, ROCK, Myosin], one row per
% cell SLOT -- including unused headroom beyond the live cell count, which
% must be excluded from the mean (see Nc below), same convention as every
% other per-cell panel in this toolkit. Meaningful when if_RhoROCK (or
% if_active_contractility, which also writes into the Myosin array) was on
% for the run being analyzed; flat at/near 0 otherwise.
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

meanRho = zeros(size(itList));
meanROCK = zeros(size(itList));
meanMyosin = zeros(size(itList));
n = numel(itList);

for k = 1:n
    it = itList(k);
    [~, ~, ~, ~, num, ~, biochemdata] = LoadData(it, nrun);

    Nc = find(num ~= 0, 1, 'last');
    meanRho(k) = mean(biochemdata(1:Nc, 1));
    meanROCK(k) = mean(biochemdata(1:Nc, 2));
    meanMyosin(k) = mean(biochemdata(1:Nc, 3));

    if ~isempty(progressFcn); progressFcn(k, n); end
end

end

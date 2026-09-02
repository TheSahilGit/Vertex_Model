function [time, maxForce, meanForce] = compute_ForceStats(itList, nrun, progressFcn)
% COMPUTE_FORCESTATS  Max and mean per-vertex force magnitude vs. time,
% read from the force_*/nrun2_force_* dumps written every it_dump.
%
%   [time, maxForce, meanForce] = compute_ForceStats(itList, nrun, progressFcn)
%
% itList : vector of Fortran timesteps to evaluate (see GetSnapshotItList
%          -- each must have an actual dump file on disk).
% nrun   : 1 or 2 (LoadData convention).
% progressFcn : optional function handle progressFcn(k, n), called after
%          each of the n requested frames is processed (used by
%          PlotAnalysis.m to print live progress); omit or pass [] for none.
%
% Force magnitude here is sqrt(fxx^2 + fyy^2) -- the same deterministic
% elastic/motility force ComputeCellColorData.m's 'Force' coloring option
% uses (forces' columns 1-2; columns 3-8 are the random/ABP/polar terms,
% which vertexmain.f90's update rule scales by sqrt(dt) rather than dt, so
% they aren't on the same footing as fxx/fyy and aren't folded into this
% magnitude).
%
% Restricted to vertices actually referenced by a live cell (num~=0):
% v_dim2 reserves extra headroom for future divisions/growth
% (Generate_Initial_Mesh.f90/Main.m), and including those unused slots
% would dilute the mean (they'd count as exact zeros) and, worse, could
% silently corrupt the max if they ever hold stale/uninitialized data.
%
% This is a useful blow-up/instability diagnostic: a sudden spike in
% maxForce (while meanForce stays flat) usually means one vertex is
% getting an unphysically large force -- e.g. from a near-degenerate T1/T2
% edge or a division artifact -- worth checking before it destabilizes dt.
%
% time(k) = itList(k)*dt, with dt read from para_Simulation.dat (not
% hardcoded) so the time axis stays correct regardless of which run this
% is pointed at.

if nargin < 3
    progressFcn = [];
end

p1 = ReadPara1Params("../para_Simulation.dat");
time = itList * p1.dt;

maxForce = zeros(size(itList));
meanForce = zeros(size(itList));
n = numel(itList);

for k = 1:n
    it = itList(k);
    [~, ~, ~, inn, num, forces] = LoadData(it, nrun);

    Nc = find(num ~= 0, 1, 'last');
    usedVerts = unique(inn(1:Nc, :));
    usedVerts = usedVerts(usedVerts > 0);   % drop each cell row's zero padding

    Fmag = sqrt(forces(:,1).^2 + forces(:,2).^2);
    maxForce(k) = max(Fmag(usedVerts));
    meanForce(k) = mean(Fmag(usedVerts));

    if ~isempty(progressFcn); progressFcn(k, n); end
end

end

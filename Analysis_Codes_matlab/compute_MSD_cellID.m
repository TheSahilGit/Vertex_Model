function [time, MSD] = compute_MSD_cellID(itList, nrun, progressFcn)
% COMPUTE_MSD_CELLID  Mean-squared cell displacement vs. time, tracked by
% persistent cell identity (NOT raw cell index).
%
%   [time, MSD] = compute_MSD_cellID(itList, nrun, progressFcn)
%
% progressFcn : optional function handle progressFcn(k, n), called after
%          each of the n requested frames is processed (used by
%          PlotAnalysis.m to print live progress); omit or pass [] for none.
%
% Cell INDICES are not a stable label across a run: T2 extrusion events
% remove a cell and shift every higher-indexed cell down by one
% (T2_swap.f90), and cell division inserts a new cell at Nc+1. cell_identity
% (a persistent string like "cell_17", written every dump) is the only
% label safe to track across time. This function tracks every cell present
% at the reference time itList(1); a cell whose identity later stops
% appearing (removed by a T2 event) is simply excluded from the average
% from that point on (`count` shrinks accordingly) rather than causing an
% error or silently mixing in the wrong cell. Daughter cells created by
% division after itList(1) are not part of the reference cohort and are
% correctly never included (standard convention: MSD tracks a fixed
% initial population).
%
% Extracted from Analysis_MSD_cellID.m's core loop, generalized so the
% reference cohort is Nc = find(num~=0,1,'last') at itList(1) rather than a
% hardcoded Lx*Ly (matters if the reference time itself is after some
% divisions have already happened).

if nargin < 3
    progressFcn = [];
end

p1 = ReadPara1Params("../para_Simulation.dat");
time = itList * p1.dt;

% ---- Reference configuration ----
it0 = itList(1);
[~, ~, v0, inn0, num0, ~, ~, cell_identity0] = LoadData(it0, nrun);
[cmX0, cmY0] = calculate_cellCentre(v0, inn0, num0);

initPos = containers.Map;
Nc0 = find(num0 ~= 0, 1, 'last');
for i = 1:Nc0
    initPos(cell_identity0(i)) = [cmX0(i), cmY0(i)];
end

% ---- MSD at each requested time ----
MSD = zeros(size(itList));
n = numel(itList);

for k = 1:n
    it = itList(k);
    [~, ~, v, inn, num, ~, ~, cell_identity] = LoadData(it, nrun);
    [cmX, cmY] = calculate_cellCentre(v, inn, num);

    % BUGFIX (log.txt): cell_identity is always num_dim-long (LoadData.m
    % reads the full reserved capacity, not just the live cells), but
    % cmX/cmY (calculate_cellCentre, below) are only Nc-long -- looping to
    % length(cell_identity) indexes cmX/cmY out of bounds the moment
    % num_dim > Nc (i.e. essentially always, since num_dim = Nc +
    % headroom). The reference-cohort loop above already gets this right
    % (loops to Nc0, not length(cell_identity0)); this loop just needs the
    % same fix for the current frame's own Nc.
    Nc = find(num ~= 0, 1, 'last');
    msd_sum = 0;
    count = 0;
    for i = 1:Nc
        key = cell_identity(i);
        if isKey(initPos, key)
            r0 = initPos(key);
            dx = cmX(i) - r0(1);
            dy = cmY(i) - r0(2);
            msd_sum = msd_sum + (dx^2 + dy^2);
            count = count + 1;
        end
    end

    if count > 0
        MSD(k) = msd_sum / count;
    else
        MSD(k) = NaN;
        warning('compute_MSD_cellID:noSurvivors', ...
            'None of the reference cells (from it=%d) still exist at it=%d.', it0, it);
    end

    if ~isempty(progressFcn); progressFcn(k, n); end
end

end


function [cmX, cmY] = calculate_cellCentre(v, inn, num)
Nc = find(num ~= 0, 1, 'last');
cmX = zeros(Nc,1); cmY = zeros(Nc,1);
for i = 1:Nc
    vx = v(inn(i, 1:num(i)), 1);
    vy = v(inn(i, 1:num(i)), 2);
    cmX(i) = mean(vx);
    cmY(i) = mean(vy);
end
end

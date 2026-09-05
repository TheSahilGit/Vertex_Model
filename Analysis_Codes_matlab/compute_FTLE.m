function FTLE = compute_FTLE(it0, it1, nrun, v0, inn0, num0, cell_identity0)
% COMPUTE_FTLE  Spatial finite-time Lyapunov exponent (FTLE) field, one
% value per live cell at the reference time `it0`.
%
%   FTLE = compute_FTLE(it0, it1, nrun)
%   FTLE = compute_FTLE(it0, it1, nrun, v0, inn0, num0, cell_identity0)
%
% For a discrete, cell-based tissue (no continuous velocity field to
% differentiate), the standard way to get a local deformation gradient is
% to fit the best affine map between a cell's neighbors' relative
% positions at the two times (the same idea as Falk--Langer/D2min
% non-affine-deformation analysis): for cell $i$ with edge-sharing
% neighbors $j$,
%   r0_j = pos0(j) - pos0(i),   r1_j = pos1(j) - pos1(i),
% stacked as columns of R0, R1 (2 x n_neighbors); the least-squares local
% deformation gradient is F = R1 * pinv(R0), and
%   FTLE(i) = (1/T) * log(sqrt(lambda_max(F' * F))),
% the standard Cauchy-Green/FTLE formula, T = (it1-it0)*dt the physical
% time between the two snapshots (dt from para_Simulation.dat).
%
% Cells are matched between the two snapshots by their PERSISTENT
% cell_identity (not raw index, which drifts across T2/division events --
% same convention as compute_MSD_cellID.m/compute_Qt.m); a cell whose
% identity doesn't exist at it1 (extruded in between), or whose
% edge-sharing neighborhood is too small/degenerate for a well-posed fit,
% gets NaN.
%
% Positions are cell centroids, PBC-unwrapped per cell before averaging
% (a cell straddling the periodic wrap would otherwise get a centroid in
% the middle of the box -- same fix as ComputeCellColorData.m/
% Analysis_COM.m), and neighbor RELATIVE positions are themselves
% minimum-imaged (a neighbor across the periodic wrap is physically
% close, not ~box-width away).
%
% it0, it1 : reference and target Fortran timesteps (it1 can be either
%            side of it0; T's sign follows accordingly -- typically
%            it1 > it0, a forward-time FTLE).
% nrun     : 1 or 2 (LoadData convention).
% v0,inn0,num0,cell_identity0 : optional -- if the caller (e.g.
%            Movie_Code.m's per-frame loop) has already loaded it0's data
%            for its own rendering, pass it here to skip a redundant
%            LoadData(it0,nrun) call. If omitted, loaded internally.
%
% FTLE is returned indexed 1..Nc0 (it0's live cell count, i.e. the same
% indexing as inn0/num0), ready to hand straight to TisuePlot alongside
% v0/inn0/num0.

if nargin < 4
    [Lx, Ly, v0, inn0, num0, ~, ~, cell_identity0] = LoadData(it0, nrun);
else
    para2 = load('../para_MeshDims.dat');
    Lx = para2(1);
    Ly = para2(2);
end
[~, ~, v1, inn1, num1, ~, ~, cell_identity1] = LoadData(it1, nrun);

p1 = ReadPara1Params('../para_Simulation.dat');
T = (it1 - it0) * p1.dt;
if T == 0
    warning('compute_FTLE:zeroTimeSeparation', ...
        'it0 == it1 -- FTLE is undefined for zero time separation; returning all-NaN.');
    Nc0 = find(num0 ~= 0, 1, 'last');
    FTLE = nan(Nc0, 1);
    return;
end

Nc0 = find(num0 ~= 0, 1, 'last');
Nc1 = find(num1 ~= 0, 1, 'last');
cell_identity0 = cell_identity0(1:Nc0);
cell_identity1 = cell_identity1(1:Nc1);

[cmX0, cmY0] = cellCentroidsPBC(v0, inn0, num0, Nc0, Lx, Ly);
[cmX1, cmY1] = cellCentroidsPBC(v1, inn1, num1, Nc1, Lx, Ly);

neigh = buildNeighbors(inn0, num0, Nc0);

% LoadData.m returns cell_identity as a `string` array, not a cell array
% of char -- parenthesis indexing (matching compute_MSD_cellID.m's proven
% pattern), and a plain containers.Map (no explicit KeyType) so it
% infers/accepts string keys directly rather than requiring char.
id1_index = containers.Map;
for k = 1:Nc1
    id1_index(cell_identity1(k)) = k;
end

FTLE = nan(Nc0, 1);

for i = 1:Nc0
    key_i = cell_identity0(i);
    if ~isKey(id1_index, key_i)
        continue;   % this cell was extruded before it1
    end
    ii1 = id1_index(key_i);

    nb = neigh{i};
    if numel(nb) < 3
        continue;   % too few neighbors for a well-posed 2D affine fit
    end

    R0 = zeros(2, numel(nb));
    R1 = zeros(2, numel(nb));
    n_valid = 0;
    for k = 1:numel(nb)
        j = nb(k);
        key_j = cell_identity0(j);
        if ~isKey(id1_index, key_j)
            continue;   % this neighbor was extruded before it1
        end
        jj1 = id1_index(key_j);

        dx0 = cmX0(j) - cmX0(i); dx0 = dx0 - Lx * round(dx0 / Lx);
        dy0 = cmY0(j) - cmY0(i); dy0 = dy0 - Ly * round(dy0 / Ly);
        dx1 = cmX1(jj1) - cmX1(ii1); dx1 = dx1 - Lx * round(dx1 / Lx);
        dy1 = cmY1(jj1) - cmY1(ii1); dy1 = dy1 - Ly * round(dy1 / Ly);

        n_valid = n_valid + 1;
        R0(:, n_valid) = [dx0; dy0];
        R1(:, n_valid) = [dx1; dy1];
    end
    if n_valid < 3
        continue;
    end
    R0 = R0(:, 1:n_valid);
    R1 = R1(:, 1:n_valid);

    F = R1 * pinv(R0);
    C = F' * F;
    lam = max(real(eig(C)));
    if lam <= 0
        continue;
    end
    FTLE(i) = (1 / T) * log(sqrt(lam));
end

end


function [cmX, cmY] = cellCentroidsPBC(v, inn, num, Nc, Lx, Ly)
% PBC-aware per-cell centroid (log.txt: same unwrap-relative-to-first-
% vertex convention as ComputeCellColorData.m/Analysis_COM.m). Harmless
% no-op for a non-periodic mesh.
cmX = zeros(Nc, 1);
cmY = zeros(Nc, 1);
for i = 1:Nc
    vx = v(inn(i, 1:num(i)), 1);
    vy = v(inn(i, 1:num(i)), 2);
    if numel(vx) > 1
        dx = vx(2:end) - vx(1); dx = dx - Lx .* round(dx ./ Lx);
        vx(2:end) = vx(1) + dx;
        dy = vy(2:end) - vy(1); dy = dy - Ly .* round(dy ./ Ly);
        vy(2:end) = vy(1) + dy;
    end
    cmX(i) = mean(vx);
    cmY(i) = mean(vy);
end
end


function neigh = buildNeighbors(inn, num, Nc)
% Edge-sharing (>=2 shared vertices) neighbor list, built from a
% vertex->incident-cells map in O(total vertex references) rather than
% the O(Nc^2) all-pairs-with-intersect scan this was originally written
% with -- matters once Nc is in the thousands.
maxV = max(max(inn(1:Nc, :)));
vertex_cells = cell(maxV, 1);
for i = 1:Nc
    vids = inn(i, 1:num(i));
    for k = 1:numel(vids)
        vertex_cells{vids(k)}(end+1) = i;
    end
end

I = zeros(0, 1);
J = zeros(0, 1);
for v = 1:maxV
    cs = vertex_cells{v};
    nk = numel(cs);
    if nk < 2
        continue;
    end
    for a = 1:nk-1
        for b = a+1:nk
            I(end+1) = cs(a); %#ok<AGROW>
            J(end+1) = cs(b); %#ok<AGROW>
        end
    end
end

S = sparse(I, J, 1, Nc, Nc);
S = S + S';
[ii, jj] = find(S >= 2);

neigh = cell(Nc, 1);
for k = 1:numel(ii)
    neigh{ii(k)}(end+1) = jj(k);
end
end

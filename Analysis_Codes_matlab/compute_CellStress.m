function [Pressure, ShearStress] = compute_CellStress(v, inn, num, biochemdata, Lx, Ly)
% COMPUTE_CELLSTRESS  Per-cell Pressure/ShearStress, one value per live
% cell, for EVERY cell in the tissue (no radius restriction, no
% area-weighted aggregation) -- for coloring a full-tissue snapshot
% (Movie_Code.m's colorBy = 'Pressure'/'ShearStress'). This is distinct
% from Calculate_Total_Stress.m/compute_StressTensor_series.m, which
% compute a single radius-restricted, area-weighted-average scalar per
% frame for PlotAnalysis.m's "Pressure"/"ShearStress" time-series panels.
%
%   [Pressure, ShearStress] = compute_CellStress(v, inn, num, biochemdata, Lx, Ly)
%
% Uses the same per-cell virial stress tensor as Stress.f90 (Fortran) /
% Calculate_Total_Stress.m (MATLAB):
%   sigma = 2*lambda*(area-A0)*I + [(2*beta*perimeter+gamma)/(2*area)] * sum_edges(edge (x) edge)/|edge|
%   Pressure(i)    = -(sigma(1,1) + sigma(2,2)) / 2
%   ShearStress(i) =  sigma(1,2)
% matching compute_StressTensor_series.m's sign convention exactly. The
% formula (including the beta/gamma rescaling below) was verified against
% a standalone Fortran port of Stress.f90's exact computation on a test
% polygon -- matched to 13+ significant figures (log.txt).
%
% BUGFIX vs. Calculate_Total_Stress.m (log.txt): that function reads
% beta/gamma directly from para_Simulation.dat -- the file's RAW values.
% But allocation.f90::read_input rescales them exactly ONCE right after
% reading (beta = beta/(lambda*Ao); gamm = gamm/(lambda*Ao^1.5)) before
% Force.f90/Stress.f90 ever use them -- so the values the simulation's own
% stress tensor (Stress.f90) actually uses are the RESCALED ones, not the
% file's raw numbers (confirmed directly against allocation.f90's
% read_input). This function applies the same rescaling.
% Calculate_Total_Stress.m/compute_StressTensor_series.m have the same
% pre-existing bug and were NOT touched here (out of scope for this
% change -- flagged separately; for the repo's default para_Simulation.dat
% it happens to be a no-op since lambda=Ao=1, so existing Pressure/
% ShearStress panels haven't shown it).
%
% Also reproduces Force.f90's per-cell DYNAMIC beta when if_RhoROCK or
% if_active_contractility is active -- every cell then has its OWN beta,
% not a single global scalar (both derived from Myosin(ic), biochemdata's
% 3rd column). Force.f90 evaluates if_RhoROCK first and
% if_active_contractility second, unconditionally overwriting beta if both
% are on (a pre-existing, documented flag conflict -- see FEATURE_IDEAS.txt
% item B2) -- this function reproduces that same overwrite order so its
% result matches whatever the simulation actually used, rather than
% silently resolving the conflict. if_RhoROCK's if_coupling_noise term
% adds a per-timestep random draw Fortran never records, so it can't be
% reproduced post-hoc -- the noise-free formula is used in that case,
% the only inexactness here.
%
% biochemdata : numdim x 3 array from LoadData (Rho, ROCK, Myosin); only
%               needed when if_RhoROCK or if_active_contractility is on.
% Lx, Ly      : mesh box dimensions (para_MeshDims.dat), for PBC-aware
%               per-cell unwrap (same convention as ComputeCellColorData.m).
%
% Returns Nc x 1 arrays (Nc = number of currently-live cells), ready to
% hand straight to TisuePlot alongside v/inn/num.

Nc = find(num ~= 0, 1, 'last');

p1 = ReadPara1Params('../para_Simulation.dat');
A0     = p1.Ao;
lambda = p1.lambda;
beta_internal  = p1.beta / (lambda * A0);
gamma_internal = p1.gamm / (lambda * A0^1.5);

beta_local = beta_internal * ones(Nc, 1);
if p1.if_RhoROCK
    beta_local = p1.Myosin_Coupling_Strength * biochemdata(1:Nc, 3) / (lambda * A0);
end
if p1.if_active_contractility
    beta_local = biochemdata(1:Nc, 3) / (lambda * A0);
end

Pressure = zeros(Nc, 1);
ShearStress = zeros(Nc, 1);

for i = 1:Nc
    n = num(i);
    vx = v(inn(i, 1:n), 1);
    vy = v(inn(i, 1:n), 2);

    % PBC-aware unwrap (log.txt) -- same convention as
    % ComputeCellColorData.m/Calculate_Total_Stress.m.
    if n > 1
        dx = vx(2:end) - vx(1); dx = dx - Lx .* round(dx ./ Lx);
        vx(2:end) = vx(1) + dx;
        dy = vy(2:end) - vy(1); dy = dy - Ly .* round(dy ./ Ly);
        vy(2:end) = vy(1) + dy;
    end

    xp = vx([2:n 1]);
    yp = vy([2:n 1]);
    area = abs(sum(vx .* yp - xp .* vy)) / 2;

    dx_e = xp - vx;
    dy_e = yp - vy;
    len_e = sqrt(dx_e.^2 + dy_e.^2);
    perimeter = sum(len_e);

    term1 = 2.0 * lambda * (area - A0);
    term2 = (2.0 * beta_local(i) * perimeter + gamma_internal) / (2.0 * area);

    s11 = term1 + term2 * sum(dx_e.^2 ./ len_e);
    s12 = term2 * sum(dx_e .* dy_e ./ len_e);
    s22 = term1 + term2 * sum(dy_e.^2 ./ len_e);

    Pressure(i)    = -(s11 + s22) / 2;
    ShearStress(i) = s12;
end

end

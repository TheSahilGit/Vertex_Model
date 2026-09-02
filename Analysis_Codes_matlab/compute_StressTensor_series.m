function [time, ShearStress_local, Pressure] = compute_StressTensor_series(itList, nrun, radius, progressFcn)
% COMPUTE_STRESSTENSOR_SERIES  Locally-recomputed shear stress and pressure
% vs. time, within a given radius of the tissue's center of mass.
%
%   [time, ShearStress_local, Pressure] = compute_StressTensor_series(itList, nrun, radius, progressFcn)
%
% progressFcn : optional function handle progressFcn(k, n), called after
%          each of the n requested frames is processed (used by
%          PlotAnalysis.m to print live progress); omit or pass [] for none.
%
% Neither of these is stored by the Fortran code as a time series over
% arbitrary cells-within-a-radius, so both are recomputed here, once per
% requested frame, by reusing Calculate_Total_Stress.m's per-cell stress
% tensor (the same one Stress.f90 computes on the Fortran side) -- this
% is the "recalculate if needed" path. (For the CHEAP, already-computed,
% whole-tissue ShearStress vs. time -- one value per Fortran timestep,
% no geometry recompute needed -- use LoadData's `all_end_data`/ShearStress
% output instead; that's what PlotAnalysis.m's "ShearStress" panel uses by
% default. This function is for the spatially-localized, pressure-capable
% version.)
%
% ShearStress_local(k) = total_stress_tensor(1,2) (the shear/off-diagonal
%   component), matching Stress.f90's ShearStress = TotalSigma(1,2).
% Pressure(k) = -(total_stress_tensor(1,1) + total_stress_tensor(2,2))/2
%   (negative of the mean normal stress -- positive Pressure means net
%   compression, the usual continuum-mechanics sign convention).

if nargin < 3 || isempty(radius)
    radius = 10;
end
if nargin < 4
    progressFcn = [];
end

p1 = ReadPara1Params("../para1_in.dat");
time = itList * p1.dt;

ShearStress_local = zeros(size(itList));
Pressure = zeros(size(itList));
n = numel(itList);

for k = 1:n
    it = itList(k);
    [Lx, Ly, v, inn, num] = LoadData(it, nrun);

    [total_stress_tensor, ~] = Calculate_Total_Stress(Lx, Ly, v, inn, num, radius);

    ShearStress_local(k) = total_stress_tensor(1,2);
    Pressure(k) = -(total_stress_tensor(1,1) + total_stress_tensor(2,2)) / 2;

    if ~isempty(progressFcn); progressFcn(k, n); end
end

end

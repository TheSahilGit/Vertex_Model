function [time, Energy, ShearStress, T1_count, T2_count, cumsum_T1, cumsum_T2] = LoadGlobalTimeSeries(nrun)
% LOADGLOBALTIMESERIES  Read the whole-tissue scalar series written every
% Fortran timestep (Energy, ShearStress, T1_count, T2_count), plus their
% T1/T2 cumulative sums. One entry per timestep (dt apart) -- these are
% NOT tied to it_dump, unlike the full mesh snapshots LoadData reads.
%
%   [time, Energy, ShearStress, T1_count, T2_count, cumsum_T1, cumsum_T2] = LoadGlobalTimeSeries(nrun)
%
% nrun (optional, default 1): matches allocation.f90's write_output, which
% writes these under a "nrun2_" prefix for nrun==2 (a restart run must not
% overwrite the original nrun==1 run's own Energy/ShearStress/T1_count/
% T2_count files -- see log.txt).
%
% These four files are now written incrementally as a run progresses (see
% allocation.f90's summary_dump_interval), not just once at completion, so
% this is safe to call on a still-running simulation -- each array will
% simply be shorter than totT until the run catches up or finishes.

if nargin < 1 || isempty(nrun)
    nrun = 1;
end

p1 = ReadPara1Params("../para_Simulation.dat");

if nrun == 1
    prefix = '';
elseif nrun == 2
    prefix = 'nrun2_';
else
    error('LoadGlobalTimeSeries:badNrun', 'nrun must be 1 or 2, got %d', nrun);
end

Energy      = readSeries(['../data/' prefix 'Energy.dat']);
ShearStress = readSeries(['../data/' prefix 'ShearStress.dat']);
T1_count    = readSeries(['../data/' prefix 'T1_count.dat']);
T2_count    = readSeries(['../data/' prefix 'T2_count.dat']);

cumsum_T1 = cumsum(T1_count);
cumsum_T2 = cumsum(T2_count);

n = max([numel(Energy), numel(ShearStress), numel(T1_count), numel(T2_count)]);
time = (1:n) * p1.dt;

end


function data = readSeries(fname)
data = [];
if isfile(fname)
    fid = fopen(fname, 'r');
    fread(fid, 1, 'float32');   % dummy Fortran unformatted record-length header
    data = fread(fid, Inf, 'float64');
    fclose(fid);
else
    warning('LoadGlobalTimeSeries:missingFile', 'Missing file: %s', fname);
end
end

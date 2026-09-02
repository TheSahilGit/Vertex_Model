function [time, Energy, ShearStress, T1_count, T2_count, cumsum_T1, cumsum_T2] = LoadGlobalTimeSeries()
% LOADGLOBALTIMESERIES  Read the whole-tissue scalar series written every
% Fortran timestep (Energy, ShearStress, T1_count, T2_count), plus their
% T1/T2 cumulative sums. One entry per timestep (dt apart) -- these are
% NOT tied to it_dump, unlike the full mesh snapshots LoadData reads.
%
%   [time, Energy, ShearStress, T1_count, T2_count, cumsum_T1, cumsum_T2] = LoadGlobalTimeSeries()
%
% These four files are now written incrementally as a run progresses (see
% allocation.f90's summary_dump_interval), not just once at completion, so
% this is safe to call on a still-running simulation -- each array will
% simply be shorter than totT until the run catches up or finishes.

p1 = ReadPara1Params("../para1_in.dat");

Energy      = readSeries('../data/Energy.dat');
ShearStress = readSeries('../data/ShearStress.dat');
T1_count    = readSeries('../data/T1_count.dat');
T2_count    = readSeries('../data/T2_count.dat');

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

function itList = GetSnapshotItList(it_dump, itStart, itEnd, itInterval)
% GETSNAPSHOTITLIST  Build a list of Fortran timesteps that actually have a
% full mesh-snapshot dump file on disk (inn_*/num_*/v_*/force_*/Myosin_*/
% cell_identity_*), given the write cadence.
%
%   itList = GetSnapshotItList(it_dump, itStart, itEnd, itInterval)
%
% vertexmain.f90 only calls write_output when modulo(it,it_dump)==0, or
% it==1, or it==2 -- so those are the ONLY it values with a dump file.
% Requesting anything else (e.g. an itInterval that isn't a multiple of
% it_dump) would try to load a file that was never written.
%
% itInterval (optional, default it_dump) is how coarsely to sample --
% must be a multiple of it_dump; if it isn't, it is rounded UP to the
% nearest multiple (with a warning), since sampling every raw it_dump can
% be far more points than needed for a plot over a long run.

if nargin < 4 || isempty(itInterval)
    itInterval = it_dump;
end

if mod(itInterval, it_dump) ~= 0
    itInterval = ceil(itInterval / it_dump) * it_dump;
    warning('GetSnapshotItList:roundedInterval', ...
        ['itInterval rounded up to %d to stay a multiple of it_dump=%d ' ...
         '-- only multiples of it_dump (plus it=1,2) are ever written to disk.'], ...
        itInterval, it_dump);
end

firstDump = max(it_dump, ceil(itStart / it_dump) * it_dump);
if firstDump <= itEnd
    itList = firstDump:itInterval:itEnd;
else
    itList = [];
end

% it==1 and it==2 always get a dump too (see vertexmain.f90); include them
% if they fall in the requested range and aren't already covered.
extra = intersect([1 2], itStart:itEnd);
itList = unique([extra, itList]);

end

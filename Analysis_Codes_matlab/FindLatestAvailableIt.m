function itLatest = FindLatestAvailableIt(nrun, it_dump, itUpperBound)
% FINDLATESTAVAILABLEIT  Find the largest it <= itUpperBound that actually
% has a full mesh-snapshot dump file on disk, so a still-running
% simulation's data can be plotted without knowing exactly how far it has
% progressed (or without guessing wrong and hitting a missing-file error).
%
%   itLatest = FindLatestAvailableIt(nrun, it_dump, itUpperBound)
%
% Returns 0 if no snapshot at all is found (e.g. the run hasn't reached
% it==1 or it==2 yet -- extremely unlikely in practice, but handled).

if nrun == 1
    fmt = '../data/inn_%08d.dat';
elseif nrun == 2
    fmt = '../data/nrun2_inn_%08d.dat';
else
    error('FindLatestAvailableIt:badNrun', 'nrun must be 1 or 2, got %d', nrun);
end

itLatest = 0;

firstCandidate = floor(itUpperBound / it_dump) * it_dump;
for it = firstCandidate:-it_dump:it_dump
    if isfile(sprintf(fmt, it))
        itLatest = it;
        return;
    end
end

for it = [2 1]
    if it <= itUpperBound && isfile(sprintf(fmt, it))
        itLatest = it;
        return;
    end
end

end

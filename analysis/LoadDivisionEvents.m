function [Division_it, Division_x, Division_y, Division_id1, Division_id2] = LoadDivisionEvents(nrun)
% LOADDIVISIONEVENTS  Read the cell-division spatial/cell-identity event
% log written incrementally by Proliferation_Core (Proliferation.f90;
% log.txt) -- same binary event-log convention as LoadT1T2Events.m (one
% Fortran unformatted record per actual event, not one entry per
% timestep).
%
%   [Division_it, Division_x, Division_y, Division_id1, Division_id2] = LoadDivisionEvents(nrun)
%
% Division_events.dat, one record per division, 5 real*8 fields:
%   it  x  y  id1  id2
% (x,y) is the PBC-wrapped midpoint of the two new vertices created by the
% split (mirrors T1's own "flipping edge midpoint" convention); id1 is the
% identity of the daughter that kept the mother's own index/identity, id2
% the identity of the brand-new daughter cell -- so both post-division
% cells can be found later purely by identity, the same as T1/T2 events.
%
% A still-running (or if_cell_division-free) simulation simply hasn't
% written the file yet, or has written an empty one; both return empty
% arrays, not an error.
%
% Division_id1/Division_id2 are Nx1 cell arrays of 'cell_<N>' strings
% (reconstructed from the numeric IDs here, so the rest of this toolkit --
% main_MovieCode.m's strcmp against cell_identity -- never needs to know about
% the binary encoding).

if nrun == 1
    fname = '../data/Division_events.dat';
else
    fname = '../data/nrun2_Division_events.dat';
end

[Division_it, Division_x, Division_y, ids] = readEventRecords(fname, 2);
Division_id1 = numToIdentity(ids(:, 1));
Division_id2 = numToIdentity(ids(:, 2));

end


function [it, x, y, idnum] = readEventRecords(fname, n_ids)
% Read every Fortran unformatted record from `fname` as one flat
% "3 + n_ids" real*8 fields (it, x, y, then n_ids identity numbers), fully
% vectorized -- identical technique to LoadT1T2Events.m's own helper of
% the same name (kept as a private copy here rather than shared, matching
% that file's own self-contained-loader convention).
it = []; x = []; y = []; idnum = zeros(0, n_ids);
if ~isfile(fname)
    return;
end

n_fields = 3 + n_ids;
payload_bytes = 8 * n_fields;
rec_bytes = 4 + payload_bytes + 4;

fid = fopen(fname, 'r');
raw = fread(fid, Inf, 'uint8=>uint8');
fclose(fid);

if isempty(raw)
    return;
end
if mod(numel(raw), rec_bytes) ~= 0
    warning('LoadDivisionEvents:badFile', ...
        '%s size (%d bytes) is not a multiple of the expected record size (%d) -- ignoring.', ...
        fname, numel(raw), rec_bytes);
    return;
end

n = numel(raw) / rec_bytes;
raw = reshape(raw, rec_bytes, n);
payload = raw(5:4+payload_bytes, :);
vals = typecast(payload(:), 'double');
vals = reshape(vals, n_fields, n)';

it = vals(:, 1);
x = vals(:, 2);
y = vals(:, 3);
idnum = vals(:, 4:end);

end


function c = numToIdentity(idnum)
% Reconstruct 'cell_<N>' strings from the numeric identity encoding --
% identical convention to LoadT1T2Events.m's own helper of the same name.
c = cell(size(idnum));
for k = 1:numel(idnum)
    n = round(idnum(k));
    if n <= 0
        c{k} = 'none';
    else
        c{k} = sprintf('cell_%d', n);
    end
end
end

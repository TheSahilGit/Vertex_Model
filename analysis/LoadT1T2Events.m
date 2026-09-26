function [T1_it, T1_x, T1_y, T1_ids, T2_it, T2_x, T2_y, T2_extruded_id, T2_nbr_ids] = LoadT1T2Events(nrun)
% LOADT1T2EVENTS  Read the T1/T2 spatial/cell-identity event logs written
% incrementally by Do_T1/Do_T2 (T1_swap.f90/T2_swap.f90; log.txt).
%
%   [T1_it, T1_x, T1_y, T1_ids, T2_it, T2_x, T2_y, T2_extruded_id, T2_nbr_ids] = LoadT1T2Events(nrun)
%
% Unlike every other time series this toolkit reads (Energy, T1_count,
% ...), these files have ONE FORTRAN UNFORMATTED RECORD PER ACTUAL EVENT
% (not one entry per timestep, and not one single big record for the
% whole run) -- sparse and irregular, since T1/T2 events don't happen
% every step. Each record is all-real*8 (matching this codebase's habit
% of using real*8 even for integer-valued quantities, so no record ever
% mixes types), with each cell identity stored as the NUMERIC SUFFIX of
% its persistent 'cell_<N>' string (see CellIdNum, allocation.f90) rather
% than the string itself, to keep the log compact -- 0.0 means 'cell_0'
% or an empty/padding slot (this codebase's usual "0 = empty" idiom).
% A still-running (or T1/T2-free) simulation simply hasn't written the
% file yet, or has written an empty one; both return empty arrays, not
% an error.
%
% T1_events.dat, one record per T1 flip, 7 real*8 fields:
%   it  x  y  id1  id2  id3  id4
% (x,y) is the flipping edge's midpoint (PBC-wrapped); id1..id4 are the
% persistent identity of the (up to 4) affected cells -- the 2 "losing"
% cells that shared the collapsing edge and the 2 "gaining" cells that
% only touched one endpoint -- 0.0 padded if fewer.
%
% T2_events.dat, one record per T2 extrusion, 10 real*8 fields:
%   it  x  y  extruded_id  nbr1  nbr2  nbr3  nbr4  nbr5  nbr6
% (x,y) is the extruded triangular cell's centroid; extruded_id is its own
% persistent identity (captured before removal); nbr1..nbr6 are its
% neighbors' identities (typically 3, the triangle's edge-sharing
% neighbors; up to 3 more vertex-only neighbors), 0.0 padded.
%
% T1_ids is an Nx4 cell array of 'cell_<N>'/'none' strings (reconstructed
% from the numeric IDs here, so the rest of this toolkit -- Movie_Code.m's
% strcmp against cell_identity -- never needs to know about the binary
% encoding); T2_nbr_ids is an Nx6 cell array; T2_extruded_id is an Nx1
% cell array.

if nrun == 1
    fnameT1 = '../data/T1_events.dat';
    fnameT2 = '../data/T2_events.dat';
else
    fnameT1 = '../data/nrun2_T1_events.dat';
    fnameT2 = '../data/nrun2_T2_events.dat';
end

[T1_it, T1_x, T1_y, T1_idnum] = readEventRecords(fnameT1, 4);
T1_ids = numToIdentity(T1_idnum);

[T2_it, T2_x, T2_y, T2_idnum] = readEventRecords(fnameT2, 7);
T2_extruded_id = numToIdentity(T2_idnum(:, 1));
T2_nbr_ids = numToIdentity(T2_idnum(:, 2:7));

end


function [it, x, y, idnum] = readEventRecords(fname, n_ids)
% Read every Fortran unformatted record from `fname` as one flat
% "3 + n_ids" real*8 fields (it, x, y, then n_ids identity numbers),
% fully vectorized: every record has exactly the same byte size, so the
% whole file is read once as raw bytes, reshaped into one column per
% record, and the leading/trailing 4-byte gfortran record-length markers
% are sliced off before typecasting the remaining bytes to double.
it = []; x = []; y = []; idnum = zeros(0, n_ids);
if ~isfile(fname)
    return;
end

n_fields = 3 + n_ids;
payload_bytes = 8 * n_fields;
rec_bytes = 4 + payload_bytes + 4;   % leading marker + payload + trailing marker

fid = fopen(fname, 'r');
raw = fread(fid, Inf, 'uint8=>uint8');
fclose(fid);

if isempty(raw)
    return;
end
if mod(numel(raw), rec_bytes) ~= 0
    warning('LoadT1T2Events:badFile', ...
        '%s size (%d bytes) is not a multiple of the expected record size (%d) -- ignoring.', ...
        fname, numel(raw), rec_bytes);
    return;
end

n = numel(raw) / rec_bytes;
raw = reshape(raw, rec_bytes, n);
payload = raw(5:4+payload_bytes, :);          % drop the two 4-byte markers
vals = typecast(payload(:), 'double');        % n_fields*n doubles, field-major within each record
vals = reshape(vals, n_fields, n)';           % n x n_fields

it = vals(:, 1);
x = vals(:, 2);
y = vals(:, 3);
idnum = vals(:, 4:end);

end


function c = numToIdentity(idnum)
% Reconstruct 'cell_<N>' strings from the numeric identity encoding used
% in the event log; 0 (or anything that isn't a real cell number) becomes
% 'none', matching the padding convention documented above.
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

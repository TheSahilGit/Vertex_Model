clear; clc; 
%close all;

% ==================== options ====================
nrun = 1;
itList = (1000);              % list of Fortran timesteps to render as frames
outFile = "Movie_test.avi";
frameRate = 1;

% Which per-cell field to color the tissue by. One of:
%   'Force' (default), 'Motility', 'Myosin', 'Rho', 'ROCK', 'Area',
%   'Perimeter', 'ShapeFactor', 'NumVertices', 'Pressure', 'ShearStress',
%   'FTLE', 'Eta'
% -- see ComputeCellColorData.m for what each one computes, except 'FTLE'
% which is handled separately below (see ftle_lookahead) since it's
% inherently a two-snapshot quantity, not a single-frame field.
colorBy = 'Force';

norm_flag = 'data';   % 'data' | '01' | 'custom'
norm_range = [];      % only used when norm_flag == 'custom', e.g. [0 2]

% For FTLE
% norm_flag = 'data';   % 'data' | '01' | 'custom'
% norm_range = [];      % only used when norm_flag == 'custom', e.g. [0 2]

% ---- T1/T2 event overlay (log.txt: spatial/cell-identity tracking, read
% from T1_events.dat/T2_events.dat -- see LoadT1T2Events.m) ----
show_T1T2_events = false;   % overlay recent event markers + highlight the cells involved
t1t2_fade_window = 10000;   % "recent" = within this many `it` units of the frame being
                           % drawn; markers fade linearly to invisible over this window,
                           % cell highlights are on/off (not faded) within it, for speed.

% ---- Division event overlay (log.txt: spatial/cell-identity tracking,
% read from Division_events.dat -- see LoadDivisionEvents.m) ----
show_division_events = false;   % overlay recent division markers + highlight both daughter cells
division_fade_window = 10000;   % same "recent"/fade convention as t1t2_fade_window above.

% ---- FTLE (colorBy = 'FTLE' only; log.txt, see compute_FTLE.m) ----
% Each frame at time `it` shows the spatial FTLE field computed forward
% to `it + ftle_lookahead` -- i.e. this many `it` units of "look-ahead"
% per frame, not a fixed reference time. Needs a real, non-tiny gap
% (comparable to how long it takes cells to actually rearrange) to be
% meaningful -- too small and every cell's local neighborhood is still
% near-identical, giving a near-uniform, uninformative field.
ftle_lookahead = 10000;

% ---- Renderer (log.txt) ----
% 'opengl' (default): hardware-accelerated, memory-stable, ~3x faster --
% measured directly on a 10,000-cell/150-frame stress test: 'painters'
% leaked ~490MB of RSS over the run and crashed intermittently (getframe:
% "A valid figure or axes handle must be specified"); opengl stayed flat
% with no crashes, and both render visually identically (transparency/
% edges). Switch to 'painters' only if opengl misbehaves on a given
% machine (e.g. no GPU/driver, headless/software-only display).
rendererMode = 'opengl';   % 'opengl' (default) | 'painters'

% ---- Partial-lattice crop (log.txt: replaces the old MovieCode_halflatt.m
% variant, which despite its name never actually did this -- it was a
% near-duplicate of this script with colorBy/nrun/itList hardcoded and
% some dead unused trailing code; this crop is a genuinely new capability,
% not a restoration of one). Cells whose (PBC-unwrapped) centroid has
% y > plottill are excluded ENTIRELY -- not drawn, and not counted in the
% colorbar's 'data' min/max -- so the kept region is treated as if it
% were the whole tissue, exactly like re-running this on a smaller mesh,
% not like zooming the camera on an unchanged color scale. e.g.
% plottill = Ly/2 keeps only the bottom half. Leave empty ([]) to keep
% the whole tissue (default, same as before this option existed).
plottill = [];
% ===================================================

para2 = load("../mesh/para_MeshDims.dat");
Lx = para2(1);
Ly = para2(2);

% Motility (etas) is only loaded if colorBy actually needs it -- an extra
% file read every other colorBy option doesn't need -- and is now loaded
% PER FRAME, inside the loop below (LoadMotility), not once here.
%
% BUGFIX (log.txt): allocation.f90's motility output used to be written
% ONCE, at it==1, to a single motility_store.dat -- a frozen snapshot of
% the INITIAL per-vertex field. Loading it once here, outside the frame
% loop, matched that (there was only ever one snapshot to load) but meant
% every frame of a movie showed the SAME t=1 motility regardless of which
% it it was actually rendering -- wrong for if_motility_decay/
% if_motility_Eulerian (which evolve mot over time) and for
% if_cell_division (new cells born after t=1 showed motility exactly 0,
% confirmed directly: 598/598 division-created vertices had frozen value
% 0.0 vs. nonzero live mot, in a real hotspot+division run). Fixed on the
% Fortran side first (allocation.f90 now writes motility_<it>.dat every
% dump, same convention as v/inn/num/force/Myosin/cell_identity) --
% LoadMotility below reads THAT per-frame file, falling back to the old
% static motility_store.dat only for data/ directories from before this
% fix (same staleness limitation those always had).
etas = [];

% eta (friction/Langevin coefficient, log.txt) -- same lazy/per-frame
% loading convention as etas above (LoadEta reads data/eta_<it>.dat,
% written every dump since the array was introduced -- no legacy
% frozen-snapshot file exists for this one, unlike motility's
% motility_store.dat, since eta was never written before that change).
eta_field = [];

% Loaded once, reused every frame -- same pattern as `etas` above. Empty
% arrays if the flag is off, or if the simulation never had if_Do_T1/
% if_Do_T2 on (LoadT1T2Events.m returns empty for a missing file).
T1_it = []; T1_x = []; T1_y = []; T1_ids = {};
T2_it = []; T2_x = []; T2_y = []; T2_extruded_id = {}; T2_nbr_ids = {};
if show_T1T2_events
    [T1_it, T1_x, T1_y, T1_ids, T2_it, T2_x, T2_y, T2_extruded_id, T2_nbr_ids] = ...
        LoadT1T2Events(nrun);
end

Division_it = []; Division_x = []; Division_y = []; Division_id1 = {}; Division_id2 = {};
if show_division_events
    [Division_it, Division_x, Division_y, Division_id1, Division_id2] = LoadDivisionEvents(nrun);
end

% Which colormap to render colorBy with -- one central lookup (log.txt),
% picked once (not per frame: colorBy doesn't change frame-to-frame, and
% slanCM.m's interpolation isn't free). Edit GetFieldColormap.m to change
% any field's colormap later.
[cmap, isDiverging] = GetFieldColormap(colorBy);

%cmap  = slanCM("jet");

fig = figure("Position", [800 800 1000 1000], 'Color','w');
% BUGFIX (log.txt): VideoWriter requires every frame to be EXACTLY the
% same pixel size -- 'Position' above only sets the figure's *nominal*
% size, and getframe(gcf) can still come back a few pixels off between
% calls (colorbar tick-label width changing with the data range each
% frame, OS/window-manager nudging the window, etc.), which throws
% "Frame must be H by W" on whichever frame first differs from the
% first one VideoWriter locked in on. 'Resize','off' stops MATLAB/the
% window manager from resizing the figure mid-loop; imresize below is
% the actual guarantee -- every frame handed to writeVideo is forced to
% the exact size of the first one, regardless of the cause.
set(fig, 'Resize', 'off');

mov = VideoWriter(outFile);
mov.FrameRate = frameRate;
open(mov);

frameSize = [];   % [height width], locked in from the first frame

for it = itList

    clf

    [Lx, Ly, v, inn, num, forces, biochemdata, cell_identity] = LoadData(it, nrun);

    if strcmp(colorBy, 'Motility')
        etas = LoadMotility(it, nrun);
    end

    if strcmp(colorBy, 'Eta')
        eta_field = LoadEta(it, nrun);
    end

    if strcmp(colorBy, 'FTLE')
        % Two-snapshot quantity -- bypasses ComputeCellColorData.m's
        % single-frame dispatch entirely; reuses this frame's
        % already-loaded v/inn/num/cell_identity as the reference
        % snapshot instead of having compute_FTLE.m reload it.
        colordata = compute_FTLE(it, it + ftle_lookahead, nrun, v, inn, num, cell_identity);
        colorbar_string = sprintf('FTLE (look-ahead %d)', ftle_lookahead);
    else
        [colordata, colorbar_string] = ComputeCellColorData( ...
            colorBy, v, inn, num, forces, biochemdata, etas, eta_field, Lx, Ly);
    end

    % Apply plottill (see option above): zero out num for any cell whose
    % (PBC-unwrapped) centroid falls outside the kept y-range. TisuePlot.m
    % skips num(i)==0 cells both when building faces (an all-NaN face
    % row -- not drawn) and, since this change, when computing the 'data'
    % colorbar range -- so the kept region is treated as if it were the
    % whole tissue, not just a zoomed-in view of an unchanged scale.
    num_eff = num;
    Nc_frame = find(num ~= 0, 1, 'last');
    if ~isempty(plottill)
        for ic = 1:Nc_frame
            vy = v(inn(ic, 1:num(ic)), 2);
            if numel(vy) > 1
                dy = vy(2:end) - vy(1); dy = dy - Ly .* round(dy ./ Ly);
                vy(2:end) = vy(1) + dy;
            end
            if mean(vy) > plottill
                num_eff(ic) = 0;
            end
        end
    end

    % Diverging fields (Pressure/ShearStress -- decided by GetFieldColormap.m
    % based on the quantity itself, not a flag threaded through TisuePlot.m)
    % get their colorbar centered symmetrically about zero instead of the
    % field's raw (usually asymmetric) min/max, folded into an ordinary
    % 'custom' norm_flag/norm_range here -- TisuePlot.m never needs to know
    % "diverging" is a concept (log.txt). Uses num_eff so an active
    % plottill excludes hidden cells from this range too, same as
    % TisuePlot.m's own 'data' case.
    if isDiverging && strcmp(norm_flag, 'data')
        live = num_eff(1:Nc_frame) ~= 0;
        L = max(abs(colordata(live)));
        [frame_norm_flag, frame_norm_range] = deal('custom', [-L L]);
    else
        [frame_norm_flag, frame_norm_range] = deal(norm_flag, norm_range);
    end

    TisuePlot(Lx, Ly, v, inn, num_eff, colordata, colorbar_string, frame_norm_flag, frame_norm_range, cmap, rendererMode);

    % Must come AFTER TisuePlot (which unconditionally sets its own
    % axis([-4 Lx+4 -4 Ly+4]) and pbaspect([Lx/Ly 1 1]) as its last lines)
    % -- otherwise this would just get overwritten. pbaspect is
    % recomputed against plottill (not Ly) so the kept cells still render
    % at their true aspect ratio instead of looking vertically stretched.
    if ~isempty(plottill)
        ylim([-4 plottill+4])
        pbaspect([Lx/plottill 1 1])
    end

    if show_T1T2_events
        hold on;
        Overlay_T1T2_Events(it, t1t2_fade_window, Lx, Ly, v, inn, num, cell_identity, ...
            T1_it, T1_x, T1_y, T1_ids, T2_it, T2_x, T2_y, T2_nbr_ids);
    end

    if show_division_events
        hold on;
        Overlay_Division_Events(it, division_fade_window, Lx, Ly, v, inn, num, cell_identity, ...
            Division_it, Division_x, Division_y, Division_id1, Division_id2);
    end

    title(num2str(it))
    drawnow;
    F = getframe(fig);

    img = F.cdata;
    if isempty(frameSize)
        frameSize = [size(img,1), size(img,2)];
    elseif ~isequal([size(img,1), size(img,2)], frameSize)
        img = imresize(img, frameSize);
    end
    writeVideo(mov, img);

    hold off;


end

close(mov)


function Overlay_T1T2_Events(it, fade_window, Lx, Ly, v, inn, num, cell_identity, ...
    T1_it, T1_x, T1_y, T1_ids, T2_it, T2_x, T2_y, T2_nbr_ids)
% OVERLAY_T1T2_EVENTS  Draw recent T1 ('x', magenta) / T2 (filled square,
% orange) event-location markers on top of the current TisuePlot, fading
% linearly to invisible over `fade_window` (it) units, and outline every
% currently-live cell that was involved in one of those recent events
% (green, not faded -- kept binary for simplicity/speed: a cell is either
% "recently involved" or it isn't).
%
% Looks up each event's persistent cell_identity string against the
% CURRENT frame's cell_identity array to find that cell's present-day
% index (which drifts over time as T2 removes/renumbers cells) -- a cell
% that has since been extruded itself simply finds no match and is
% skipped, exactly like compute_MSD_cellID.m's cohort-tracking convention.

Nc = find(num ~= 0, 1, 'last');
highlighted = false(Nc, 1);

% ---- T1 markers ----
% "Fade" is a blend toward white as the event ages, not alpha transparency
% -- MarkerFaceColor (needed below for T2's filled square) doesn't support
% the 4-element RGBA extension that plain line Color sometimes does, so a
% single consistent technique (plain 3-element RGB, both markers) is used
% for both, rather than relying on that inconsistent MATLAB behavior.
T1_color = [1 0 1];    % magenta
T2_color = [1 0.5 0];  % orange
for k = 1:numel(T1_it)
    age = it - T1_it(k);
    if age < 0 || age > fade_window
        continue;
    end
    fade = min(0.85, age / fade_window);
    c = T1_color * (1 - fade) + [1 1 1] * fade;
    plot(T1_x(k), T1_y(k), 'x', 'Color', c, 'LineWidth', 2, 'MarkerSize', 10);
    for jj = 1:size(T1_ids, 2)
        idx = find(strcmp(cell_identity(1:Nc), T1_ids{k, jj}), 1);
        if ~isempty(idx)
            highlighted(idx) = true;
        end
    end
end

% ---- T2 markers ----
for k = 1:numel(T2_it)
    age = it - T2_it(k);
    if age < 0 || age > fade_window
        continue;
    end
    fade = min(0.85, age / fade_window);
    c = T2_color * (1 - fade) + [1 1 1] * fade;
    plot(T2_x(k), T2_y(k), 's', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 9);
    for jj = 1:size(T2_nbr_ids, 2)
        idx = find(strcmp(cell_identity(1:Nc), T2_nbr_ids{k, jj}), 1);
        if ~isempty(idx)
            highlighted(idx) = true;
        end
    end
end

% ---- highlight outline on every currently-live, recently-affected cell ----
% One combined patch (not one per cell) -- same NaN-padded Faces/Vertices
% technique as TisuePlot.m, for the same reason: cheap even for many cells.
hi = find(highlighted);
if ~isempty(hi)
    maxN = max(num(hi));
    F = NaN(numel(hi), maxN);
    Vexp = zeros(sum(num(hi)), 2);
    row = 0;
    for ii = 1:numel(hi)
        i = hi(ii);
        n = num(i);
        vids = inn(i, 1:n);
        x0 = v(vids(1), 1);
        y0 = v(vids(1), 2);
        for k = 1:n
            row = row + 1;
            dx = v(vids(k), 1) - x0; dx = dx - Lx * round(dx / Lx);
            dy = v(vids(k), 2) - y0; dy = dy - Ly * round(dy / Ly);
            Vexp(row, 1) = x0 + dx;
            Vexp(row, 2) = y0 + dy;
            F(ii, k) = row;
        end
    end
    patch('Faces', F, 'Vertices', Vexp, 'FaceColor', 'none', ...
        'EdgeColor', [0.15 0.8 0.15], 'LineWidth', 3);
end

end


function etas = LoadMotility(it, nrun)
% LOADMOTILITY  Per-vertex motility field for frame `it` (log.txt). Reads
% the per-timestep dump allocation.f90 now writes every it_dump
% (motility_<it8digit>.dat, same convention as v/inn/num/force/Myosin/
% cell_identity); falls back to the old single frozen-at-it=1
% motility_store.dat only for data/ directories from before that fix
% existed (isfile-gated, so a fixed-up Fortran run's per-frame files are
% always preferred when present).
if nrun == 1
    motFile = sprintf('../data/motility_%08d.dat', it);
    motFileLegacy = '../data/motility_store.dat';
else
    motFile = sprintf('../data/nrun2_motility_%08d.dat', it);
    motFileLegacy = '../data/nrun2_motility_store.dat';
end

if isfile(motFile)
    fid = fopen(motFile);
elseif isfile(motFileLegacy)
    fid = fopen(motFileLegacy);
else
    etas = [];
    return;
end
fread(fid, 1, 'float32');
etas = fread(fid, Inf, 'float64');
fclose(fid);
end


function eta_field = LoadEta(it, nrun)
% LOADETA  Per-vertex Langevin/friction field for frame `it` (log.txt).
% Reads the per-timestep dump allocation.f90 writes every it_dump, same
% convention as motility -- but with no legacy frozen-snapshot fallback
% (unlike LoadMotility's motility_store.dat): eta was never written to
% data/ before this array existed, so there's no old-format file to fall
% back to -- a data/ directory from before this change just has no
% eta_<it>.dat at all, and this returns [] for it.
if nrun == 1
    etaFile = sprintf('../data/eta_%08d.dat', it);
else
    etaFile = sprintf('../data/nrun2_eta_%08d.dat', it);
end

if ~isfile(etaFile)
    eta_field = [];
    return;
end
fid = fopen(etaFile);
fread(fid, 1, 'float32');
eta_field = fread(fid, Inf, 'float64');
fclose(fid);
end


function Overlay_Division_Events(it, fade_window, Lx, Ly, v, inn, num, cell_identity, ...
    Division_it, Division_x, Division_y, Division_id1, Division_id2)
% OVERLAY_DIVISION_EVENTS  Draw recent division-location markers (circle,
% cyan) on top of the current TisuePlot, fading linearly to invisible over
% `fade_window` (it) units -- same fade/lookup technique as
% Overlay_T1T2_Events, distinct marker shape and color (cyan, not magenta/
% orange) so both overlays can be shown together without confusion -- and
% outline every currently-live daughter cell from a recent division
% (cyan, not faded).
%
% Looks up each event's persistent cell_identity strings (Division_id1,
% the daughter that kept the mother's own index/identity; Division_id2,
% the brand-new daughter) against the CURRENT frame's cell_identity array,
% same cohort-tracking convention as Overlay_T1T2_Events/
% compute_MSD_cellID.m -- a daughter that has since divided again or been
% extruded simply finds no match and is skipped.

Nc = find(num ~= 0, 1, 'last');
highlighted = false(Nc, 1);

Division_color = [0 0.75 0.75];  % cyan -- distinct from T1 (magenta), T2
                                  % (orange), and the T1/T2 highlight outline (green)
for k = 1:numel(Division_it)
    age = it - Division_it(k);
    if age < 0 || age > fade_window
        continue;
    end
    fade = min(0.85, age / fade_window);
    c = Division_color * (1 - fade) + [1 1 1] * fade;
    plot(Division_x(k), Division_y(k), 'o', 'Color', c, 'LineWidth', 2, 'MarkerSize', 9);

    idx1 = find(strcmp(cell_identity(1:Nc), Division_id1{k}), 1);
    if ~isempty(idx1)
        highlighted(idx1) = true;
    end
    idx2 = find(strcmp(cell_identity(1:Nc), Division_id2{k}), 1);
    if ~isempty(idx2)
        highlighted(idx2) = true;
    end
end

% ---- highlight outline on every currently-live, recently-divided cell ----
% Same combined-patch technique as Overlay_T1T2_Events/TisuePlot.m.
hi = find(highlighted);
if ~isempty(hi)
    maxN = max(num(hi));
    F = NaN(numel(hi), maxN);
    Vexp = zeros(sum(num(hi)), 2);
    row = 0;
    for ii = 1:numel(hi)
        i = hi(ii);
        n = num(i);
        vids = inn(i, 1:n);
        x0 = v(vids(1), 1);
        y0 = v(vids(1), 2);
        for k = 1:n
            row = row + 1;
            dx = v(vids(k), 1) - x0; dx = dx - Lx * round(dx / Lx);
            dy = v(vids(k), 2) - y0; dy = dy - Ly * round(dy / Ly);
            Vexp(row, 1) = x0 + dx;
            Vexp(row, 2) = y0 + dy;
            F(ii, k) = row;
        end
    end
    patch('Faces', F, 'Vertices', Vexp, 'FaceColor', 'none', ...
        'EdgeColor', [0 0.6 0.6], 'LineWidth', 3);
end

end

clear; clc; close all;

% ==================== options ====================
nrun = 1;
itList = (300000);              % list of Fortran timesteps to render as frames
outFile = "Movie_test.avi";
frameRate = 1;

% Which per-cell field to color the tissue by. One of:
%   'Force' (default), 'Motility', 'Myosin', 'Rho', 'ROCK', 'Area',
%   'Perimeter', 'ShapeFactor', 'NumVertices'
% -- see ComputeCellColorData.m for what each one computes.
colorBy = 'ShapeFactor';

norm_flag = 'data';   % 'data' | '01' | 'custom'
norm_range = [];      % only used when norm_flag == 'custom', e.g. [0 2]

% ---- T1/T2 event overlay (log.txt: spatial/cell-identity tracking, read
% from T1_events.dat/T2_events.dat -- see LoadT1T2Events.m) ----
show_T1T2_events = false;   % overlay recent event markers + highlight the cells involved
t1t2_fade_window = 10000;   % "recent" = within this many `it` units of the frame being
                           % drawn; markers fade linearly to invisible over this window,
                           % cell highlights are on/off (not faded) within it, for speed.
% ===================================================

para2 = load("../para_MeshDims.dat");
Lx = para2(1);
Ly = para2(2);

% Only load motility_store.dat if it's actually going to be used --
% it's a separate file read that every other colorBy option doesn't need.
%
% BUGFIX (log.txt): this was hardcoded to 'motility_store.dat' regardless
% of nrun. allocation.f90 writes it under a "nrun2_" prefix for nrun==2
% (a restart run must not overwrite the original nrun==1 run's own
% motility_store.dat) -- so for nrun=2 this was silently reading the
% ORIGINAL run's motility field, not the current one. Concretely: this
% run's own if_motility_hotspot pattern (nrun2_motility_store.dat, 2177
% distinct values) was masked by a stale nrun=1 file that happened to be
% completely uniform (motility_store.dat, exactly etas_max everywhere --
% from an earlier run with if_motility_hotspot off), which is exactly why
% the 'Motility' coloring showed no spatial variation at all while
% 'Force' coloring (already correctly nrun-aware via LoadData.m) showed
% the hotspots in the right place.
etas = [];
if strcmp(colorBy, 'Motility')
    if nrun == 1
        motFile = '../data/motility_store.dat';
    else
        motFile = '../data/nrun2_motility_store.dat';
    end
    fid = fopen(motFile);
    fread(fid, 1, 'float32');
    etas = fread(fid, Inf, 'float64');
    fclose(fid);
end

% Loaded once, reused every frame -- same pattern as `etas` above. Empty
% arrays if the flag is off, or if the simulation never had if_Do_T1/
% if_Do_T2 on (LoadT1T2Events.m returns empty for a missing file).
T1_it = []; T1_x = []; T1_y = []; T1_ids = {};
T2_it = []; T2_x = []; T2_y = []; T2_extruded_id = {}; T2_nbr_ids = {};
if show_T1T2_events
    [T1_it, T1_x, T1_y, T1_ids, T2_it, T2_x, T2_y, T2_extruded_id, T2_nbr_ids] = ...
        LoadT1T2Events(nrun);
end

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

    [colordata, colorbar_string] = ComputeCellColorData( ...
        colorBy, v, inn, num, forces, biochemdata, etas, Lx, Ly);

    TisuePlot(Lx, Ly, v, inn, num, colordata, colorbar_string, norm_flag, norm_range);

    if show_T1T2_events
        hold on;
        Overlay_T1T2_Events(it, t1t2_fade_window, Lx, Ly, v, inn, num, cell_identity, ...
            T1_it, T1_x, T1_y, T1_ids, T2_it, T2_x, T2_y, T2_nbr_ids);
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

clear; clc; close all;

% ==================== options ====================
nrun = 1;
itList = (500000);              % list of Fortran timesteps to render as frames
outFile = "Movie_test.avi";
frameRate = 1;

% Which per-cell field to color the tissue by. One of:
%   'Force' (default), 'Motility', 'Myosin', 'Rho', 'ROCK', 'Area',
%   'Perimeter', 'ShapeFactor', 'NumVertices'
% -- see ComputeCellColorData.m for what each one computes.
colorBy = 'Area';

norm_flag = 'data';   % 'data' | '01' | 'custom'
norm_range = [];      % only used when norm_flag == 'custom', e.g. [0 2]
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
    
    [Lx, Ly, v, inn, num, forces, biochemdata] = LoadData(it, nrun);

    [colordata, colorbar_string] = ComputeCellColorData( ...
        colorBy, v, inn, num, forces, biochemdata, etas);

    TisuePlot(Lx, Ly, v, inn, num, colordata, colorbar_string, norm_flag, norm_range);

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

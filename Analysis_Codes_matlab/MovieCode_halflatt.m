clear; clc; close all;

% Renders a Motility-colored movie. Was a standalone script with its own
% local TisuePlot() (raw, non-PBC-aware per-cell polyshape/plot loop) --
% migrated to the shared, PBC-aware ComputeCellColorData.m/TisuePlot.m
% toolkit (see Movie_Code.m, which already made this move) so this stops
% rendering the "wired" straddling-cell artifact on a periodic mesh
% (log.txt).

% ==================== options ====================
nrun = 2;
itList = (620000);              % list of Fortran timesteps to render as frames
outFile = "Movie1_test.avi";
frameRate = 1;
colorBy = 'Motility';
norm_flag = 'data';
norm_range = [];
% ===================================================

para2 = load("../para_MeshDims.dat");
Lx = para2(1);
Ly = para2(2);

% BUGFIX (log.txt, see Movie_Code.m): motility_store.dat is written under
% a "nrun2_" prefix for nrun==2 (a restart run must not overwrite the
% original nrun==1 run's own file).
if nrun == 1
    motFile = '../data/motility_store.dat';
else
    motFile = '../data/nrun2_motility_store.dat';
end
fid = fopen(motFile);
fread(fid, 1, 'float32');
etas = fread(fid, Inf, 'float64');
fclose(fid);

fig = figure('Position', [100 100 800 800], 'Color', 'w');
set(fig, 'Resize', 'off');   % BUGFIX (log.txt): see Movie_Code.m

mov = VideoWriter(outFile);
mov.FrameRate = frameRate;
open(mov);

frameSize = [];
for it = itList

    clf

    [Lx, Ly, v, inn, num, forces, biochemdata] = LoadData(it, nrun);

    [colordata, colorbar_string] = ComputeCellColorData( ...
        colorBy, v, inn, num, forces, biochemdata, etas, Lx, Ly);

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


fname_energy = sprintf('../data/Energy.dat');
fid = fopen(fname_energy);
dum4 = fread(fid,1,'float32');
energy = fread(fid,100000000,'float64');


fname_ST = sprintf('../data/ShearStress.dat');
fid = fopen(fname_ST);
dum4 = fread(fid,1,'float32');
ShearStress = fread(fid,100000000,'float64');


fname_T1 = sprintf('../data/T1_count.dat');
fid = fopen(fname_T1);
dum4 = fread(fid,1,'float32');
T1_count = fread(fid,100000000,'float64');
cumsum_T1 = cumsum(T1_count);

fname_T2 = sprintf('../data/T2_count.dat');
fid = fopen(fname_T2);
dum4 = fread(fid,1,'float32');
T2_count = fread(fid,100000000,'float64');
cumsum_T2 = cumsum(T2_count);

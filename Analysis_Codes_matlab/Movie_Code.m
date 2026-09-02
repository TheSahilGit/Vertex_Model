clear; clc; close all;

% ==================== options ====================
nrun = 1;
itList = 100000; %10000000;              % list of Fortran timesteps to render as frames
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

para2 = load("../para2_in.dat");
Lx = para2(1);
Ly = para2(2);

% Only load motility_store.dat if it's actually going to be used --
% it's a separate file read that every other colorBy option doesn't need.
etas = [];
if strcmp(colorBy, 'Motility')
    fid = fopen('../data/motility_store.dat');
    fread(fid, 1, 'float32');
    etas = fread(fid, Inf, 'float64');
    fclose(fid);
end

figure("Position", [200 200 800 800])

mov = VideoWriter(outFile);
mov.FrameRate = frameRate;
open(mov);

for it = itList

    [Lx, Ly, v, inn, num, forces, biochemdata] = LoadData(it, nrun);

    [colordata, colorbar_string] = ComputeCellColorData( ...
        colorBy, v, inn, num, forces, biochemdata, etas);

    TisuePlot(Lx, Ly, v, inn, num, colordata, colorbar_string, norm_flag, norm_range);

    title(num2str(it))
    F = getframe(gcf);
    writeVideo(mov, F);

    hold off;

end

close(mov)

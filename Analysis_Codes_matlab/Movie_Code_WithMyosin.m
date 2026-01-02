clear; clc ;
close all;


para2 = load("../para2_in.dat");
%etas_in = load("../motility_in.dat");



fname_etas = sprintf('../data/motility_store.dat');
fid = fopen(fname_etas);
dum4 = fread(fid,1,'float32');
etas = fread(fid,100000000,'float64');


%%

Lx = para2(1);
Ly = para2(2);
numdim  = para2(3);
vdim1 = para2(4);
vdim2 = para2(5);
inndim1 = para2(6);
inndim2 = para2(7);


mov = VideoWriter("Movie1_test.avi");
mov.FrameRate = 1;
open(mov); % Open the video file    before entering the loop

%fig = figure('Position', [100, 100, Lx*50, Ly*50]);
ct = 1;
%nt = length(borderver(:,1));
interv = 50000;

nrun = 2;

%filepath = '../';
boundary_ = [      1           2           3           4           5           6           7           8           9          10          11          12          13          14          15          16         241         242         243         244         245         246         247         248         249         250         251         252         253         254         255         256          16          32          48          64          80          96         112         128         144         160         176         192         208         224         240         256           1          17          33          49          65          81          97         113         129         145         161         177         193         209         225         241           1          16         241         256
];

for it = 1000000

    [Lx, Ly, v, inn, num, ~, biochemdata] = LoadData(it, nrun);
    Myosin = biochemdata(:,3);


    TisuePlot(Lx,Ly,v,inn,num,etas,boundary_, Myosin)

  

    title(num2str(it))
    F = getframe(gcf); % Get the frame
    writeVideo(mov, F); % Add frame to the video file

    hold off;

    ct = ct + 1;
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




function [cmX, cmY] = calculate_cellCentre(Lx,Ly,v,inn,num)
idx = find(inn(:,1) == 0, 1, 'first');
for i = 1:idx-1
    vx=v(inn(i,1:num(i)),1);
    vy=v(inn(i,1:num(i)),2);

    cmX(i) = mean(vx); cmY(i) = mean(vy);
end
end



function TisuePlot(Lx,Ly,v,inn,num,etas,boundary_, Myosin)
etas = etas*1000;
etasmax = max(etas);
etasmin = min(etas);



% plottill = Ly;
% for j = 1:Ly:(Lx*Ly)
%     for i = j:j+plottill-1

idx = find(num ~= 0, 1, 'last')

minMyosin = min(Myosin(1:idx));
maxMyosin = max(Myosin(1:idx));

% minMyosin = 0;
% maxMyosin = 1;


for i = 1:idx

        vx=v(inn(i,1:num(i)),1);
        vy=v(inn(i,1:num(i)),2);
        pl = polyshape(vx,vy);


        %r = (Myosin(i) - min(Myosin))/(max(Myosin) - min(Myosin));
        r = Myosin(i);
        
        cmap = jet(256);
        
        domain = linspace(minMyosin, maxMyosin, size(cmap,1));
        rgb = interp1(domain, cmap, r, 'linear', 'extrap');
        plot(pl, FaceColor=rgb, FaceAlpha=0.5, LineWidth=1.5)

        %plot(pl,EdgeColor='b',FaceColor='b', FaceAlpha=0.1, LineWidth=1.5)
        hold on;
    %end
end


cb = colorbar();
cb.Label.String = 'Myosin Intensity';

colormap(cmap)
clim([minMyosin maxMyosin]);

pbaspect([Lx/Ly 1 1])
axis([-8 Lx+8 -8 Ly+8])

%pbaspect([Lx/Ly 1 1])
axis("off")
set(gca, "FontName", "Serif", "FontSize", 30)
end

function [area,perimeter] = AverageAreaPerimeter(Lx,Ly,v,inn,num)

for i = 1:Lx*Ly
    % Get the vertices for the current polygon
    vx = v(inn(i, 1:num(i)), 1);
    vy = v(inn(i, 1:num(i)), 2);

    % Calculate area
    areaE(i) = polyarea(vx, vy);

    perimeterE(i) = calculatePerimeter(vx, vy);

end
area = areaE; perimeter = perimeterE;
%area = mean(areaE); perimeter = mean(perimeterE);
end

function perimeter = calculatePerimeter(vertices_x, vertices_y)
% Ensure that the number of x and y vertices match
assert(length(vertices_x) == length(vertices_y), 'Number of x and y vertices must be the same');

% Number of vertices
num_vertices = length(vertices_x);

% Initialize perimeter
perimeter = 0;

% Calculate distance between consecutive vertices
for i = 1:num_vertices
    % Calculate index of next vertex considering cyclic boundary
    next_index = mod(i, num_vertices) + 1;

    % Calculate distance between current and next vertices
    dx = vertices_x(next_index) - vertices_x(i);
    dy = vertices_y(next_index) - vertices_y(i);

    % Add distance to perimeter
    perimeter = perimeter + sqrt(dx^2 + dy^2);
end
end

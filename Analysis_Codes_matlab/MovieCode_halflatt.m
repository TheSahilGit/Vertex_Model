clear; clc ;
%close all;


para2 = load("../para2_in.dat");
etas = load("../motility_in.dat");



%%

Lx = para2(1);
Ly = para2(2);
numdim  = para2(3);
vdim1 = para2(4);
vdim2 = para2(5);
inndim1 = para2(6);
inndim2 = para2(7);


mov = VideoWriter("Movie1_test.avi");
mov.FrameRate = 24;
open(mov); % Open the video file    before entering the loop

fig = figure('Position', [100, 100, Lx*50, Ly*50]);
ct = 1;
%nt = length(borderver(:,1));
interv = 50000;

nrun = 1;

%filepath = '../';
boundary_ = [      1           2           3           4           5           6           7           8           9          10          11          12          13          14          15          16         241         242         243         244         245         246         247         248         249         250         251         252         253         254         255         256          16          32          48          64          80          96         112         128         144         160         176         192         208         224         240         256           1          17          33          49          65          81          97         113         129         145         161         177         193         209         225         241           1          16         241         256
];
for it = 232
    [Lx, Ly, v,inn,num, forces] = LoadData(it, nrun);

   
    f(ct) = sqrt(max(forces(:,1))^2 + max(forces(:,2))^2);      

    TisuePlot(Lx,Ly,v,inn,num,etas,boundary_)

    num(119)

    

    % hold on
    % vxx=v(inn(240,1:num(240)),1);
    % vyy=v(inn(240,1:num(240)),2);
    % pl = polyshape(vxx,vyy);
    % 
    % plot(pl, FaceColor='g')
    % 
    [cmX, cmY] = calculate_cellCentre(Lx,Ly,v,inn,num);

    globalMeanX = mean(cmX);
    globalmeanY = mean(cmY);

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
for i = 1:Lx*Ly
    vx=v(inn(i,1:num(i)),1);
    vy=v(inn(i,1:num(i)),2);

    cmX(i) = mean(vx); cmY(i) = mean(vy);
end
end



function TisuePlot(Lx,Ly,v,inn,num,etas,boundary_)
etas = etas*1000;
etasmax = max(etas);
etasmin = min(etas);

plottill = Ly;
for j = 1:Ly:(Lx*Ly)
    for i = j:j+plottill-1

        vx=v(inn(i,1:num(i)),1);
        vy=v(inn(i,1:num(i)),2);
        pl = polyshape(vx,vy);

     %   polyarea(vx,vy)

        % r = mean(etas(inn(i,1:num(i))));
        % if etasmax ~=0
        %     r = r/etasmax;
        % end
        % index = round(r * 99) + 1; % 100 colors in the colormap
        % % If index is out of bounds, set it to the first or last color index
        % index = max(min(index, 100), 1);
        % 
        % cmap = [flipud(jet(100)); 0 0 0];
        % % Interpolate the RGB value from the colormap
        % rgb = interp1(linspace(1, 0, size(cmap, 1)), cmap, r);
        %plot(pl, FaceColor=rgb, FaceAlpha=0.01, LineWidth=1.5)
        plot(pl,EdgeColor='b',FaceColor='r', FaceAlpha=0.01, LineWidth=1.5)
        hold on;
    end
end


%clim([0 1]); % Set the limits of the colorbar to the range of r
colormap([flipud(jet(100)); 0 0 0]); % Define colormap, including black for r = 0
cb = colorbar;
cb.Ticks = linspace(0, 1, 5);
cb.TickLabels = linspace(etasmax, etasmin, 5)./1000;
cb.Label.String = 'Motility';


hold on;


for ii=1:length(boundary_)
    i = boundary_(ii)
    vx=v(inn(i,1:num(i)),1);
    vy=v(inn(i,1:num(i)),2);
    pl = polyshape(vx,vy);
    plot(pl,FaceColor='g', FaceAlpha=1, LineWidth=1.5)
    hold on;
end


% for i=Ly
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='g', FaceAlpha=0.2, LineWidth=1.5)
%     hold on;
% end
%
% for i=Lx*Ly
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='g', FaceAlpha=0.2, LineWidth=1.5)
%     hold on;
% end
%
% for i=Lx*Ly-Ly+1
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='g', FaceAlpha=0.2, LineWidth=1.5)
%     hold on;
% end




% for i=65
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='g', FaceAlpha=0.2, LineWidth=1.5)
%     hold on;
% end
%
% for i=68
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='g', FaceAlpha=0.2, LineWidth=1.5)
%     hold on;
% end
%
% for i=71
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='g', FaceAlpha=0.2, LineWidth=1.5)
%     hold on;
% end
% for i=185
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='g', FaceAlpha=0.2, LineWidth=1.5)
%     hold on;
% end



% for ii = 1:borderver(im,1)
%     vbx(ii) = v(borderver(im ,ii+1),1);
%     vby(ii) = v(borderver(im ,ii+1),2);
% end
% scatter(vbx,vby, 'black', 'filled')


% for i=3
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='k', FaceAlpha=0.5, LineWidth=1.5)
%     hold on;
%     scatter(vx(2),vy(2),100, 'red')
% end
%
%
% for i=Lx*Ly-3*Ly+1
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='k', FaceAlpha=0.5, LineWidth=1.5)
%     hold on;
%     scatter(vx(4),vy(4),100, 'red')
% end
%
% for i=1
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='k', FaceAlpha=0.5, LineWidth=1.5)
%     hold on;
%     scatter(vx(1),vy(1), 'red', 'filled')
% end
% for i=2*Ly + 4
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='k', FaceAlpha=0.5, LineWidth=1.5)
%     hold on;
%     scatter(vx(1),vy(1), 'red', 'filled')
% end
% for i=Lx*Ly-2*Ly + 1
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='k', FaceAlpha=0.5, LineWidth=1.5)
%     hold on;
%     scatter(vx(1),vy(1), 'b', 'filled')
% end

% clim([0 1]); % Set the limits of the colorbar to the range of r
% colormap([flipud(jet(100)); 0 0 0]); % Define colormap, including black for r = 0
% cb = colorbar;
% cb.Ticks = linspace(0, 1, 5);
% cb.TickLabels = linspace(etasmax, etasmin, 5);
% cb.Label.String = 'Motility';

pbaspect([Lx/plottill 1 1])
%axis([-6 Lx+6 -6 plottill+6])
axis("off")

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







clear; clc; close all;


para2 = load("../para2_in.dat");
Lx = para2(1);
Ly = para2(2);
numdim  = para2(3);
vdim1 = para2(4);
vdim2 = para2(5);
inndim1 = para2(6);
inndim2 = para2(7);

para1 = readtable("../para1_in.dat");
dt = table2array(para1(8,1));


nrun = 2;




figure("Position",[200 200 800 800])

mov = VideoWriter("Movie_test.avi");
mov.FrameRate = 1;
open(mov); % Open the video file    before entering the loop


ct = 1; 
for it = 50000

    [Lx, Ly, v, inn, num, forces, ~, ~, stiched_data] = LoadData(it, nrun);

    
    %colordata = rand(1,numdim);
    %colorbar_string = "Nothing";

    Fcell =  MeanVertexForceMagnitude_Cell(forces(:,1:2), inn, num, numdim); 
    %Fcell =  MeanVertexForceMagnitude_Cell(forces(:,9:10), inn, num, numdim); 

    colordata = Fcell; 
    colorbar_string = "Force"; 


    norm_flag = 'data';
    norm_range = [];

    % norm_flag = '01';
    % norm_range = []; 
    % 
    % norm_flag = 'custom';
    % norm_range = [0 2]; 

    TisuePlot(Lx,Ly,v,inn,num, colordata, colorbar_string, norm_flag, norm_range)

    
    title(num2str(it))
    F = getframe(gcf); % Get the frame
    writeVideo(mov, F); % Add frame to the video file

    hold off;

    ct = ct + 1

end

close(mov)



function TisuePlot(Lx,Ly,v,inn,num,colordata,colorbar_string,...
                   norm_flag, norm_range)


Nc = find(num ~= 0, 1, 'last')
if isempty(Nc); return; end

% ----- normalization -----
switch norm_flag
    case 'data'
        cmin = min(colordata(1:Nc));
        cmax = max(colordata(1:Nc));
    case '01'
        cmin = 0;
        cmax = 1;
    case 'custom'
        cmin = norm_range(1);
        cmax = norm_range(2);
    otherwise
        error('Unknown norm_flag')
end

cmap   = jet(256);
domain = linspace(cmin, cmax, size(cmap,1));

%hold on
for i = 1:Nc
    vx = v(inn(i,1:num(i)),1);
    vy = v(inn(i,1:num(i)),2);

    pl = polyshape(vx,vy);

    rgb = interp1(domain, cmap, colordata(i), 'linear', 'extrap');

    plot(pl, ...
        FaceColor = rgb, ...
        FaceAlpha = 0.5, ...
        LineWidth = 1.5);

 % patch(vx, vy, rgb, ...
 %      'FaceAlpha', 0.5, ...
 %      'EdgeColor', 'k', ...
 %      'LineWidth', 1.0);

 hold on


end

cb = colorbar;
cb.Label.String = colorbar_string;
colormap(cmap)
clim([cmin cmax])

pbaspect([Lx/Ly 1 1])
axis([-8 Lx+8 -8 Ly+8])
axis off
set(gca,"FontName","Serif","FontSize",30)
set(gcf,"Renderer","opengl")
%hold off

end




function Fcell = MeanVertexForceMagnitude_Cell(forces, inn, num, numdim)

% forces(v,1) = Fx(v)
% forces(v,2) = Fy(v)

% 1) vertex-wise force magnitude
Fmag_v = sqrt(forces(:,1).^2 + forces(:,2).^2);

% 2) cell-wise averaging
Fcell = zeros(numdim,1);

Nc = find(num ~= 0, 1, 'last');

for i = 1:Nc
    vids = inn(i,1:num(i));      % vertices of cell i
    Fcell(i) = mean(Fmag_v(vids));
end

end






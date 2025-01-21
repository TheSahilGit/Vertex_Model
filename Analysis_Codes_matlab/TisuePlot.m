function TisuePlot(Lx,Ly,v,inn,num,etas)
etas = etas*1000;
etasmax = max(etas);
etasmin = min(etas);

plottill = Ly;
for j = 1:Ly:(Lx*Ly)
    for i = j:j+plottill-1

        vx=v(inn(i,1:num(i)),1);
        vy=v(inn(i,1:num(i)),2);
        pl = polyshape(vx,vy);
        r = mean(etas(inn(i,1:num(i))));
        if etasmax ~=0
            r = r/etasmax;
        end
        index = round(r * 99) + 1; % 100 colors in the colormap
        % If index is out of bounds, set it to the first or last color index
        index = max(min(index, 100), 1);

        cmap = [flipud(jet(100)); 0 0 0];
        % Interpolate the RGB value from the colormap
        rgb = interp1(linspace(1, 0, size(cmap, 1)), cmap, r);
        %plot(pl, FaceColor=rgb, FaceAlpha=0.01, LineWidth=1.5)
        plot(pl,EdgeColor='r' ,FaceColor='r', FaceAlpha=0.01, LineWidth=1.5)
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


% for i=1
%     vx=v(inn(i,1:num(i)),1);
%     vy=v(inn(i,1:num(i)),2);
%     pl = polyshape(vx,vy);
%     plot(pl,FaceColor='g', FaceAlpha=0.2, LineWidth=1.5)
%     hold on;
% end
%
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
axis([-6 Lx+6 -6 plottill+6])
axis("off")

end
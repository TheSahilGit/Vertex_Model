


plottill = Ly;
for j = 1:Ly:(Lx*Ly)
    for i = j:j+plottill-1

        vx=v(inn(i,1:num(i)),1);
        vy=v(inn(i,1:num(i)),2);
        pl = polyshape(vx,vy);

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

        
        plot(pl,EdgeColor='r',FaceColor='g', FaceAlpha=0.01, LineWidth=1.5)
        hold on;
        scatter(mean(vx), mean(vy), 5, 'MarkerFaceColor','g', 'MarkerEdgeColor','k')
    end
end

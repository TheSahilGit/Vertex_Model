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

     hold on;

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
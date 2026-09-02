function [colordata, colorbar_string] = ComputeCellColorData(colorBy, v, inn, num, forces, biochemdata, etas)
% COMPUTECELLCOLORDATA  Dispatch: compute a per-cell scalar field to color
% the tissue by, given a name. Generalizes what used to be three separate,
% nearly-identical scripts (Movie_Code.m: Force, MovieCode_halflatt.m:
% Motility, Movie_Code_WithMyosin.m: Myosin) into one selectable option.
%
%   [colordata, colorbar_string] = ComputeCellColorData(colorBy, v, inn, num, forces, biochemdata, etas)
%
% colorBy : one of 'Force' (default), 'Motility', 'Myosin', 'Rho', 'ROCK',
%           'Area', 'Perimeter', 'ShapeFactor', 'NumVertices'.
% forces      : vdim2 x 8 array from LoadData (fxx,fyy,fxx_ran,fyy_ran,
%               fxx_ABP,fyy_ABP,fxx_Polar,fyy_Polar), or [] if not needed.
% biochemdata : numdim x 3 array from LoadData (Rho, ROCK, Myosin), or []
%               if not needed.
% etas        : per-vertex motility field (from data/motility_store.dat),
%               or [] if not needed (only 'Motility' uses it).

Nc = find(num ~= 0, 1, 'last');

switch colorBy
    case 'Force'
        Fmag_v = sqrt(forces(:,1).^2 + forces(:,2).^2);
        colordata = perCellMeanOfVertexField(Fmag_v, inn, num, Nc);
        colorbar_string = 'Mean vertex force magnitude';

    case 'Motility'
        colordata = perCellMeanOfVertexField(etas, inn, num, Nc);
        colorbar_string = 'Motility';

    case 'Myosin'
        colordata = biochemdata(1:Nc, 3);
        colorbar_string = 'Myosin';

    case 'Rho'
        colordata = biochemdata(1:Nc, 1);
        colorbar_string = 'Rho';

    case 'ROCK'
        colordata = biochemdata(1:Nc, 2);
        colorbar_string = 'ROCK';

    case 'Area'
        colordata = perCellAreaPerimeter(v, inn, num, Nc, 'area');
        colorbar_string = 'Area';

    case 'Perimeter'
        colordata = perCellAreaPerimeter(v, inn, num, Nc, 'perimeter');
        colorbar_string = 'Perimeter';

    case 'ShapeFactor'
        area = perCellAreaPerimeter(v, inn, num, Nc, 'area');
        perim = perCellAreaPerimeter(v, inn, num, Nc, 'perimeter');
        colordata = perim ./ sqrt(area);
        colorbar_string = 'Shape factor (P/\surdA)';

    case 'NumVertices'
        colordata = double(num(1:Nc));
        colorbar_string = 'Number of vertices';

    otherwise
        error('ComputeCellColorData:unknownField', ...
            ['Unknown colorBy option "%s". Valid options: Force, Motility, ' ...
             'Myosin, Rho, ROCK, Area, Perimeter, ShapeFactor, NumVertices.'], colorBy);
end

end


function cellVal = perCellMeanOfVertexField(vertexField, inn, num, Nc)
cellVal = zeros(Nc, 1);
for i = 1:Nc
    vids = inn(i, 1:num(i));
    cellVal(i) = mean(vertexField(vids));
end
end


function val = perCellAreaPerimeter(v, inn, num, Nc, which)
val = zeros(Nc, 1);
for i = 1:Nc
    vx = v(inn(i, 1:num(i)), 1);
    vy = v(inn(i, 1:num(i)), 2);
    if strcmp(which, 'area')
        val(i) = polyarea(vx, vy);
    else
        dx = diff([vx; vx(1)]);
        dy = diff([vy; vy(1)]);
        val(i) = sum(sqrt(dx.^2 + dy.^2));
    end
end
end

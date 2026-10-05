function TisuePlot(Lx, Ly, v, inn, num, colordata, colorbar_string, norm_flag, norm_range, cmap_colors, renderer_mode)
% TISUEPLOT  Draw the whole tissue as a single colored patch object.
%
%   TisuePlot(Lx, Ly, v, inn, num, colordata, colorbar_string, norm_flag, norm_range)
%   TisuePlot(..., norm_range, cmap_colors)
%   TisuePlot(..., cmap_colors, renderer_mode)
%
% renderer_mode : optional, 'opengl' (default) or 'painters' (log.txt).
%                 opengl is hardware-accelerated, memory-stable and ~3x
%                 faster for the repeated clf+patch+getframe movie loop --
%                 measured directly on a 10,000-cell stress test: painters
%                 leaked ~490MB over 150 frames and crashed intermittently
%                 (getframe: "A valid figure or axes handle must be
%                 specified"), opengl stayed flat with no crashes, and both
%                 rendered visually identically (transparency/edges).
%                 painters is kept as an opt-in fallback (Movie_Code.m's
%                 rendererMode option) for cases opengl doesn't suit, e.g.
%                 a headless/software-only display with no GPU/driver.
%
% cmap_colors  : optional Nx3 RGB colormap (e.g. from slanCM.m via
%                GetFieldColormap.m) -- see log.txt. Omit/[] for the
%                original default (jet(256)); every existing caller that
%                doesn't pass this keeps its exact previous look. Whether
%                a field needs diverging treatment (e.g. Pressure/
%                ShearStress: a symmetric-about-zero clim, not raw
%                min/max) is decided by the caller (Movie_Code.m, via
%                GetFieldColormap.m) and folded into an ordinary
%                norm_flag/norm_range before it ever reaches here --
%                TisuePlot itself has no notion of "diverging".
%
% OPTIMIZATION: previously drew each cell as its own polyshape+plot() call
% in a loop -- one separate graphics object per cell. For a mesh of
% thousands of cells (this codebase's production runs use 32x32=1024+
% cells, more once cell division is on) that is one of the classic slow
% paths in MATLAB plotting: each plot() call has real per-object overhead,
% and building a movie means paying that cost every single frame. This now
% builds ONE NaN-padded Faces/Vertices matrix for the whole tissue and
% draws it with a single patch() call -- MATLAB's per-face coloring
% (FaceVertexCData + 'flat' + colormap/clim) replaces the old per-cell
% manual RGB interpolation, and rendering thousands of faces via one patch
% object is typically an order of magnitude faster than one plot() per
% cell. Output is visually identical to the old per-cell version.
%
% colordata        : per-cell scalar (or per-vertex, mean-pooled elsewhere)
%                     values to color each cell by, length >= Nc.
% colorbar_string  : label for the colorbar.
% norm_flag        : 'data' (use min/max of colordata(1:Nc)), '01', or
%                     'custom' (use norm_range).
% norm_range       : [cmin cmax], used only when norm_flag=='custom'.

Nc = find(num ~= 0, 1, 'last');
if isempty(Nc); return; end

if nargin < 10 || isempty(cmap_colors)
    cmap_colors = jet(256);
end
if nargin < 11 || isempty(renderer_mode)
    renderer_mode = 'opengl';
end

% ----- normalization -----
% 'data' only looks at LIVE cells (num(i)~=0) -- normally every cell in
% 1:Nc is live (Nc is defined as the last live index), so this is a
% no-op for ordinary callers. A caller that zeroes out num for a SUBSET
% of cells within 1:Nc (e.g. Movie_Code.m's plottill, to render/scale as
% if only that subset were the whole tissue) relies on this to exclude
% the hidden cells from the colorbar range too, not just from the drawn
% faces (which were already correctly skipped below -- an all-NaN face
% row for any num(i)==0).
switch norm_flag
    case 'data'
        live = num(1:Nc) ~= 0;
        cmin = min(colordata(live));
        cmax = max(colordata(live));
    case '01'
        cmin = 0;
        cmax = 1;
    case 'custom'
        cmin = norm_range(1);
        cmax = norm_range(2);
    otherwise
        error('TisuePlot:unknownNormFlag', 'Unknown norm_flag: %s', norm_flag);
end
if cmax == cmin
    % BUGFIX (log.txt): plain `eps` is eps(1) (~2.22e-16) regardless of
    % cmin's own magnitude -- for any |cmin| bigger than ~1, that's smaller
    % than a double's actual precision (ULP) at cmin, so `cmin + eps` was
    % silently rounding straight back down to cmin, leaving clim degenerate
    % anyway (clim() then errors: "Limits must be ... increasing"). eps(cmin)
    % scales the increment to cmin's own magnitude, always representable.
    cmax = cmin + eps(cmin);   % avoid a degenerate clim range
end

% ----- build one NaN-padded Faces matrix for the whole tissue -----
% PBC-aware (log.txt): a periodic mesh (Generate_Initial_Mesh.f90's
% if_periodic=true) has cells whose vertex list can legitimately span the
% periodic wrap -- a shared vertex ID with an opposite-edge neighbor cell,
% stored far away in raw (wrapped) coordinates. Indexing directly into
% the shared v array (as before) draws a long straight line connecting
% those raw, far-apart coordinates -- a "wired" look spanning the box.
% Instead, give each face its OWN local vertex list, each vertex unwrapped
% (minimum-image) relative to that face's own first vertex -- mirroring
% the Fortran-side Gather_Cell_Vertices_PBC technique -- so a cell that
% straddles the wrap is drawn as its true, compact shape. Harmless no-op
% for a non-periodic mesh: no real vertex pair within one cell is ever
% more than half the box apart, so the correction below always evaluates
% to zero there, and rendering is pixel-identical to before.
maxN = max(num(1:Nc));
F = NaN(Nc, maxN);
Vexp = zeros(sum(num(1:Nc)), 2);
row = 0;
for i = 1:Nc
    n = num(i);
    % A caller that zeroes out num for a subset of cells within 1:Nc
    % (Movie_Code.m's plottill) leaves a gap TisuePlot never otherwise
    % sees -- ordinarily Nc IS the last live index, so num(i) is never 0
    % for i<=Nc. F(i,:) is already all-NaN (patch draws nothing for it)
    % and no Vexp rows are allocated to it (sum(num(1:Nc)) already
    % excludes it) -- just skip the vertex math that would otherwise
    % index inn(i,1:0) with nothing in it.
    if n == 0; continue; end
    ids = inn(i, 1:n);
    x0 = v(ids(1), 1);
    y0 = v(ids(1), 2);
    for k = 1:n
        row = row + 1;
        dx = v(ids(k), 1) - x0;
        dx = dx - Lx * round(dx / Lx);
        dy = v(ids(k), 2) - y0;
        dy = dy - Ly * round(dy / Ly);
        Vexp(row, 1) = x0 + dx;
        Vexp(row, 2) = y0 + dy;
        F(i, k) = row;
    end
end

patch('Faces', F, 'Vertices', Vexp, ...
    'FaceVertexCData', colordata(1:Nc), ...
    'FaceColor', 'flat', ...
    'FaceAlpha', 0.5, ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.5);

cb = colorbar;
cb.Label.String = colorbar_string;
colormap(cmap_colors);
clim([cmin cmax]);

pbaspect([Lx/Ly 1 1])
axis([-4 Lx+4 -4 Ly+4])
axis off
set(gca, "FontName", "Serif", "FontSize", 30)
set(gcf, "Renderer", renderer_mode)

end

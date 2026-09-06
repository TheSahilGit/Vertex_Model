function [cmap, is_diverging] = GetFieldColormap(colorBy)
% GETFIELDCOLORMAP  Central place to pick which colormap (via slanCM.m)
% each Movie_Code.m colorBy quantity is rendered with. Movie_Code.m/
% TisuePlot.m don't hardcode any of this -- they just ask this function
% for "the right colormap for X" -- so to change a field's colormap later,
% edit ONE line in the COLORMAP_TABLE/DIVERGING_FIELDS block below; every
% name slanCM.m understands is fair game (see its own header comment for
% the full list pulled from matplotlib/scicomap/cmasher/colorcet).
%
%   [cmap, is_diverging] = GetFieldColormap(colorBy)
%
% cmap         : 256 x 3 RGB array, ready for TisuePlot.m's cmap_colors.
% is_diverging : true for fields with a signed, physically meaningful
%                zero (Pressure, ShearStress) -- passed through to
%                TisuePlot.m so norm_flag=='data' centers the colorbar on
%                that zero instead of the field's raw (usually
%                asymmetric) min/max, which would otherwise misalign a
%                diverging colormap's neutral middle color with the
%                actual sign-change point.
%
% Scientific-rigor convention (log.txt, user request): DIVERGING colormaps
% only for genuinely signed/zero-centered quantities; SEQUENTIAL for every
% nonnegative/magnitude-like quantity. Every colormap picked below is a
% simple, established, publication-standard choice (matplotlib's own
% perceptually-uniform defaults, ColorBrewer, or literature convention for
% FTLE specifically) -- deliberately NOT rainbow/jet (this codebase's old
% default), which is well known to be non-perceptually-uniform, visually
% distort data, and fail on greyscale printing.

persistent COLORMAP_TABLE DIVERGING_FIELDS
if isempty(COLORMAP_TABLE)
    COLORMAP_TABLE = containers.Map();

    % ---- sequential (nonnegative / magnitude-like quantities) ----
    COLORMAP_TABLE('Force')       = 'viridis';   % perceptually-uniform default
    COLORMAP_TABLE('Motility')    = 'plasma';
    COLORMAP_TABLE('Myosin')      = 'inferno';   % warm ramp -- contractile activity
    COLORMAP_TABLE('Rho')         = 'magma';
    COLORMAP_TABLE('ROCK')        = 'cividis';   % colorblind-safe, distinct from Rho
    COLORMAP_TABLE('Area')        = 'Blues';     % ColorBrewer, simple geometric field
    COLORMAP_TABLE('Perimeter')   = 'Greens';
    COLORMAP_TABLE('ShapeFactor') = 'YlOrBr';    % warm ramp, common in shape-index figures
    COLORMAP_TABLE('NumVertices') = 'PuBuGn';    % discrete count, visually distinct
    COLORMAP_TABLE('FTLE')        = 'YlOrRd';    % literature convention for FTLE fields

    % ---- diverging (signed, zero-centered quantities) ----
    COLORMAP_TABLE('Pressure')    = 'RdBu';      % classic compression/tension diverging map
    COLORMAP_TABLE('ShearStress') = 'PuOr';      % diverging, visually distinct from Pressure

    DIVERGING_FIELDS = {'Pressure', 'ShearStress'};
end

if isKey(COLORMAP_TABLE, colorBy)
    name = COLORMAP_TABLE(colorBy);
else
    name = 'viridis';   % fallback for any future/unlisted colorBy
end
cmap = slanCM(name, 256);

is_diverging = any(strcmp(colorBy, DIVERGING_FIELDS));

end

function PlotAnalysis(varargin)
% PLOTANALYSIS  One figure, one subplot per requested diagnostic, each
% independently toggleable -- the single entry point tying together
% compute_Circularity, compute_Qt, compute_MSD_cellID,
% compute_StressTensor_series, and the global Energy/ShearStress/cumsumT1/
% cumsumT2 series (LoadGlobalTimeSeries).
%
% USAGE (all arguments are optional name-value pairs):
%
%   PlotAnalysis()                                  % Energy only, defaults
%   PlotAnalysis('doEnergy', true, 'doCircularity', true, 'doQt', true, ...
%                'doMSD', true, 'doShearStress', true, 'doPressure', true)
%   PlotAnalysis('doCumsumT1', true, 'doCumsumT2', true)   % only once the run has finished
%   PlotAnalysis('doCircularity', true, 'onlyPlot', true, 'fontSize', 20)  % replot from cache, restyled
%
% Name-value options:
%   nrun         (1)      1 or 2, matching LoadData's convention.
%   itStart      (1)      first Fortran timestep to include for the
%                         per-frame panels (Circularity/Qt/MSD/Pressure).
%   itEnd        ([])     last timestep to include for those panels.
%                         [] = auto-detect the latest snapshot actually on
%                         disk (FindLatestAvailableIt) -- safe to leave
%                         empty for a still-running simulation.
%   itInterval   ([])     sampling stride, in units of it (must be a
%                         multiple of it_dump; rounded up otherwise).
%                         [] = use it_dump itself (every available dump).
%   doEnergy         (true)   Energy vs. time            (global file, live)
%   doShearStress    (false)  whole-tissue ShearStress vs. time (global file, live)
%   doPressure       (false)  locally-recomputed Pressure vs. time (always
%                             recomputed -- Pressure is never stored on disk)
%   doCircularity    (false)  tissue-boundary circularity vs. time
%   doQt             (false)  self-overlap order parameter Q(t)
%   doMSD            (false)  mean-squared displacement, tracked by
%                             cell_identity (robust to T2 cell removal)
%   doCumsumT1       (false)  cumulative T1 count -- ONLY plotted once the
%                             run has reached it==totT (see below); a
%                             partial cumulative count mid-run would be
%                             misleading, not just incomplete.
%   doCumsumT2       (false)  cumulative T2 count -- same restriction.
%   ac           (1.0)    Qt cage/persistence length scale.
%   radius       (10)     radius (from tissue COM) used by the Pressure
%                         panel's local stress-tensor recomputation.
%
% ---- re-plotting without re-running the analysis ----
%   onlyPlot     (false)  if true, skip every compute_*.m call entirely and
%                         instead load each requested panel's data from
%                         cacheFile (whatever was cached the last time it
%                         was computed with onlyPlot=false). Use this to
%                         iterate on styling (fontSize/lineWidth/
%                         markerSize below) without repaying the cost of
%                         reloading mesh snapshots. A panel that was never
%                         computed/cached yet is skipped with a warning
%                         (run once with onlyPlot=false first).
%   cacheFile    ('PlotAnalysis_cache.mat')
%                         where computed panel data is cached. Every
%                         non-onlyPlot run MERGES its freshly computed
%                         panels into whatever is already in this file, so
%                         you can build up a full cache across several
%                         separate runs with different flag combinations,
%                         then later do one onlyPlot=true run with
%                         everything enabled to replot it all together.
%
% ---- plot styling ----
%   fontSize     (14)     axes/tick/label font size, applied to every panel.
%   lineWidth    (2)      line width, applied to every panel.
%   markerSize   (6)      marker size, for the panels that use markers
%                         (Pressure/Circularity/Qt/MSD -- Energy/
%                         ShearStress/cumsumT1/cumsumT2 are plain lines,
%                         one value per Fortran timestep, too dense for
%                         markers to be useful).
%   showTitles   (true)   set false to draw every panel without its title
%                         (axis labels are unaffected either way).
%
% ---- axis scales (fixed per panel, not user-configurable -- these match
% the physics each quantity is conventionally plotted with) ----
%   Energy: semilogy.  cumsumT1/cumsumT2: loglog.  Q(t): semilogx.
%   MSD: loglog, with t^1/t^2 reference guide lines overlaid (anchored to
%   pass through the actual data's last point, not a hardcoded prefactor,
%   so they stay meaningful regardless of which run this is pointed at).
%   ShearStress/Pressure/Circularity: linear.
%
% PROGRESS: prints "[k/N] Running <panel>..." as each panel starts
% computing (or "Loading cached <panel>..." under onlyPlot), and for the
% expensive per-frame panels (Pressure/Circularity/Qt/MSD, each of which
% reloads a full mesh snapshot per requested it) an in-place-updating
% "<panel>: XX% (k/n)" line as frames are processed.
%
% WHY cumsumT1/cumsumT2 are final-only: Total_T1_count/Total_T2_count
% (T1_count.dat/T2_count.dat) are running counts filled in as the
% simulation proceeds; their CUMULATIVE SUM is only a meaningful "total T1
% events so far" once every entry up to it==totT has actually been written.
% Requesting them before that would silently plot a truncated, misleading
% curve -- so this function checks length(T1_count) (or T2_count) against
% totT (from para1_in.dat) and skips the panel with a warning if the run
% isn't finished yet, rather than plotting something wrong. (This check is
% skipped under onlyPlot -- a cached cumsum panel was already validated as
% final when it was first computed.)

p = inputParser;
addParameter(p, 'nrun', 1);
addParameter(p, 'itStart', 1);
addParameter(p, 'itEnd', []);
addParameter(p, 'itInterval', []);
addParameter(p, 'doEnergy', true);
addParameter(p, 'doShearStress', false);
addParameter(p, 'doPressure', false);
addParameter(p, 'doCircularity', false);
addParameter(p, 'doQt', false);
addParameter(p, 'doMSD', false);
addParameter(p, 'doCumsumT1', false);
addParameter(p, 'doCumsumT2', false);
addParameter(p, 'ac', 1.0);
addParameter(p, 'radius', 10);
addParameter(p, 'onlyPlot', false);
addParameter(p, 'cacheFile', 'PlotAnalysis_cache.mat');
addParameter(p, 'fontSize', 14);
addParameter(p, 'lineWidth', 2);
addParameter(p, 'markerSize', 6);
addParameter(p, 'showTitles', true);
parse(p, varargin{:});
opt = p.Results;

style = struct('fontSize', opt.fontSize, 'lineWidth', opt.lineWidth, 'markerSize', opt.markerSize);

p1 = ReadPara1Params("../para1_in.dat");
it_dump = p1.it_dump;
totT    = p1.totT;

% itList is only needed to actually compute the per-frame panels -- skip
% deriving it (and the disk probing FindLatestAvailableIt does) entirely
% under onlyPlot, since nothing here will call compute_*.m.
itList = [];
if ~opt.onlyPlot
    if isempty(opt.itEnd)
        opt.itEnd = FindLatestAvailableIt(opt.nrun, it_dump, totT);
        if opt.itEnd == 0
            error('PlotAnalysis:noData', 'No snapshot data found on disk yet for nrun=%d.', opt.nrun);
        end
    end
    itList = GetSnapshotItList(it_dump, opt.itStart, opt.itEnd, opt.itInterval);
end

% ---- load whatever's already cached, so a fresh compute of some panels
% doesn't clobber previously cached OTHER panels ----
cache = struct();
if isfile(opt.cacheFile)
    s = load(opt.cacheFile, 'cache');
    cache = s.cache;
end

% ---- Build the list of requested panels, in a fixed, sensible order.
% Each entry pairs a cache key + title with a computeFcn (data source) and
% a plotFcn (how to draw it) -- computing and drawing are deliberately
% kept separate so onlyPlot can skip straight to plotFcn using cached data. ----
requested = {};

if opt.doEnergy
    requested{end+1} = reqPanel('Energy', 'Energy vs. time', @() computeEnergy(), ...
        @(ax,t,y,s) plotSeries(ax,t,y,s,'Energy','semilogy',false), false, '');
end
if opt.doShearStress
    requested{end+1} = reqPanel('ShearStress', 'ShearStress vs. time', @() computeShearStress(), ...
        @(ax,t,y,s) plotSeries(ax,t,y,s,'Shear Stress','linear',false), false, '');
end
if opt.doPressure
    requested{end+1} = reqPanel('Pressure', 'Pressure vs. time', @() computePressure(itList, opt), ...
        @(ax,t,y,s) plotSeries(ax,t,y,s,'Pressure','linear',true), false, '');
end
if opt.doCircularity
    requested{end+1} = reqPanel('Circularity', 'Circularity vs. time', @() computeCircularity(itList, opt), ...
        @plotCircularityAx, false, '');
end
if opt.doQt
    requested{end+1} = reqPanel('Qt', 'Q(t)', @() computeQt(itList, opt), ...
        @(ax,t,y,s) plotSeries(ax,t,y,s,'Q(t)','semilogx',true), false, '');
end
if opt.doMSD
    requested{end+1} = reqPanel('MSD', 'MSD (cell-ID tracked)', @() computeMSD(itList, opt), ...
        @plotMSDAx, false, '');
end
if opt.doCumsumT1
    requested{end+1} = reqPanel('CumsumT1', 'Cumulative T1 count', @() computeCumsum('T1'), ...
        @(ax,t,y,s) plotSeries(ax,t,y,s,'Cumulative T1 count','loglog',false), true, 'T1');
end
if opt.doCumsumT2
    requested{end+1} = reqPanel('CumsumT2', 'Cumulative T2 count', @() computeCumsum('T2'), ...
        @(ax,t,y,s) plotSeries(ax,t,y,s,'Cumulative T2 count','loglog',false), true, 'T2');
end

% ---- Resolve each requested panel to (time, y) data, either by computing
% it fresh (and caching the result) or, under onlyPlot, by loading it from
% the cache -- then hand off to layout/drawing. ----
panels = {};   % each entry: {title, time, y, plotFcn}
nRequested = numel(requested);
cacheDirty = false;

for i = 1:nRequested
    r = requested{i};

    if r.finalOnly && ~opt.onlyPlot && ~checkFinal(totT, r.which)
        continue;   % checkFinal already warned
    end

    if opt.onlyPlot
        if isfield(cache, r.key)
            fprintf('[%d/%d] Loading cached %s...\n', i, nRequested, r.title);
            time = cache.(r.key).time;
            y = cache.(r.key).y;
        else
            warning('PlotAnalysis:notCached', ...
                ['No cached data for "%s" -- skipping. Run once with ' ...
                 'onlyPlot=false to compute and cache it first.'], r.title);
            continue;
        end
    else
        fprintf('[%d/%d] Running %s...\n', i, nRequested, r.title);
        [time, y] = r.computeFcn();
        cache.(r.key) = struct('time', time, 'y', y);
        cacheDirty = true;
    end

    panels{end+1} = {r.title, time, y, r.plotFcn}; %#ok<AGROW>
end

if cacheDirty
    save(opt.cacheFile, 'cache');
end

if isempty(panels)
    warning('PlotAnalysis:nothingToPlot', 'No panels requested (or all requested panels were skipped) -- nothing to plot.');
    return;
end

% ---- Lay out and draw ----
figure('Position', [100 100 500*min(numel(panels),3) 450*ceil(numel(panels)/3)]);
set(gcf, 'Renderer', 'Painters');   % crisp vector output for line/marker plots
tl = tiledlayout(ceil(numel(panels)/3), min(numel(panels),3), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:numel(panels)
    [ttl, time, y, plotFcn] = panels{k}{:};
    ax = nexttile(tl);
    plotFcn(ax, time, y, style);
    if opt.showTitles
        title(ax, ttl);
    end
    set(ax, 'FontSize', style.fontSize);
    axis(ax, 'square');
end
fprintf('PlotAnalysis: all %d panel(s) drawn.\n', numel(panels));

end


function r = reqPanel(key, title, computeFcn, plotFcn, finalOnly, which)
r = struct('key', key, 'title', title, 'computeFcn', computeFcn, ...
    'plotFcn', plotFcn, 'finalOnly', finalOnly, 'which', which);
end


% ==================== per-panel data sources (compute_*.m callers) ====================

function [time, Energy] = computeEnergy()
[time, Energy] = LoadGlobalTimeSeries();
end

function [time, ShearStress] = computeShearStress()
[time, ~, ShearStress] = LoadGlobalTimeSeries();
end

function [time, Pressure] = computePressure(itList, opt)
progressFcn = makeProgressPrinter('Pressure');
[time, ~, Pressure] = compute_StressTensor_series(itList, opt.nrun, opt.radius, progressFcn);
end

function [time, circularity] = computeCircularity(itList, opt)
progressFcn = makeProgressPrinter('Circularity');
[time, circularity] = compute_Circularity(itList, opt.nrun, progressFcn);
end

function [time, Qt] = computeQt(itList, opt)
progressFcn = makeProgressPrinter('Q(t)');
[time, Qt] = compute_Qt(itList, opt.nrun, opt.ac, progressFcn);
end

function [time, MSD] = computeMSD(itList, opt)
progressFcn = makeProgressPrinter('MSD');
[time, MSD] = compute_MSD_cellID(itList, opt.nrun, progressFcn);
end

function [time, y] = computeCumsum(which)
[time, ~, ~, ~, ~, cumsum_T1, cumsum_T2] = LoadGlobalTimeSeries();
if strcmp(which, 'T1')
    y = cumsum_T1;
else
    y = cumsum_T2;
end
end


% ==================== per-panel plotting (styling applied here) ====================

function plotFn = axisPlotFcn(axtype)
switch axtype
    case 'semilogy'
        plotFn = @semilogy;
    case 'semilogx'
        plotFn = @semilogx;
    case 'loglog'
        plotFn = @loglog;
    otherwise
        plotFn = @plot;
end
end

function plotSeries(ax, time, y, style, ylab, axtype, useMarker)
plotFn = axisPlotFcn(axtype);
if useMarker
    plotFn(ax, time, y, '-o', 'LineWidth', style.lineWidth, 'MarkerSize', style.markerSize);
else
    plotFn(ax, time, y, 'LineWidth', style.lineWidth);
end
xlabel(ax, 'Time'); ylabel(ax, ylab);
end

function plotCircularityAx(ax, time, y, style)
plot(ax, time, y, '-o', 'LineWidth', style.lineWidth, 'MarkerSize', style.markerSize);
yline(ax, 1, '--');
xlabel(ax, 'Time'); ylabel(ax, 'Circularity');
end

function plotMSDAx(ax, time, y, style)
loglog(ax, time, y, '-o', 'LineWidth', style.lineWidth, 'MarkerSize', style.markerSize, ...
    'DisplayName', 'MSD');
hold(ax, 'on');

% t^1 (diffusive) / t^2 (ballistic) reference guides. Anchored to pass
% through the LAST valid (finite, positive) data point rather than a
% hardcoded prefactor (the old AnalysisMSD.m script used 5e-4, tuned for
% one specific dataset's magnitude/timescale and meaningless for any
% other) -- this way the guides land on the actual curve regardless of
% which run/parameters this is pointed at.
%
% BUGFIX: time(1)/y(1) is the reference frame compared to itself, so
% y(1)==0 exactly by construction (compute_MSD_cellID.m) -- loglog can't
% plot a zero, so it silently OMITS that point, meaning the real MSD
% trace visually starts at time(2), not time(1). The guide-line formula
% below is never exactly zero, so drawing it over the full unfiltered
% time array made it visibly extend further left (to time(1)) than the
% real curve. Fix: draw the guides over the same valid (positive,
% finite) domain the real curve actually gets plotted over, `tRef`, not
% the raw `time`.
valid = isfinite(time) & isfinite(y) & time > 0 & y > 0;
if any(valid)
    tRef = time(valid); yRef = y(valid);
    %t0 = tRef(end); y0 = yRef(end);
    t0 = tRef(1); y0 = yRef(1);

    guideWidth = max(1, style.lineWidth - 1);
    loglog(ax, tRef, y0 * (tRef / t0), '--', 'LineWidth', guideWidth, 'DisplayName', 't^1');
    loglog(ax, tRef, y0 * (tRef / t0).^2, ':', 'LineWidth', guideWidth, 'DisplayName', 't^2');
end

hold(ax, 'off');
xlabel(ax, 'Time'); ylabel(ax, 'MSD');
legend(ax, 'Location', 'best');
end


% ==================== progress printing ====================

function fcn = makeProgressPrinter(label)
% Returns a function handle progressFcn(k, n) that prints an
% in-place-updating "label: XX% (k/n)" line, ending with a newline once
% k==n. Passed down into the compute_*.m functions' per-frame loops so
% the actual print statements stay centralized here in PlotAnalysis.m
% rather than duplicated in every compute_ function.
fcn = @(k, n) printProgress(label, k, n);
end

function printProgress(label, k, n)
fprintf('\r  %s: %3d%% (%d/%d)', label, round(100*k/n), k, n);
if k == n
    fprintf('\n');
end
end


function ok = checkFinal(totT, which)
[~, ~, ~, T1_count, T2_count] = LoadGlobalTimeSeries();
if strcmp(which, 'T1')
    n = numel(T1_count);
else
    n = numel(T2_count);
end
ok = (n >= totT);
if ~ok
    warning('PlotAnalysis:notFinal', ...
        ['Cumulative %s count requires a completed run (has %d of %d ' ...
         'timesteps'' worth of data) -- skipping this panel. A partial ' ...
         'cumulative count would be misleading, not just incomplete.'], ...
        which, n, totT);
end
end

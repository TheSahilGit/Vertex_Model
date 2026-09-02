clear; clc; close all;
% Thin driver script for PlotAnalysis.m -- edit the options below and run
% this file. PlotAnalysis itself is a function (see PlotAnalysis.m); this
% script just calls it with whatever you set here, the same pattern
% Movie_Code.m uses for its own "options" block.

% ==================== options ====================
nrun = 1;

itStart    = 100;          % first Fortran timestep to include
itEnd      = [];         % last timestep to include; [] = auto-detect the
                          % latest snapshot actually on disk (safe to leave
                          % empty for a still-running simulation)
itInterval = 200000;     % sampling stride, in units of it -- does NOT have
                          % to be it_dump; use a bigger number for a
                          % faster/coarser check, smaller for more detail.
                          % Must be a multiple of it_dump (only multiples of
                          % it_dump -- plus it=1,2 -- ever have a dump file
                          % on disk); anything else is rounded up
                          % automatically, with a warning.

% ---- which panels to draw ----
doEnergy       = true;
doShearStress  = false;
doPressure     = false;
doCircularity  = true;
doQt           = true;
doMSD          = true;
doCumsumT1     = true;   % only actually plotted once the run has reached
doCumsumT2     = true;   % it==totT -- see PlotAnalysis.m's header comment

% ---- extra parameters used by specific panels ----
ac     = 0.8;   % Qt cage/persistence length scale
radius = 10;    % region (from tissue COM) used by the Pressure panel

% ---- re-plotting without re-running the analysis ----
onlyPlot  = true;                        % true = skip all computation, just
                                           % replot from cacheFile (e.g. to
                                           % try different styling below)
cacheFile = 'PlotAnalysis_cache.mat';     % every non-onlyPlot run merges its
                                           % freshly computed panels into
                                           % whatever's already in this file

% ---- plot styling ----
fontSize   = 40;
lineWidth  = 6;
markerSize = 10;
showTitles = false;   % set false to draw every panel without its title
% ===================================================

PlotAnalysis('nrun', nrun, 'itStart', itStart, 'itEnd', itEnd, 'itInterval', itInterval, ...
    'doEnergy', doEnergy, 'doShearStress', doShearStress, 'doPressure', doPressure, ...
    'doCircularity', doCircularity, 'doQt', doQt, 'doMSD', doMSD, ...
    'doCumsumT1', doCumsumT1, 'doCumsumT2', doCumsumT2, ...
    'ac', ac, 'radius', radius, ...
    'onlyPlot', onlyPlot, 'cacheFile', cacheFile, 'showTitles', showTitles, ...
    'fontSize', fontSize, 'lineWidth', lineWidth, 'markerSize', markerSize);

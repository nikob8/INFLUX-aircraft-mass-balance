%% Collect_Solutions.m
%  Gathers <Date>/Solution-MethodA.mat ... Solution-MethodD.mat from every
%  flight-date folder into ONE file:
%
%     All_Solutions.mat   with structs  SolA, SolB, SolC, SolD  (one element per
%                         flight day, all saved variables except FlightList)
%     All_Solutions.csv   one row per transect: date, transect, time, height,
%                         ABL, CosTheta and the flux of each method
%
%  Run it from the folder that holds the date folders (your MassB folder) after
%  the four Method scripts have finished.  Missing files are skipped and listed.

clear; clc;
ROOT = pwd;                                   % folder with 20110301/, 20110429/, ...
d = dir(ROOT);
dates = sort({d([d.isdir] & ~cellfun(@isempty, regexp({d.name}, '^\d{8}$', 'once'))).name});
fprintf('%d date folders found in %s\n', numel(dates), ROOT);

methods = {'A','B','C','D'};
for m = 1:numel(methods)
    S = struct([]); missing = {};
    for k = 1:numel(dates)
        f = fullfile(ROOT, dates{k}, ['Solution-Method' methods{m} '.mat']);
        if ~isfile(f), missing{end+1} = dates{k}; continue; end   %#ok<SAGROW>
        x = load(f);
        if isfield(x,'FlightList'), x = rmfield(x,'FlightList'); end
        x.Date = dates{k};
        if isempty(S), S = x; else, S = merge_struct(S, x); end
    end
    eval(['Sol' methods{m} ' = S;']);
    fprintf('Method %s: %d days collected', methods{m}, numel(S));
    if ~isempty(missing), fprintf(', missing: %s', strjoin(missing, ' ')); end
    fprintf('\n');
end
save(fullfile(ROOT,'All_Solutions.mat'), 'SolA', 'SolB', 'SolC', 'SolD');

%% one table, one row per downwind transect
Date = {}; Transect = []; TimeUTC = {}; Height_m = []; ABL_m = []; CosTheta = [];
FluxA = []; FluxB = []; FluxC = []; FluxD = []; WSPD = []; WDIR = [];
for i = 1:numel(SolA)
    a = SolA(i); n = numel(a.Flux);
    Date     = [Date; repmat({a.Date}, n, 1)];                          %#ok<AGROW>
    Transect = [Transect; (1:n)'];                                       %#ok<AGROW>
    t = a.DownwindTime(:);
    TimeUTC  = [TimeUTC; arrayfun(@(s) sprintf('%02d:%02d', floor(s/3600), round(mod(s,3600)/60)), t, 'UniformOutput', false)]; %#ok<AGROW>
    Height_m = [Height_m; a.LegHeight(:)];                               %#ok<AGROW>
    ABL_m    = [ABL_m; a.DownwindABL(:)];                                %#ok<AGROW>
    CosTheta = [CosTheta; getf(a,'CosTheta',n)];                         %#ok<AGROW>
    WSPD     = [WSPD; repmat(a.MeanWSPD, n, 1)];                         %#ok<AGROW>
    WDIR     = [WDIR; repmat(a.MeanWDIR, n, 1)];                         %#ok<AGROW>
    FluxA    = [FluxA; a.Flux(:)];                                       %#ok<AGROW>
    FluxB    = [FluxB; fluxOf(SolB, a.Date, n)];                         %#ok<AGROW>
    FluxC    = [FluxC; fluxOf(SolC, a.Date, n)];                         %#ok<AGROW>
    FluxD    = [FluxD; fluxOf(SolD, a.Date, n)];                         %#ok<AGROW>
end
T = table(Date, Transect, TimeUTC, Height_m, ABL_m, WSPD, WDIR, CosTheta, FluxA, FluxB, FluxC, FluxD);
writetable(T, fullfile(ROOT,'All_Solutions.csv'));
fprintf('Saved All_Solutions.mat and All_Solutions.csv (%d transects)\n', height(T));
fprintf('CosTheta: min %.3f, median %.3f, max %.3f\n', min(CosTheta), median(CosTheta,'omitnan'), max(CosTheta));

%% ------------------------------------------------------------ helpers
function S = merge_struct(S, x)
% append x to struct array S even if the field sets differ
    fs = fieldnames(S); fx = fieldnames(x);
    for q = setdiff(fx, fs)', [S.(q{1})] = deal([]); end
    for q = setdiff(fs, fx)', x.(q{1}) = []; end
    S(end+1) = orderfields(x, S(1));
end

function v = getf(a, name, n)
    if isfield(a, name) && numel(a.(name)) == n, v = a.(name)(:); else, v = NaN(n,1); end
end

function v = fluxOf(Sol, date, n)
    v = NaN(n,1);
    if isempty(Sol), return; end
    k = find(strcmp({Sol.Date}, date), 1);
    if ~isempty(k) && numel(Sol(k).Flux) == n, v = Sol(k).Flux(:); end
end

%% Diagnostic_ratio.m
%  WRF-Chem diagnostic ratio of each flight day (Section 3.2 of the paper): the
%  mean CH4 enhancement of the Marion County (city) tracers divided by the
%  standard deviation of the regional background tracers along an idealized
%  downwind transect, at model level 8 (about 300 m above ground), averaged
%  over 17 to 20 UTC.  Days with a ratio above 2 are the "good" background days
%  of Table 3 and Figures 13 and 14.
%
%  The idealized transects are straight lines of grid points on the 3 km grid
%  about 30 km upwind and downwind of the Marion County centre, chosen by the
%  mean model wind direction of the day (eight sectors; the lines are drawn in
%  Figure 7).  They are defined by grid indices inside the 35 x 36 point box of
%  the original analysis (full-grid indices 143:177 west to east and 132:167
%  south to north), so that the published values are reproduced.  Both the
%  wind used for the sector and the lines are grid-relative.
%
%  Input: the WRF_subset_<yyyymmdd>.mat files written by Extract_WRF_subset.m
%  (released as a Zenodo data set; they hold the tracers at 16 to 20 UTC) and,
%  optionally, All_Solutions.mat for the Downwind Edge emission rates.
%  Output: diagnostic_ratios.csv (one row per flight day) and Diag_ratio.mat
%  with the matrix Diag_All2 [wind speed, wind direction, ratio, Downwind Edge
%  daily emission rate, yyyymmdd], the layout read by plot_Figure13_diagnostic.m,
%  plot_Figure14_bootstrap.m and Trajectory_check_R2_7.m (point their
%  DIAG_FILE at Diag_ratio.mat).  If DIAG_OLD exists, the script prints the
%  published values next to the new ones.
%
%  This script replaces the original diagnostic-ratio script (N. Balashov),
%  which read the wrfout files directly; the computation is the same,
%  with one correction: for south-easterly winds the original concatenated the
%  city values into the background vector (a typo), so the ratio of 11/21/2014
%  is recomputed here.

clear; clc;

%% ------------------------------------------------------------ settings
SUB_DIR  = fullfile(pwd, 'WRF_subsets');                    % WRF_subset_<yyyymmdd>.mat files
SOL_FILE = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB/All_Solutions.mat';   % optional
DIAG_OLD = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/Diag4.mat';   % optional: published values for comparison
OUT_DIR  = pwd;
LEVEL    = 8;            % model level (about 300 m above ground)
HOURS    = 16:20;        % UTC hours read; the ratio uses the downwind line at 17 to 20 UTC
BOX_I    = 143:177;      % analysis box in full d02 grid indices (west to east)
BOX_J    = 132:167;      %                                        (south to north)
THRESH   = 2;

%% ------------------------------------------------------------ idealized transects (box indices [x y])
% One pair of lines per wind sector: 'dn' is the downwind line (city plume
% and background are evaluated on it), 'up' the upwind line (background
% only, one hour earlier).  N-S and E-W lines are grid rows/columns, the
% diagonal sectors use 45 degree lines through the corners of the box.
L = struct();
L.S  = tlines(6:28, 8,            6:28, 28);
L.N  = tlines(6:28, 28,           6:28, 8);
L.W  = tlines(7, 6:32,            28, 6:32);
L.E  = tlines(31, 6:33,           7, 6:32);
L.SW = tlines(16:-1:3, 3:16,      33:-1:18, 18:33);
L.NE = tlines(33:-1:19, 18:32,    16:-1:3, 3:16);
L.NW = tlines(3:17, 18:32,        19:33, 2:16);
L.SE = tlines(19:33, 4:18,        2:17, 18:33);
sectorNames = {'N','NE','E','SE','S','SW','W','NW'};         % centred on 0, 45, ... 315 degrees

%% ------------------------------------------------------------ loop over the days
files = dir(fullfile(SUB_DIR, 'WRF_subset_*.mat'));
assert(~isempty(files), 'no WRF_subset files in %s', SUB_DIR);
n = numel(files);
Date = strings(n,1); sector = strings(n,1); ratio = NaN(n,1); sdRatio = NaN(n,1); upRatio = NaN(n,1);
wsDay = NaN(n,1); wdDay = NaN(n,1); cityMean = NaN(n,1); bgSD = NaN(n,1);
for k = 1:n
    S = load(fullfile(SUB_DIR, files(k).name), 'grid', 'hours', 'tracerHours', 'U', 'V', 'TRACER');
    Date(k) = regexp(files(k).name, '\d{8}', 'match', 'once');
    ib = BOX_I - S.grid.i_range(1) + 1;  jb = BOX_J - S.grid.j_range(1) + 1;     % box in subset indices
    assert(ib(1) >= 1 && ib(end) <= size(S.U,1) && jb(1) >= 1 && jb(end) <= size(S.U,2), '%s: analysis box outside the subset', Date(k));
    [ok1, hi] = ismember(HOURS, S.hours);  [ok2, ti] = ismember(HOURS, S.tracerHours);
    assert(all(ok1) && all(ok2), '%s: hours %s not all in the subset', Date(k), mat2str(HOURS));

    % ---- grid-relative wind at LEVEL (the subset stores earth-relative U, V)
    ca = double(S.grid.COSALPHA(ib,jb)); sa = double(S.grid.SINALPHA(ib,jb));
    Ue = double(squeeze(S.U(ib,jb,LEVEL,hi))); Ve = double(squeeze(S.V(ib,jb,LEVEL,hi)));
    Ug = Ue.*ca + Ve.*sa;  Vg = Ve.*ca - Ue.*sa;
    ws = sqrt(Ug.^2 + Vg.^2);  wd = mod(atan2d(Ug, Vg) + 180, 360);            % direction the wind blows from
    wsDay(k) = mean(ws(:), 'omitnan');
    wdPt = atan2d(mean(sind(wd), 3), mean(cosd(wd), 3));                       % per point, mean over the hours
    wdDay(k) = mod(atan2d(mean(sind(wdPt(:))), mean(cosd(wdPt(:)))), 360);    % then over the box
    sec = sectorNames{mod(round(wdDay(k)/45), 8) + 1}; sector(k) = sec;

    % ---- tracers (ppb) at LEVEL: regional background and city
    bg = double(squeeze(sum(S.TRACER(ib,jb,LEVEL,ti,1:5), 5)));                 % nx x ny x nHours
    ct = double(squeeze(sum(S.TRACER(ib,jb,LEVEL,ti,6:8), 5)));
    dn = L.(sec).dn;  up = L.(sec).up;
    nh = numel(HOURS) - 1; stn = NaN(1,nh); sdr = stn; upr = stn; cm = stn; bs = stn;
    for h = 1:nh                      % downwind line at HOURS(h+1), upwind line one hour earlier
        c = pick(ct, dn, h+1);  b = pick(bg, dn, h+1);  u = pick(bg, up, h);
        cm(h) = mean(c); bs(h) = std(b);
        stn(h) = mean(c)/std(b);      % the diagnostic ratio (signal to noise) of this hour
        sdr(h) = std(b)/std(c);       % secondary: background SD over city SD on the downwind line
        upr(h) = std(u)/std(c);       % secondary: upwind background SD over city SD
    end
    ratio(k) = mean(stn, 'omitnan'); sdRatio(k) = mean(sdr, 'omitnan'); upRatio(k) = mean(upr, 'omitnan');
    cityMean(k) = mean(cm); bgSD(k) = mean(bs);
    fprintf('%2d %s  wind %4.1f m/s %5.1f deg (%2s)  city %5.1f ppb  bg SD %5.1f ppb  ratio %5.2f\n', ...
        k, Date(k), wsDay(k), wdDay(k), sec, cityMean(k), bgSD(k), ratio(k));
end

%% ------------------------------------------------------------ Downwind Edge emission rates (optional)
flux = NaN(n,1);
if isfile(SOL_FILE)
    SA = load(SOL_FILE, 'SolA'); SolA = SA.SolA;
    for k = 1:n
        j = find(strcmp({SolA.Date}, char(Date(k))), 1);
        if ~isempty(j), flux(k) = mean(SolA(j).Flux); end
    end
end

%% ------------------------------------------------------------ output
T = table(Date, sector, ratio, sdRatio, upRatio, cityMean, bgSD, wsDay, wdDay, flux, 'VariableNames', ...
    {'Date', 'Sector', 'Diagnostic_ratio', 'SD_ratio_background_city', 'SD_ratio_upwind_city', ...
     'City_enhancement_mean_ppb', 'Background_SD_ppb', 'Mean_wind_speed_ms', 'Mean_wind_direction_deg', 'Flux_DownwindEdge_mols'});
T.Date = arrayfun(@(d) sprintf('%s-%s-%s', d{1}(1:4), d{1}(5:6), d{1}(7:8)), cellstr(Date), 'UniformOutput', false);
writetable(T, fullfile(OUT_DIR, 'diagnostic_ratios.csv'));
Diag_All2 = [wsDay, wdDay, ratio, flux, str2double(Date)];      % layout of Diag4.mat
save(fullfile(OUT_DIR, 'Diag_ratio.mat'), 'Diag_All2', 'T');
fprintf('\n%d days; ratio > %g on %d days\n', n, THRESH, nnz(ratio > THRESH));

if isfile(DIAG_OLD)
    O = load(DIAG_OLD); O = O.Diag_All2;
    fprintf('\ncomparison with %s (old / new):\n', DIAG_OLD);
    for k = 1:n
        j = find(O(:,5) == str2double(Date(k)), 1); if isempty(j), continue; end
        flag = ''; if (O(j,3) > THRESH) ~= (ratio(k) > THRESH), flag = '  <-- changes side of the threshold'; end
        fprintf('  %s  ratio %5.2f / %5.2f   wind %4.1f / %4.1f m/s  %5.1f / %5.1f deg%s\n', Date(k), O(j,3), ratio(k), O(j,1), wsDay(k), O(j,2), wdDay(k), flag);
    end
end
fprintf('wrote diagnostic_ratios.csv and Diag_ratio.mat in %s\n', OUT_DIR);

%% ------------------------------------------------------------ helpers
function s = tlines(xu, yu, xd, yd)
% upwind and downwind lines as lists of [x y] box indices (scalars are expanded)
    nu = max(numel(xu), numel(yu)); nd = max(numel(xd), numel(yd));
    s.up = [expand(xu, nu), expand(yu, nu)];  s.dn = [expand(xd, nd), expand(yd, nd)];
end
function v = expand(v, n)
    v = v(:); if numel(v) == 1, v = repmat(v, n, 1); end
end
function v = pick(A, pts, h)
% values of the nx x ny x nHours array A at the box points pts (rows [x y]) and hour index h
    v = A(sub2ind(size(A), pts(:,1), pts(:,2), repmat(h, size(pts,1), 1)));
end

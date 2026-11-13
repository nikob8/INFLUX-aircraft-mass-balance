%% ========================================================================
%  Get_Flight_Wind_Directions.m
%  ------------------------------------------------------------------------
%  Determine the mean wind direction (and speed) used for each of the 37
%  INFLUX mass balance flight days and write it into a copy of
%  Indy_Mass_Balance_Table.xlsx as a new column on the date rows.
%
%  Two independent estimates are made for every day:
%    (1) MB   : MeanWDIR / MeanWSPD stored in Solution-MethodA.mat, i.e. the
%               wind the mass balance calculation actually used (if present)
%    (2) OBS  : vector mean of the aircraft winds over ALL downwind legs of
%               the day, legs found exactly as in Compare_WRF_Obs_v5.m
%               (DownwindTime anchors + LegHeight).  This mirrors how the
%               mass balance wind was defined: one wind per day from all
%               downwind transects.
%  The OBS speed is compared with column E of the spreadsheet (the speed
%  used in the paper) as a check that the legs and averaging reproduce the
%  mass balance wind.  The value written to the spreadsheet is MB when it
%  exists, otherwise OBS.
%
%  Also reported: leg-to-leg spread of wind direction within each day
%  (circular SD), which is useful for Reviewer 2 (Line 250, representativeness
%  of a single mean wind direction).
%
%  Outputs (in OUT_DIR):
%     Wind_by_day.csv     one row per flight day
%     Wind_by_leg.csv     one row per downwind leg
%     Indy_Mass_Balance_Table_withWD.xlsx   copy of the table + new columns
%
%  Helper functions are copied verbatim from Compare_WRF_Obs_v5.m.
%  Nikolay Balashov -- INFLUX aircraft paper, reviewer revision
%  ========================================================================

clc; clear; close all;

%% ---------------- settings (same paths as Compare_WRF_Obs_v5.m) ---------
ALL_DATES = {'20110301','20110429','20110601','20110630','20110712', ...
             '20120606','20121004','20121108','20121218','20130314', ...
             '20130319','20130814','20131107','20140224','20140725', ...
             '20140814','20141110','20141113','20141114','20141117', ...
             '20141119','20141120','20141121','20141201','20141203', ...
             '20150706','20150714','20150720','20150811','20150820', ...
             '20150825','20150902','20160229','20160306','20160307', ...
             '20160317','20160318'};

cfg.GUI_ROOT    = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0';
cfg.OBS_DIR     = fullfile(cfg.GUI_ROOT,'Observations','INFLUX');
cfg.VTABLE_FILE = fullfile(cfg.OBS_DIR,'VTable_INFLUX');
cfg.MB_ROOT     = '/storage/work/nvb5011/INFLUX_Paper_2/MassBalanceResults2';
cfg.MB_SOLUTION = 'Solution-MethodA.mat';
cfg.XLS_IN      = 'Indy_Mass_Balance_Table.xlsx';     % in the current folder
cfg.OUT_DIR     = fullfile(pwd,'wind_dir');
cfg.MISSING_VAL = -99;
cfg.MIN_LEG_KEEP = 240;      % s, drop climbs/descents (as in v5)
cfg.FLAG_DEG    = 15;        % flag days where MB and OBS directions differ by more

legOpt = struct('alttol',40,'tolmax',150,'minsec',300,'zwin',60, ...
                'search',900,'gap',20,'trim',true,'fence',true);   % = v5 settings

if ~isfolder(cfg.OUT_DIR), mkdir(cfg.OUT_DIR); end

%% ---------------- loop over flights -------------------------------------
nD = numel(ALL_DATES);
DAY = table('Size',[nD 13], ...
    'VariableTypes',{'string','double','double','double','double','double', ...
                     'double','double','double','double','double','double','string'}, ...
    'VariableNames',{'date','nLegs','wd_MB','ws_MB','wd_obs','ws_obs_vec', ...
                     'ws_obs_scalar','wd_legSD','wd_legRange','dWD_MB_obs', ...
                     'wd_final','ws_xls','note'});
LEG = table();

for f = 1:nD
    d = ALL_DATES{f};
    DAY.date(f) = string(d);  DAY{f,2:12} = nan(1,11);  DAY.note(f) = "";
    fprintf('\n=== %s/%s/%s (%d of %d) ===\n', d(5:6), d(7:8), d(1:4), f, nD);
    try
        % ---- mass balance solution: stored wind, if any
        mbFile = fullfile(cfg.MB_ROOT, d, cfg.MB_SOLUTION);
        MB = struct();
        if isfile(mbFile), MB = load(mbFile); end
        if isfield(MB,'MeanWDIR'), DAY.wd_MB(f) = MB.MeanWDIR(1); end
        if isfield(MB,'MeanWSPD'), DAY.ws_MB(f) = MB.MeanWSPD(1); end

        % ---- aircraft winds on the downwind legs
        obs = load_influx_obs(cfg.OBS_DIR, cfg.VTABLE_FILE, d, cfg.MISSING_VAL);
        if max(obs.alt) < 20, obs.alt = obs.alt*1000; end
        [TR, ~] = build_transects_legheight(obs, mbFile, legOpt);
        TR = TR(arrayfun(@(t) t.s1 - t.s0, TR) >= cfg.MIN_LEG_KEEP);
        if isempty(TR), error('no usable downwind legs'); end

        allIdx = [];  wdLeg = nan(numel(TR),1);
        for k = 1:numel(TR)
            ii = TR(k).idx;  allIdx = [allIdx; ii]; %#ok<AGROW>
            [wdLeg(k), wsK, nK] = mean_wind_from_wdws(obs.wdir(ii), obs.wspd(ii));
            LEG = [LEG; table(string(d), string(TR(k).label), string(hhmmss(TR(k).s0)), ...
                   string(hhmmss(TR(k).s1)), mean(obs.alt(ii),'omitnan'), wdLeg(k), wsK, nK, ...
                   'VariableNames',{'date','leg','t0','t1','altASL','wd_obs','ws_obs','n_wind'})]; %#ok<AGROW>
        end
        % one wind per day from ALL downwind legs pooled (as in the mass balance)
        [DAY.wd_obs(f), DAY.ws_obs_vec(f)] = mean_wind_from_wdws(obs.wdir(allIdx), obs.wspd(allIdx));
        DAY.ws_obs_scalar(f) = mean(obs.wspd(allIdx),'omitnan');
        DAY.nLegs(f)       = numel(TR);
        DAY.wd_legSD(f)    = circ_std_deg(wdLeg);
        dl = wrap180(wdLeg - DAY.wd_obs(f));
        DAY.wd_legRange(f) = max(dl) - min(dl);
        if ~isnan(DAY.wd_MB(f))
            DAY.dWD_MB_obs(f) = wrap180(DAY.wd_obs(f) - DAY.wd_MB(f));
        end
    catch ME
        DAY.note(f) = "FAILED: " + string(ME.message);
        fprintf(2, '*** %s failed: %s\n', d, ME.message);
    end
    % value to use
    if ~isnan(DAY.wd_MB(f)), DAY.wd_final(f) = DAY.wd_MB(f);
    else,                    DAY.wd_final(f) = DAY.wd_obs(f);
        if DAY.note(f) == "", DAY.note(f) = "no MeanWDIR in .mat; aircraft value used"; end
    end
    if abs(DAY.dWD_MB_obs(f)) > cfg.FLAG_DEG
        DAY.note(f) = strtrim(DAY.note(f) + " CHECK: MB and aircraft differ by " + ...
                      sprintf('%.0f deg', DAY.dWD_MB_obs(f)));
    end
end

%% ---------------- compare with the speed in the spreadsheet -------------
C = readcell(cfg.XLS_IN, 'Sheet','FlightMasterSheet');
dateRow = nan(nD,1);
for r = 2:size(C,1)
    lab = C{r,1};
    if isa(lab,'datetime'), s = datestr(lab,'yyyymmdd');
    elseif ischar(lab) || isstring(lab), s = strrep(strtrim(char(lab)),'/','');
    else, continue; end
    j = find(strcmp(ALL_DATES, s), 1);
    if ~isempty(j), dateRow(j) = r; end
end
for f = 1:nD
    if ~isnan(dateRow(f)) && isnumeric(C{dateRow(f),5}), DAY.ws_xls(f) = C{dateRow(f),5}; end
end

fprintf('\n%-9s %4s %7s %7s %7s | %6s %6s %6s %6s | %6s  %s\n', 'date','legs', ...
    'WD_MB','WD_obs','dWD','WS_xls','WS_MB','WSvec','WSscl','legSD','note');
for f = 1:nD
    fprintf('%-9s %4d %7.0f %7.0f %+7.0f | %6.2f %6.2f %6.2f %6.2f | %6.0f  %s\n', DAY.date(f), ...
        DAY.nLegs(f), DAY.wd_MB(f), DAY.wd_obs(f), DAY.dWD_MB_obs(f), DAY.ws_xls(f), ...
        DAY.ws_MB(f), DAY.ws_obs_vec(f), DAY.ws_obs_scalar(f), DAY.wd_legSD(f), DAY.note(f));
end
ok = ~isnan(DAY.ws_xls);
fprintf('\nSpeed check vs spreadsheet column E (which aircraft average reproduces it?):\n');
fprintf('   |WS_MB     - WS_xls| median %.2f m/s\n', median(abs(DAY.ws_MB(ok)        - DAY.ws_xls(ok)),'omitnan'));
fprintf('   |WS_vector - WS_xls| median %.2f m/s\n', median(abs(DAY.ws_obs_vec(ok)   - DAY.ws_xls(ok)),'omitnan'));
fprintf('   |WS_scalar - WS_xls| median %.2f m/s\n', median(abs(DAY.ws_obs_scalar(ok)- DAY.ws_xls(ok)),'omitnan'));
fprintf('Within-day leg-to-leg wind direction SD: median %.0f deg, max %.0f deg\n', ...
    median(DAY.wd_legSD,'omitnan'), max(DAY.wd_legSD));

writetable(DAY, fullfile(cfg.OUT_DIR,'Wind_by_day.csv'));
writetable(LEG, fullfile(cfg.OUT_DIR,'Wind_by_leg.csv'));

%% ---------------- write into a COPY of the spreadsheet -----------------
xlsOut = fullfile(cfg.OUT_DIR, 'Indy_Mass_Balance_Table_withWD.xlsx');
copyfile(cfg.XLS_IN, xlsOut);
% first empty columns after column O (15): P = direction, Q = leg-to-leg SD
writecell({sprintf('Mean Flight\nWind Direction (deg)'), sprintf('Leg-to-leg Wind\nDirection SD (deg)')}, ...
          xlsOut, 'Sheet','FlightMasterSheet', 'Range','P1');
for f = 1:nD
    if isnan(dateRow(f)) || isnan(DAY.wd_final(f)), continue; end
    writematrix([round(DAY.wd_final(f)) round(DAY.wd_legSD(f))], xlsOut, ...
        'Sheet','FlightMasterSheet', 'Range', sprintf('P%d', dateRow(f)));
end
fprintf('\nWrote %s\n(columns P and Q on the date rows). Check flagged days before using it.\n', xlsOut);
fprintf('Done.\n');


%% ======================================================================
%  Helper functions (copied from Compare_WRF_Obs_v5.m)
%  ======================================================================


function obs = load_influx_obs(obsDir, vtableFile, caseDate, missingVal)
% INFLUX_YYYYMMDD is plain ASCII, no header, columns per VTable_INFLUX:
%   TIME 1  LAT 2  LON 3  ALT 4  CO2 5  CH4 6  H2O 7
%   CO 8    Ps 9   T 10   WSPD 11  WDIR 12
% TIME is seconds of day (UTC).  ALT is metres above sea level.

    fname = fullfile(obsDir, ['INFLUX_' caseDate]);
    if ~isfile(fname), error('Observation file not found: %s', fname); end

    col  = read_vtable(vtableFile);
    nCol = max(cell2mat(struct2cell(col)));
    A = [];
    try, A = readmatrix(fname,'FileType','text'); catch, A = []; end
    if isempty(A) || size(A,2) < nCol, A = load(fname); end
    if size(A,2) < nCol
        error('%s has %d columns, VTable expects at least %d.', fname, size(A,2), nCol);
    end

    obs.secs = A(:,col.TIME);  obs.lat = A(:,col.LAT);  obs.lon = A(:,col.LON);
    obs.alt  = A(:,col.ALT);   obs.ch4 = A(:,col.CH4);
    gv = @(f) iff(isfield(col,f), @() A(:,col.(f)), @() nan(size(A,1),1));
    obs.co2 = gv('CO2');  obs.wspd = gv('WSPD');  obs.wdir = gv('WDIR');

    obs.secs(obs.secs < 0)           = NaN;
    obs.lat (obs.lat  < -90)         = NaN;
    obs.lon (obs.lon  < -360)        = NaN;
    obs.alt (obs.alt  <= missingVal) = NaN;
    obs.ch4 (obs.ch4  <= missingVal) = NaN;
    obs.co2 (obs.co2  <= missingVal) = NaN;
    obs.wspd(obs.wspd <= missingVal) = NaN;
    obs.wdir(obs.wdir <= missingVal) = NaN;
    obs.lon(obs.lon > 180) = obs.lon(obs.lon > 180) - 360;

    ok = ~isnan(obs.secs) & ~isnan(obs.lat) & ~isnan(obs.lon) & ...
         ~isnan(obs.alt)  & ~isnan(obs.ch4);
    f = fieldnames(obs);
    for i = 1:numel(f), obs.(f{i}) = obs.(f{i})(ok); end
    [~,srt] = sort(obs.secs);
    for i = 1:numel(f), obs.(f{i}) = obs.(f{i})(srt); end
    fprintf('Read %s : %d valid of %d records.\n', fname, sum(ok), numel(ok));
end

function v = iff(c, fT, fF)
    if c, v = fT(); else, v = fF(); end
end

function col = read_vtable(vtableFile)
    def = struct('TIME',1,'LAT',2,'LON',3,'ALT',4,'CO2',5,'CH4',6, ...
                 'H2O',7,'CO',8,'Ps',9,'T',10,'WSPD',11,'WDIR',12);
    if ~isfile(vtableFile)
        warning('VTable not found -- using the default INFLUX column map.');
        col = def; return;
    end
    col = struct();
    fid = fopen(vtableFile);  fgetl(fid);
    while ~feof(fid)
        ln = strtrim(fgetl(fid));
        if isempty(ln), continue; end
        p = strsplit(ln);
        if numel(p) < 2, continue; end
        v = str2double(p{2});
        if ~isnan(v) && v >= 1, col.(matlab.lang.makeValidName(p{1})) = v; end
    end
    fclose(fid);
    need = {'TIME','LAT','LON','ALT','CH4'};
    for i = 1:numel(need)
        if ~isfield(col,need{i}), col.(need{i}) = def.(need{i}); end
    end
end

function [TR, MB] = build_transects_legheight(obs, mbFile, opt)
% Anchor on DownwindTime; select by LegHeight (AGL) + dZ, where dZ is the
% terrain offset re-derived from the data.
    if ~isfile(mbFile)
        error('Mass balance solution not found: %s', mbFile);
    end
    MB = load(mbFile);
    if ~isfield(MB,'DownwindTime'), error('%s contains no DownwindTime.', mbFile); end

    anchors = MB.DownwindTime(:);            % force columns -- shapes in the
    nA      = numel(anchors);                % .mat are not consistent
    hasLH   = isfield(MB,'LegHeight') && numel(MB.LegHeight) >= nA;
    hasFl   = isfield(MB,'Flux')      && numel(MB.Flux)      >= nA;
    LH = []; FL = [];
    if hasLH, LH = MB.LegHeight(:); end
    if hasFl, FL = MB.Flux(:);      end

    zAtAnchor = nan(nA,1);
    for k = 1:nA
        w = abs(obs.secs - anchors(k)) <= 60;
        if any(w), zAtAnchor(k) = median(obs.alt(w),'omitnan'); end
    end
    if hasLH, dZ = median(zAtAnchor - LH(1:nA),'omitnan'); else, dZ = 0; end
    if ~isscalar(dZ), error('dZ is not scalar -- shape problem in %s', mbFile); end

    fprintf('Anchors: dZ (ALT_ASL - LegHeight) = %.0f m', dZ);
    if hasLH
        resid = zAtAnchor - LH(1:nA) - dZ;
        fprintf('   residuals %s m', strjoin(compose('%+.0f', resid(:)'), ' '));
    end
    fprintf('\n');
    for k = 1:nA
        fprintf('   T%d  %s', k, hhmmss(anchors(k)));
        if hasLH, fprintf('   LegHeight %6.1f m AGL -> %6.1f m ASL', LH(k), LH(k)+dZ); end
        if hasFl, fprintf('   Flux %5.1f mol/s', FL(k)); end
        fprintf('\n');
    end

    % Fences: half way to the neighbouring anchor in time.  Without these,
    % two legs flown at similar altitudes close together in time can both
    % select the same stretch of track.
    fLo = -inf(nA,1);  fHi = inf(nA,1);
    if opt.fence && nA > 1
        [as, ord] = sort(anchors);
        mid = (as(1:end-1) + as(2:end))/2;
        for m = 2:nA,   fLo(ord(m))   = mid(m-1); end
        for m = 1:nA-1, fHi(ord(m))   = mid(m);   end
    end

    TR = struct('label',{},'s0',{},'s1',{},'idx',{},'anchor',{},'dZ',{});
    for k = 1:nA
        % Target altitude: what the aircraft was ACTUALLY flying at the
        % anchor.  LegHeight + dZ is only a fallback -- a single day-wide
        % dZ misses legs whose residual about it is large.
        target = zAtAnchor(k);
        if isnan(target) && hasLH, target = LH(k) + dZ; end
        if isnan(target), warning('T%d: no target altitude.', k); continue; end

        tLo   = max(anchors(k) - opt.search, fLo(k));
        tHi   = min(anchors(k) + opt.search, fHi(k));
        inWin = obs.secs >= tLo & obs.secs <= tHi;

        % Widen the altitude band until the leg reaches a sensible length.
        tol = opt.alttol;  idx = [];  used = tol;
        while true
            cand = find(inWin & abs(obs.alt - target) <= tol);
            if ~isempty(cand)
                brk = [0; find(diff(obs.secs(cand)) > opt.gap); numel(cand)];
                bestG = 1;  bestD = Inf;
                for g = 1:numel(brk)-1
                    seg = cand(brk(g)+1 : brk(g+1));
                    d = max(0, max(anchors(k)-obs.secs(seg(end)), ...
                                   obs.secs(seg(1))-anchors(k)));
                    if d < bestD, bestD = d; bestG = g; end
                end
                grp = cand(brk(bestG)+1 : brk(bestG+1));
                if opt.trim, grp = trim_turnaround(obs, grp, anchors(k)); end
                if numel(grp) > numel(idx), idx = grp; used = tol; end
                if ~isempty(grp) && ...
                   obs.secs(grp(end)) - obs.secs(grp(1)) >= opt.minsec
                    break
                end
            end
            tol = tol*1.5;
            if tol > opt.tolmax, break; end
        end

        if numel(idx) < 10
            warning('T%d: only %d points near %.0f m ASL at %s.', ...
                k, numel(idx), target, hhmmss(anchors(k)));
            continue;
        end
        if used > opt.alttol
            fprintf('   T%d: altitude band widened to +/-%.0f m.\n', k, used);
        end

        TR(end+1) = struct('label',sprintf('T%d',k),'s0',obs.secs(idx(1)), ...
            's1',obs.secs(idx(end)),'idx',idx(:),'anchor',anchors(k),'dZ',dZ); %#ok<AGROW>
    end

    % flag any leg that still shares points with another
    for a = 1:numel(TR)
        for b = a+1:numel(TR)
            nsh = numel(intersect(TR(a).idx, TR(b).idx));
            if nsh > 0
                warning('%s and %s share %d points -- check the windows.', ...
                    TR(a).label, TR(b).label, nsh);
            end
        end
    end
end

function idx = trim_turnaround(obs, idx, tAnchor)
% If the aircraft holds altitude through a turn, the selected stretch can
% double back on itself.  Keep only the part travelling in one direction.
%
% The direction is the leg's OWN principal axis, found by PCA on the track
% in km -- NOT longitude.  Many INFLUX legs run north-south or diagonally,
% and on those the longitude derivative is noise, which is what shredded
% those legs in v4 (07/20/2015, 11/07/2013, 03/19/2013, 11/08/2012).
    if numel(idx) < 20, return; end

    lat = obs.lat(idx);  lon = obs.lon(idx);
    la0 = mean(lat, 'omitnan');
    x = (lon - mean(lon,'omitnan')) * 111.32 * cosd(la0);   % km east
    y = (lat - la0)                 * 110.57;               % km north

    ok = ~isnan(x) & ~isnan(y);
    if nnz(ok) < 10, return; end
    [V, D] = eig(cov([x(ok) y(ok)]));
    [~, m] = max(diag(D));
    p = [x y] * V(:,m);            % distance along the principal axis, km

    d = movmean([0; diff(p)], 31, 'omitnan');
    sgn = sign(sum(d, 'omitnan'));
    if sgn == 0, return; end
    keep = (sign(d) == sgn) | (abs(d) < 1e-4);

    dd = diff([false; keep(:); false]);
    rs = find(dd == 1);  re = find(dd == -1) - 1;
    if isempty(rs), return; end
    [~, ja] = min(abs(obs.secs(idx) - tAnchor));
    hit = find(rs <= ja & re >= ja, 1);
    if isempty(hit), [~, hit] = max(re - rs); end
    idx = idx(rs(hit):re(hit));
end

function [wd, ws, n] = mean_wind_from_wdws(wdir, wspd)
% Vector mean of aircraft winds given from-direction (deg) and speed.
    ok = ~isnan(wdir) & ~isnan(wspd);
    n  = nnz(ok);
    if n == 0, wd = NaN; ws = NaN; return; end
    u = -wspd(ok).*sind(wdir(ok));      % toward-components
    v = -wspd(ok).*cosd(wdir(ok));
    [wd, ws] = mean_wind_from_uv(u, v);
end

function [wd, ws] = mean_wind_from_uv(u, v)
% Same convention as the WRFGUI: atan2d(U,V)+180 gives the from-direction.
    um = mean(u(:),'omitnan');  vm = mean(v(:),'omitnan');
    wd = mod(atan2d(um, vm) + 180, 360);
    ws = mean(sqrt(u(:).^2 + v(:).^2),'omitnan');
end

function a = wrap180(a)
    a = mod(a + 180, 360) - 180;
end

function m = circ_mean_deg(d)
    d = d(~isnan(d));
    if isempty(d), m = NaN; return; end
    m = mod(atan2d(mean(sind(d)), mean(cosd(d))), 360);
end

function s = circ_std_deg(d)
    d = d(~isnan(d));
    if numel(d) < 2, s = NaN; return; end
    R = hypot(mean(sind(d)), mean(cosd(d)));
    s = rad2deg(sqrt(-2*log(max(R,1e-9))));
end

function s = hhmmss(sec)
    if isnan(sec), s = '  --  '; return; end
    h = floor(sec/3600);  m = floor((sec-h*3600)/60);  ss = round(sec-h*3600-m*60);
    if ss == 60, ss = 0; m = m + 1; end
    if m  == 60, m  = 0; h = h + 1; end
    s = sprintf('%02d:%02d:%02d', h, m, ss);
end
%% Trajectory_check_R2_7.m
%  Reviewer 2 (R2.7, R2.8, R2.11): does the straight-line, mean-wind projection
%  of the county boundary onto the downwind transects misplace the edges of the
%  city segment, and by how much would that change the fluxes?
%
%  For every downwind transect of every flight day a back trajectory is started
%  from each 1 Hz point of the transect, at the time the aircraft was there, and
%  followed for up to TMAX hours with four wind models:
%     straight : the constant mean wind of the mass balance (Method A).  This
%                reproduces the published edges and fluxes and is the check.
%     legwind  : the mean wind of that transect only, constant in time.
%     obs      : the wind measured by the aircraft on every transect (upwind and
%                downwind), interpolated in time, uniform in space: the parcels
%                follow the observed change of the wind during the flight.
%     wrf      : the WRF-Chem 3 km hourly winds of the run used for the
%                diagnostic ratio, averaged over the observed ABL depth and
%                interpolated in space and time.
%  The points whose trajectory passes through Marion County form the city
%  segment of the transect; its two ends are the edges.  The Edge flux is then
%  recomputed exactly as in MassBalanceMethodA.m with the new edges, and the
%  Upwind Edge flux (Method D) with the new edges and the trajectory-based
%  pairing of downwind and upwind points (R2.11).  The script also compares the
%  WRF wind with the aircraft wind at every transect, measures the spatial and
%  temporal homogeneity of the modelled flow (R2.8), and computes a
%  profile-shape (trough-fit) background for comparison with the geometric one
%  (R2.7, second part).
%
%  Inputs (the same files the Method scripts use)
%     FlightList.mat, LegAndProfileList/<Date>.mat, INFLUX_<Date>
%     <Date>/Solution-MethodA.mat and Solution-MethodD.mat   (your MassB folder)
%     wrfout_d02_* files, found with the same lookup as Figure7_v4.m
%  Outputs (in Traj_results/)
%     <Date>.mat                    per-transect results of one day
%     Trajectory_check_legs.csv     one row per transect, all numbers
%     Trajectory_check_log.txt      what happened on every day
%     FigS6_trajectory_check.png    six-panel figure for the SI
%  and a printed summary with the numbers for the responses.
%
%  Run from any folder; it takes a while (it reads about 8 wrfout files per day).
%  If it stops, run it again: days already in Traj_results/ are skipped.

clear; close all; clc;

%% ------------------------------------------------------------ settings
MB_ROOT   = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';
LEG_DIR   = fullfile(MB_ROOT, 'LegAndProfileList');
OBS_DIR   = '/storage/work/zrb5027/WRFGUI_V3/Observations/INFLUX';            % INFLUX_<Date>
SOL_ROOT  = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB';
WRF_ROOT  = '/storage/group/zrb5027/default/INFLUX/nvb5011/';
DIAG_FILE = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/Diag4.mat';       % optional: tags the ratio > 2 days
OUT_DIR   = fullfile(pwd, 'Traj_results');

P.DOMAIN   = 2;   P.KMAX = 25;    % WRF domain (3 km) and number of model levels read
P.DT       = 30;                  % trajectory time step (s)
P.TMAX     = 8*3600;              % longest back trajectory (s); 8 h so that the slow-wind days reach the upwind transect
P.BOX_KM   = 160;                 % half-width of the WRF sub-domain kept in memory (km)
P.HOMOG_KM = 50;                  % radius for the spatial homogeneity of the modelled wind (R2.8)
P.EXAMPLE_DATE = '20110301';  P.EXAMPLE_HOUR = 17.33;   % transect drawn in the figure (the Figure 2 transect)
START_AT  = 1;  SKIP_DONE = true;  % SKIP_DONE = true reuses Traj_results/<Date>.mat from an earlier run (only the summary and the figure are remade)
FN = 'Arial';

%BottomRight BottomLeft TopLeft TopRight  (as in the Method scripts)
P.CC = [39.63760423971629, -85.95148383144084; 39.63337385377353, -86.32570560518701; ...
        39.92446141560437, -86.32734193809195; 39.9271330073512, -85.95582253907553];
P.LATC = mean(P.CC(:,1));  P.LONC = mean(P.CC(:,2));
P.WRF_ROOT = WRF_ROOT;  P.MODELS = {'straight','legwind','obs','wrf'};

addpath(MB_ROOT); addpath(LEG_DIR); addpath(OBS_DIR);
if ~isfolder(OUT_DIR), mkdir(OUT_DIR); end
logFile = fullfile(OUT_DIR, 'Trajectory_check_log.txt');
fid = fopen(logFile, 'a'); fprintf(fid, '\n==== %s ====\n', datestr(now)); fclose(fid);

FL = load('FlightList'); FlightList = FL.FlightList;
fprintf('%d flight days in FlightList\n', numel(FlightList));

%% ------------------------------------------------------------ loop over the days
for k = START_AT:numel(FlightList)
    Date = FlightList(k).name; Date = Date(1:8);
    outFile = fullfile(OUT_DIR, [Date '.mat']);
    if SKIP_DONE && isfile(outFile), fprintf('%2d %s: done earlier, skipped\n', k, Date); continue; end
    t0 = tic;
    try
        D = load_flight(Date, LEG_DIR, SOL_ROOT);
        W = load_wrf_day(D, P);
        R = run_day(D, W, P);
        save(outFile, '-struct', 'R');
        msg = sprintf('%d transects, %s, %.0f s', numel(R.leg), R.wrfNote, toc(t0));
    catch err
        msg = ['ERROR: ' err.message];
        if ~isempty(err.stack), msg = sprintf('%s (line %d of %s)', msg, err.stack(1).line, err.stack(1).name); end
    end
    txt = sprintf('%2d %s: %s', k, Date, msg); disp(txt);
    fid = fopen(logFile, 'a'); fprintf(fid, '%s\n', txt); fclose(fid);
end

%% ------------------------------------------------------------ collect, summarise, plot
collect_and_report(OUT_DIR, DIAG_FILE, P, FN);


%% ============================================================ one flight day
function R = run_day(D, W, P)
nLeg = size(D.LegList, 1);
kx = @(lon) (lon - P.LONC)*111.32*cosd(P.LATC);  ky = @(lat) (lat - P.LATC)*110.57;
cty.x = kx(P.CC(:,2)); cty.y = ky(P.CC(:,1));
uhat = [D.MeanU D.MeanV]/hypot(D.MeanU, D.MeanV);

% upwind transect: axis, s-coordinate of its points, wind comparison
up = []; sUpPts = [];
R.up = struct('exists', D.hasUp);
if D.hasUp
    xu = kx(D.Lon(D.UpLeg)); yu = ky(D.Lat(D.UpLeg)); oku = ~isnan(xu) & ~isnan(yu);
    Au = leg_axis(xu(oku), yu(oku), uhat);
    up = struct('c', Au.c, 'e', Au.e, 'nv', Au.nv);
    sUpPts = (xu - Au.c(1))*Au.e(1) + (yu - Au.c(2))*Au.e(2);
    R.up.s = sUpPts;  R.up.x = xu;  R.up.y = yu;  R.up.ch = D.CH4(D.UpLeg);
    R.up.t = nmean(D.Time(D.UpLeg));
    [R.up.dirAir, R.up.spdAir] = dir_spd(nmean(D.U(D.UpLeg)), nmean(D.V(D.UpLeg)));
    R.up.dirWrf = NaN; R.up.spdWrf = NaN;
    if ~isempty(W)
        [uw, vw] = wrf_wind(nmedian(xu), nmedian(yu), R.up.t, W);
        [R.up.dirWrf, R.up.spdWrf] = dir_spd(uw, vw);
    end
end

% wind models
nM = numel(P.MODELS);
wf = cell(1, nM);
wf{1} = @(x, y, t) const_wind(x, D.MeanU, D.MeanV);
wf{3} = @(x, y, t) obs_wind(t, D.obsT, D.obsU, D.obsV);
if ~isempty(W), wf{4} = @(x, y, t) wrf_wind(x, y, t, W); end

% example transect of this day?
exLeg = 0;
if strcmp(D.Date, P.EXAMPLE_DATE)
    [~, exLeg] = min(abs(D.LegTime/3600 - P.EXAMPLE_HOUR));
end

R.Date = D.Date; R.MeanWSPD = D.MeanWSPD; R.MeanWDIR = D.MeanWDIR; R.hasD = D.hasD;
R.wrfNote = 'no WRF'; if ~isempty(W), R.wrfNote = W.note; end
R.leg = struct([]);
for L = 1:nLeg
    idx = D.LegList(L,1):D.LegList(L,2);
    leg.lat = D.Lat(idx); leg.lon = D.Lon(idx); leg.ch = D.CH4(idx); leg.t = D.Time(idx);
    n = numel(idx);
    x0 = kx(leg.lon); y0 = ky(leg.lat);
    okp = ~isnan(x0) & ~isnan(y0) & ~isnan(leg.t);
    A = leg_axis(x0(okp), y0(okp), uhat);
    s = NaN(n,1); s(okp) = A.s(x0(okp), y0(okp));

    % the transect's own wind (model 2) and the model/aircraft comparison
    uL = nmean(D.U(idx)); vL = nmean(D.V(idx));
    wf{2} = @(x, y, t) const_wind(x, uL, vL);

    % analytic straight-line edges: envelope of the four corners projected along the mean wind
    sC = NaN(4,1); M = [A.e' -uhat'];
    for c = 1:4
        if abs(det(M)) > 1e-6, ab = M \ ([cty.x(c) cty.y(c)]' - A.pO'); sC(c) = ab(1); end
    end
    Q = struct('Date', D.Date, 'Leg', L, 'hourUTC', D.LegTime(L)/3600, 'nPts', n, ...
        'height_m', D.LegHeight(L), 'ABL_m', D.ABL(L), 'orient', A.orient, ...
        'sMin', nmin(s), 'sMax', nmax(s), 'sEdgeAnalytic', [nmin(sC) nmax(sC)], ...
        'FluxA_file', D.FluxA(L), 'DistA_file', D.DistA(L), 'cosTheta', D.CosTheta(L), ...
        'dirAir', NaN, 'spdAir', NaN, 'dirDay', D.MeanWDIR, 'spdDay', D.MeanWSPD);
    [Q.dirAir, Q.spdAir] = dir_spd(uL, vL);
    Q.FluxD_file = NaN; if D.hasD, Q.FluxD_file = D.FluxD(L); end
    Q.example = [];

    % trajectories with each wind model
    Q.k1 = NaN(1,nM); Q.k2 = NaN(1,nM); Q.sEdge = NaN(nM,2); Q.contig = NaN(1,nM);
    Q.Flux = NaN(1,nM); Q.FluxD = NaN(1,nM); Q.lostFrac = NaN(1,nM);
    Q.transitNear_min = NaN(1,nM); Q.transitFar_min = NaN(1,nM);
    Q.upShiftMed_km = NaN(1,nM); Q.upShiftAbsMed_km = NaN(1,nM); Q.upShiftEdges_km = NaN(nM,2);
    Q.bg = NaN(nM,2); Q.width_km = NaN(1,nM); Q.pairFrac = NaN(1,nM);
    sUpAll = NaN(n, nM); ex = struct([]);
    for m = 1:nM
        if isempty(wf{m}), continue; end
        keepEvery = 0; if L == exLeg, keepEvery = 15; end
        T = back_traj(x0(okp), y0(okp), leg.t(okp), wf{m}, P.DT, P.TMAX, cty, up, keepEvery);
        ent = false(n,1); ent(okp) = T.entered;
        tin = NaN(n,1); tin(okp) = T.tIn;  tout = NaN(n,1); tout(okp) = T.tOut;
        Q.lostFrac(m) = mean(T.lost);
        kk = find(ent);
        if isempty(kk), continue; end                      % trajectories never touch the county
        k1 = kk(1); k2 = kk(end);
        Q.k1(m) = k1; Q.k2(m) = k2; Q.contig(m) = numel(kk)/(k2 - k1 + 1);
        Q.sEdge(m,:) = sort([s(k1) s(k2)]);
        Q.transitNear_min(m) = nmedian((leg.t(ent) - tin(ent))/60);
        Q.transitFar_min(m)  = nmedian((leg.t(ent) - tout(ent))/60);
        [Q.Flux(m), info] = edge_flux(leg, k1, k2, D.Density(L), D.MeanWSPD, D.CosTheta(L), D.ABL(L));
        Q.bg(m,:) = [info.bg1 info.bg2]; Q.width_km(m) = info.Dist;
        if D.hasUp
            su = NaN(n,1); okc = false(n,1);
            su(okp) = T.sUp; okc(okp) = T.okUp;
            okc = okc & su >= nmin(sUpPts) & su <= nmax(sUpPts);   % crossing within the upwind transect
            su(~okc) = NaN; sUpAll(:,m) = su;
            Q.pairFrac(m) = mean(okc(k1:k2));
            if D.hasD
                pair = NaN(n,1);
                [ss, io] = sort(sUpPts); good = ~isnan(ss); ss = ss(good); io = io(good);
                [ss, iu] = unique(ss, 'first'); io = io(iu);
                if numel(ss) >= 2
                    pair(okc) = interp1(ss, io, su(okc), 'nearest');
                end
                Q.FluxD(m) = upwind_edge_flux(leg, D.CH4(D.UpLeg), pair, k1, k2, D.CF(L), ...
                    D.Density(L), D.MeanWSPD, D.CosTheta(L), D.ABL(L));
            end
        end
        if keepEvery > 0
            ex(m).X = T.X; ex(m).Y = T.Y; ex(m).ks = T.ks;   %#ok<AGROW>
        end
    end
    % shift of the landing point on the upwind transect relative to the straight line (R2.11)
    if D.hasUp && ~isnan(Q.k1(1))
        seg = Q.k1(1):Q.k2(1);
        for m = 2:nM
            d = sUpAll(seg,m) - sUpAll(seg,1);
            Q.upShiftMed_km(m) = nmedian(d); Q.upShiftAbsMed_km(m) = nmedian(abs(d));
            Q.upShiftEdges_km(m,:) = [sUpAll(Q.k1(1),m) - sUpAll(Q.k1(1),1), sUpAll(Q.k2(1),m) - sUpAll(Q.k2(1),1)];
        end
    end
    % profile-shape (trough-fit) background with the published edges (R2.7, second part),
    % and how much of the plume signal the geometric edges leave outside the segment
    Q.FluxShapeSeg = NaN; Q.FluxShapeFull = NaN; Q.shapeBg = NaN(1,2);
    Q.FluxShape25Seg = NaN; Q.shapeBg25 = NaN(1,2);
    Q.edgeEnh_ppb = NaN(1,2); Q.outsideFrac = NaN;
    if ~isnan(Q.k1(1))
        [Q.FluxShapeSeg, Q.FluxShapeFull, sh] = shape_flux(leg, Q.k1(1), Q.k2(1), D.Density(L), D.MeanWSPD, D.CosTheta(L), D.ABL(L), 0.10);
        Q.shapeBg = [sh.bgL sh.bgR];
        [Q.FluxShape25Seg, ~, sh] = shape_flux(leg, Q.k1(1), Q.k2(1), D.Density(L), D.MeanWSPD, D.CosTheta(L), D.ABL(L), 0.25);
        Q.shapeBg25 = [sh.bgL sh.bgR];
        [Q.edgeEnh_ppb, Q.outsideFrac] = edge_check(leg, s, Q.k1(1), Q.k2(1));
    end
    % WRF wind at the transect vs the aircraft (ABL mean and at flight level)
    Q.dirWrf = NaN; Q.spdWrf = NaN; Q.dirWrfLev = NaN; Q.spdWrfLev = NaN; Q.pblhWrf = NaN;
    if ~isempty(W)
        xm = nmedian(x0(okp)); ym = nmedian(y0(okp)); tm = nmean(leg.t);
        [uw, vw] = wrf_wind(xm, ym, tm, W);                 [Q.dirWrf, Q.spdWrf] = dir_spd(uw, vw);
        [uz, vz, Q.pblhWrf] = wrf_wind_level(xm, ym, tm, D.LegHeight(L), W);   [Q.dirWrfLev, Q.spdWrfLev] = dir_spd(uz, vz);
    end
    Q.dDirWrf = wrap180(Q.dirWrf - Q.dirAir);  Q.dDirWrfLev = wrap180(Q.dirWrfLev - Q.dirAir);
    Q.dDirLegDay = wrap180(Q.dirAir - D.MeanWDIR);
    % example transect: keep what the figure needs
    if L == exLeg
        Q.example = struct('x', x0, 'y', y0, 's', s, 'ch', leg.ch, 'traj', {ex}, 'cty', cty, ...
            'upx', [], 'upy', [], 'ups', [], 'upch', [], 'uhat', uhat, 'pO', A.pO, 'e', A.e);
        if D.hasUp, Q.example.upx = R.up.x; Q.example.upy = R.up.y; Q.example.ups = R.up.s; Q.example.upch = R.up.ch; end
    end
    if isempty(R.leg), R.leg = Q; else, R.leg(L) = Q; end
end
R.homog = struct('exists', false);
if ~isempty(W), R.homog = wrf_homogeneity(W, D, P); end
end


%% ============================================================ inputs of one day
function D = load_flight(Date, LEG_DIR, SOL_ROOT)
L = load(fullfile(LEG_DIR, Date));
x = load(['INFLUX_' Date]);                       % ASCII file on the path, as in the Method scripts
x(x < -998) = NaN;
D.Date = Date;
D.Time = x(:,1); D.Lat = x(:,2); D.Lon = x(:,3); D.Alt = x(:,4);
WSPD = x(:,11); WDIR = x(:,12);
D.U = -WSPD.*sind(WDIR);  D.V = -WSPD.*cosd(WDIR);
D.CH4 = mmean_omitnan(x(:,6), 30);
if isfield(L, 'CH4Flagged'), D.CH4 = L.CH4Flagged(:); end
D.LegList = L.LegList;
D.hasUp = isfield(L, 'UpwindLeg');
if D.hasUp, D.UpLeg = L.UpwindLeg(1):L.UpwindLeg(2); end

SA = load(fullfile(SOL_ROOT, Date, 'Solution-MethodA.mat'));
D.MeanWSPD = SA.MeanWSPD; D.MeanWDIR = SA.MeanWDIR;
D.MeanU = -SA.MeanWSPD*sind(SA.MeanWDIR);  D.MeanV = -SA.MeanWSPD*cosd(SA.MeanWDIR);
D.Density = SA.Density(:); D.ABL = SA.DownwindABL(:); D.CosTheta = SA.CosTheta(:);
D.FluxA = SA.Flux(:); D.DistA = SA.Distance(:); D.LegTime = SA.DownwindTime(:); D.LegHeight = SA.LegHeight(:);
D.ProfileTime = SA.ProfileTime(:); D.ABLHeight = SA.ABLHeight(:);
[D.ProfileTime, io] = unique(D.ProfileTime); D.ABLHeight = D.ABLHeight(io);
fD = fullfile(SOL_ROOT, Date, 'Solution-MethodD.mat');
D.hasD = D.hasUp && isfile(fD);
if D.hasD
    SD = load(fD); D.CF = SD.CorrectionFactor(:); D.CF2 = SD.CorrectionFactor2(:); D.FluxD = SD.Flux(:);
end
% the wind measured on each transect, as a function of time (model 'obs')
segs = {}; if D.hasUp, segs{end+1} = D.UpLeg; end
for i = 1:size(L.LegList,1), segs{end+1} = L.LegList(i,1):L.LegList(i,2); end %#ok<AGROW>
tt = NaN(numel(segs),1); uu = tt; vv = tt;
for i = 1:numel(segs)
    tt(i) = nmean(D.Time(segs{i})); uu(i) = nmean(D.U(segs{i})); vv(i) = nmean(D.V(segs{i}));
end
[tt, io] = unique(tt); D.obsT = tt; D.obsU = uu(io); D.obsV = vv(io);
end


%% ============================================================ WRF winds of one day
function W = load_wrf_day(D, P)
W = [];
tAll = []; segs = D.LegList;
for i = 1:size(segs,1), tAll = [tAll; D.Time(segs(i,1):segs(i,2))]; end %#ok<AGROW>
if D.hasUp, tAll = [tAll; D.Time(D.UpLeg)]; end
h0 = floor((nmin(tAll) - P.TMAX)/3600); h1 = ceil(nmax(tAll)/3600);
hours = h0:h1; nH = numel(hours);
yyyy = D.Date(1:4); mm = D.Date(5:6); dd = D.Date(7:8);
fns = cell(1,nH); have = false(1,nH);
for q = 1:nH
    if hours(q) < 0 || hours(q) > 23, continue; end
    fns{q} = find_wrfout(P.WRF_ROOT, P.DOMAIN, yyyy, mm, dd, hours(q));
    have(q) = ~isempty(fns{q});
end
if ~any(have), return; end
fn0 = fns{find(have,1)};
XLAT = double(ncread(fn0,'XLAT'));  XLAT = XLAT(:,:,1);
XLON = double(ncread(fn0,'XLONG')); XLON = XLON(:,:,1);
HGT  = double(ncread(fn0,'HGT'));   HGT  = HGT(:,:,1);
x = (XLON - P.LONC)*111.32*cosd(P.LATC); y = (XLAT - P.LATC)*110.57;
in = abs(x) <= P.BOX_KM & abs(y) <= P.BOX_KM;
[ii, jj] = find(in); i1 = min(ii); i2 = max(ii); j1 = min(jj); j2 = max(jj);
ni = i2 - i1 + 1; nj = j2 - j1 + 1;
xs = x(i1:i2, j1:j2); ys = y(i1:i2, j1:j2); hgt = HGT(i1:i2, j1:j2);
[I, J] = ndgrid(1:ni, 1:nj);
Amat = [ones(ni*nj,1) xs(:) ys(:)];
ai = Amat \ I(:); aj = Amat \ J(:);
resid = max([abs(Amat*ai - I(:)); abs(Amat*aj - J(:))]);
try
    ca = double(ncread(fn0,'COSALPHA')); ca = ca(i1:i2, j1:j2, 1);
    sa = double(ncread(fn0,'SINALPHA')); sa = sa(i1:i2, j1:j2, 1);
catch
    ca = ones(ni,nj); sa = zeros(ni,nj);
end
K = P.KMAX;
Uabl = NaN(ni,nj,nH); Vabl = Uabl; PBLH = Uabl; U3 = NaN(ni,nj,K,nH); V3 = U3; Z3 = U3; zi = NaN(1,nH);
for q = find(have)
    fn = fns{q};
    U = double(ncread(fn,'U',[i1 j1 1 1],[ni+1 nj K 1])); U = 0.5*(U(1:end-1,:,:) + U(2:end,:,:));
    V = double(ncread(fn,'V',[i1 j1 1 1],[ni nj+1 K 1])); V = 0.5*(V(:,1:end-1,:) + V(:,2:end,:));
    Ue = U.*ca - V.*sa;  Ve = V.*ca + U.*sa;            % earth-relative
    zw = (double(ncread(fn,'PH',[i1 j1 1 1],[ni nj K+1 1])) + double(ncread(fn,'PHB',[i1 j1 1 1],[ni nj K+1 1])))/9.81;
    zagl = 0.5*(zw(:,:,1:end-1) + zw(:,:,2:end)) - hgt;
    dz = diff(zw, 1, 3);
    zi(q) = abl_at(D, hours(q)*3600);
    w = dz .* (zagl < zi(q));
    w1 = w(:,:,1); w1(sum(w,3) == 0) = 1; w(:,:,1) = w1;  % at least the lowest level
    Uabl(:,:,q) = sum(Ue.*w,3)./sum(w,3);  Vabl(:,:,q) = sum(Ve.*w,3)./sum(w,3);
    U3(:,:,:,q) = Ue; V3(:,:,:,q) = Ve; Z3(:,:,:,q) = zagl;
    PBLH(:,:,q) = double(ncread(fn,'PBLH',[i1 j1 1],[ni nj 1]));
end
% hours without a file: interpolate between neighbours or copy the nearest
filled = {};
for q = find(~have)
    before = find(have(1:q-1), 1, 'last'); after = find(have(q+1:end), 1, 'first');
    if ~isempty(after), after = after + q; end
    if ~isempty(before) && ~isempty(after)
        wq = (q - before)/(after - before);
        Uabl(:,:,q) = (1-wq)*Uabl(:,:,before) + wq*Uabl(:,:,after);  Vabl(:,:,q) = (1-wq)*Vabl(:,:,before) + wq*Vabl(:,:,after);
        U3(:,:,:,q) = (1-wq)*U3(:,:,:,before) + wq*U3(:,:,:,after);  V3(:,:,:,q) = (1-wq)*V3(:,:,:,before) + wq*V3(:,:,:,after);
        Z3(:,:,:,q) = Z3(:,:,:,before); PBLH(:,:,q) = (1-wq)*PBLH(:,:,before) + wq*PBLH(:,:,after); zi(q) = abl_at(D, hours(q)*3600);
        filled{end+1} = sprintf('%02d(interp)', hours(q)); %#ok<AGROW>
    else
        src = before; if isempty(src), src = after; end
        Uabl(:,:,q) = Uabl(:,:,src); Vabl(:,:,q) = Vabl(:,:,src); U3(:,:,:,q) = U3(:,:,:,src); V3(:,:,:,q) = V3(:,:,:,src);
        Z3(:,:,:,q) = Z3(:,:,:,src); PBLH(:,:,q) = PBLH(:,:,src); zi(q) = zi(src);
        filled{end+1} = sprintf('%02d(copy of %02d)', hours(q), hours(src)); %#ok<AGROW>
    end
end
W = struct('tHours', hours*3600, 'hours', hours, 'ai', ai, 'aj', aj, 'ni', ni, 'nj', nj, 'xs', xs, 'ys', ys, ...
    'Uabl', Uabl, 'Vabl', Vabl, 'U3', U3, 'V3', V3, 'Z3', Z3, 'PBLH', PBLH, 'zi', zi, 'have', have, 'files', {fns});
[~, runDir] = fileparts(fileparts(fn0));
W.note = sprintf('WRF %s, hours %02d-%02d UTC (%d files, grid fit residual %.2f cells)', runDir, hours(1), hours(end), nnz(have), resid);
if ~isempty(filled), W.note = [W.note ', filled: ' strjoin(filled, ' ')]; end
end

function zi = abl_at(D, t)
if numel(D.ProfileTime) < 2, zi = D.ABLHeight(1); return; end
t = min(max(t, D.ProfileTime(1)), D.ProfileTime(end));
zi = interp1(D.ProfileTime, D.ABLHeight, t);
end

function fn = find_wrfout(wrfRoot, dom, yyyy, mm, dd, hh)
    % same lookup as Figure7_v4.m, with the file list indexed once
    persistent keys names root
    if isempty(keys) || ~strcmp(root, wrfRoot)
        root = wrfRoot; keys = {}; names = {};
        pat = sprintf('wrfout_d0%d_(\\d{4}-\\d{2}-\\d{2})_(\\d{2})', dom);
        dirs = dir(fullfile(wrfRoot,'WRF_*')); dirs = dirs([dirs.isdir]);
        for i = 1:numel(dirs)
            L = dir(fullfile(wrfRoot, dirs(i).name, sprintf('wrfout_d0%d_*', dom)));
            for j = 1:numel(L)
                tok = regexp(L(j).name, pat, 'tokens','once');
                if ~isempty(tok)
                    keys{end+1} = [tok{1} '_' tok{2}]; names{end+1} = fullfile(wrfRoot, dirs(i).name, L(j).name); %#ok<AGROW>
                end
            end
        end
        fprintf('Indexed %d wrfout_d0%d files under %s\n', numel(names), dom, wrfRoot);
    end
    key = sprintf('%s-%s-%s_%02d', yyyy, mm, dd, hh);
    cands = names(strcmp(keys, key));
    fn = '';
    if isempty(cands), return; end
    isdi = endsWith(cands, '_di');
    if any(isdi), fn = cands{find(isdi,1)}; else, fn = cands{1}; end
end


%% ============================================================ wind models
function [u, v] = const_wind(x, u0, v0)
u = u0*ones(size(x)); v = v0*ones(size(x));
end

function [u, v] = obs_wind(t, T, U, V)
if numel(T) < 2, u = U(1)*ones(size(t)); v = V(1)*ones(size(t)); return; end
tt = min(max(t, T(1)), T(end));
u = interp1(T, U, tt); v = interp1(T, V, tt);
end

function [u, v] = wrf_wind(x, y, t, W)
% ABL-mean WRF wind at (x, y) km and time t (s UTC), bilinear in space, linear in time
tt = min(max(t, W.tHours(1)), W.tHours(end));
q = min(floor((tt - W.tHours(1))/3600) + 1, numel(W.tHours) - 1);
wq = (tt - W.tHours(q))/3600;
i = W.ai(1) + W.ai(2)*x + W.ai(3)*y;  j = W.aj(1) + W.aj(2)*x + W.aj(3)*y;
u = NaN(size(x)); v = u;
for qq = unique(q(:))'
    sel = q == qq;
    u(sel) = (1 - wq(sel)).*interp2(W.Uabl(:,:,qq), j(sel), i(sel)) + wq(sel).*interp2(W.Uabl(:,:,qq+1), j(sel), i(sel));
    v(sel) = (1 - wq(sel)).*interp2(W.Vabl(:,:,qq), j(sel), i(sel)) + wq(sel).*interp2(W.Vabl(:,:,qq+1), j(sel), i(sel));
end
end

function [u, v, pblh] = wrf_wind_level(x, y, t, z, W)
% WRF wind at one point and height z (m AGL), interpolated in the vertical and in time
tt = min(max(t, W.tHours(1)), W.tHours(end));
q = min(floor((tt - W.tHours(1))/3600) + 1, numel(W.tHours) - 1);  wq = (tt - W.tHours(q))/3600;
i = W.ai(1) + W.ai(2)*x + W.ai(3)*y;  j = W.aj(1) + W.aj(2)*x + W.aj(3)*y;
uu = NaN(1,2); vv = uu; pp = uu;
for h = 1:2
    K = size(W.U3,3); up = NaN(K,1); vp = up; zp = up;
    for lev = 1:K
        up(lev) = interp2(W.U3(:,:,lev,q+h-1), j, i); vp(lev) = interp2(W.V3(:,:,lev,q+h-1), j, i); zp(lev) = interp2(W.Z3(:,:,lev,q+h-1), j, i);
    end
    zc = min(max(z, zp(1)), zp(end));
    uu(h) = interp1(zp, up, zc); vv(h) = interp1(zp, vp, zc); pp(h) = interp2(W.PBLH(:,:,q+h-1), j, i);
end
u = (1-wq)*uu(1) + wq*uu(2); v = (1-wq)*vv(1) + wq*vv(2); pblh = (1-wq)*pp(1) + wq*pp(2);
end


%% ============================================================ back trajectories
function R = back_traj(x0, y0, t0, windfun, dt, tmax, cty, up, keepEvery)
% backward kinematic trajectories (midpoint scheme), all parcels at once.
% x0, y0 in km, t0 in s; windfun(x, y, t) returns u, v in m/s.
n = numel(x0); x = x0(:); y = y0(:); t = t0(:);
nStep = ceil(tmax/dt);
entered = inpolygon(x, y, cty.x, cty.y);
tIn = NaN(n,1); tIn(entered) = t(entered);     % going backward: first time inside  (= when the parcel left the county)
tOut = NaN(n,1);                                % going backward: first time outside again (= when it entered the county)
inPrev = entered; lost = false(n,1);
sUp = NaN(n,1); tUp = NaN(n,1); okUp = false(n,1);
if ~isempty(up), ndPrev = (x - up.c(1))*up.nv(1) + (y - up.c(2))*up.nv(2); end
keep = keepEvery > 0;
if keep, ks = 1:keepEvery:n; X = NaN(numel(ks), nStep+1); Y = X; X(:,1) = x(ks); Y(:,1) = y(ks); end
for st = 1:nStep
    [u, v] = windfun(x, y, t);
    bad = isnan(u) | isnan(v); lost = lost | bad; u(bad) = 0; v(bad) = 0;
    xm = x - 0.5*dt*u/1000;  ym = y - 0.5*dt*v/1000;
    [um, vm] = windfun(xm, ym, t - 0.5*dt);
    bad = isnan(um) | isnan(vm); lost = lost | bad; um(bad) = 0; vm(bad) = 0;
    xn = x - dt*um/1000;  yn = y - dt*vm/1000;  tn = t - dt;
    xn(lost) = x(lost); yn(lost) = y(lost);
    in = inpolygon(xn, yn, cty.x, cty.y);
    new = in & ~entered;  entered(new) = true;  tIn(new) = tn(new);
    out = ~in & inPrev & isnan(tOut);  tOut(out) = tn(out);
    inPrev = in;
    if ~isempty(up)
        nd = (xn - up.c(1))*up.nv(1) + (yn - up.c(2))*up.nv(2);
        cr = (nd.*ndPrev < 0 | nd == 0) & ~okUp & ~lost;
        if any(cr)
            f = ndPrev(cr) ./ (ndPrev(cr) - nd(cr)); f(isnan(f)) = 1;
            xc = x(cr) + f.*(xn(cr) - x(cr));  yc = y(cr) + f.*(yn(cr) - y(cr));
            sUp(cr) = (xc - up.c(1))*up.e(1) + (yc - up.c(2))*up.e(2);
            tUp(cr) = t(cr) - f*dt;  okUp(cr) = true;
        end
        ndPrev = nd;
    end
    x = xn; y = yn; t = tn;
    if keep, X(:,st+1) = x(ks); Y(:,st+1) = y(ks); end
end
R = struct('entered', entered, 'tIn', tIn, 'tOut', tOut, 'lost', lost, 'sUp', sUp, 'tUp', tUp, 'okUp', okUp);
if keep, R.X = X; R.Y = Y; R.ks = ks; end
end

function A = leg_axis(xk, yk, uhat)
% principal axis of a transect (km), oriented east (or north), with the origin
% where the mean-wind line through the county centre crosses the transect
c = [mean(xk) mean(yk)];
[Vv, Dd] = eig(cov([xk(:) - c(1), yk(:) - c(2)])); [~, mx] = max(diag(Dd)); e = Vv(:,mx)';
if abs(e(1)) >= abs(e(2)), if e(1) < 0, e = -e; end, A.orient = 'E'; else, if e(2) < 0, e = -e; end, A.orient = 'N'; end
nv = [-e(2) e(1)];
M = [e' -uhat'];
if abs(det(M)) < 1e-6, pO = c; else, ab = M \ (-c'); pO = c + ab(1)*e; end
A.c = c; A.e = e; A.nv = nv; A.pO = pO;
A.s = @(x, y) (x - pO(1))*e(1) + (y - pO(2))*e(2);
end


%% ============================================================ fluxes (replicas of the Method scripts)
function [F, info] = edge_flux(leg, k1, k2, rho, Uw, cosT, zi)
% Edge method (MassBalanceMethodA.m) with the city segment leg points k1:k2
n = numel(leg.ch); ch = leg.ch(:);
xx = leg_x(leg);
k1 = max(1, min(k1, n)); k2 = max(k1, min(k2, n));
m10 = round(n/10);
if k1 > 1, bg1 = nmedian(ch(1:k1-1)); else, bg1 = nmedian(ch(1:m10)); end
if k2 < n, bg2 = nmedian(ch(k2+1:n)); else, bg2 = nmedian(ch(n-m10:n)); end
c = polyfit([xx(1) xx(n)], [bg1 bg2], 1); bgLine = c(1)*xx + c(2);
[Enh, Dn, Up] = segment_enhancement(xx, ch, bgLine, k1, k2);
Dist = arcdeg(leg.lat(k1), leg.lon(k1), leg.lat(k2), leg.lon(k2))*111;        % km, as Method A
F = rho*Uw*cosT*(Enh/1e6)*(Dist*1000)*zi;
info = struct('bg1', bg1, 'bg2', bg2, 'Enh', Enh, 'Dist', Dist, 'Dn', Dn, 'Up', Up);
end

function F = upwind_edge_flux(leg, chUp, pair, k1, k2, CF, rho, Uw, cosT, zi)
% Upwind Edge method (MassBalanceMethodD.m): background from the paired upwind
% points, gap filling, lag alignment, 10th-percentile scaling, then the flux
% over the city segment k1:k2.  pair(k) = index into chUp, NaN = no pairing.
n = numel(leg.ch); ch = leg.ch(:); xx = leg_x(leg);
bg = NaN(1,n); ok = ~isnan(pair); bg(ok) = chUp(pair(ok))*CF;
bg = fill_gaps(bg);
Dn = ch;
R = NaN(1,101); cnt = 0;
for Shift = -50:50
    cnt = cnt + 1;
    TB = NaN(n+100,1); TD = NaN(n+100,1);
    TD(51:end-50) = Dn;  TB(51+Shift:end-50+Shift) = bg;
    okc = ~isnan(TB) & ~isnan(TD);
    if nnz(okc) > 2, cc = corrcoef(TB(okc), TD(okc)); R(cnt) = cc(1,2); end
end
if all(isnan(R)), F = NaN; return; end
Shift = find(R == max(R), 1) - 51;
BD = NaN(n+100,1); BD(51+Shift:end-50+Shift) = bg; bg = BD(51:end-50)';
bg = fill_gaps(bg);
bg(isnan(bg)) = nmedian(bg);
Bad = find(~isnan(Dn) & ~isnan(bg'));
TempSpot = max(1, round(numel(Bad)*0.1));
s1 = sort(bg(Bad)); s2 = sort(Dn(Bad));
CF2 = s2(TempSpot)/s1(TempSpot);
bg = bg*CF2;
Enh = segment_enhancement(xx, ch, bg(:), k1, k2);
Dist = arcdeg(leg.lat(k1), leg.lon(k1), leg.lat(k2), leg.lon(k2))*111;
F = rho*Uw*cosT*(Enh/1e6)*(Dist*1000)*zi;
end

function [Fseg, Ffull, info] = shape_flux(leg, k1, k2, rho, Uw, cosT, zi, frac)
% profile-shape background: a straight line through the mean of the lowest
% fraction frac (0.10 or 0.25) of the CH4 values on either side of the plume
% maximum (no geometry at all)
n = numel(leg.ch); ch = leg.ch(:); xx = leg_x(leg);
sm = mmean_omitnan(ch, 60); [~, kp] = nmax(sm);
left = 1:kp-1; right = kp+1:n;
if numel(left) < 10 || numel(right) < 10, Fseg = NaN; Ffull = NaN; info = struct('bgL',NaN,'bgR',NaN); return; end
[bgL, xL] = low_tenth(ch(left), xx(left), frac); [bgR, xR] = low_tenth(ch(right), xx(right), frac);
c = polyfit([xL xR], [bgL bgR], 1); bgLine = c(1)*xx + c(2);
Enh = segment_enhancement(xx, ch, bgLine, k1, k2);
Dist = arcdeg(leg.lat(k1), leg.lon(k1), leg.lat(k2), leg.lon(k2))*111;
Fseg = rho*Uw*cosT*(Enh/1e6)*(Dist*1000)*zi;
Enh = segment_enhancement(xx, ch, bgLine, 1, n);
Dist = arcdeg(leg.lat(1), leg.lon(1), leg.lat(n), leg.lon(n))*111;
Ffull = rho*Uw*cosT*(Enh/1e6)*(Dist*1000)*zi;
info = struct('bgL', bgL, 'bgR', bgR, 'xL', xL, 'xR', xR);
end

function [p, xp] = low_tenth(ch, xx, frac)
ok = ~isnan(ch); ch = ch(ok); xx = xx(ok);
[cs, io] = sort(ch); m = max(1, round(frac*numel(cs)));
p = mean(cs(1:m)); xp = mean(xx(io(1:m)));
end

function [edgeEnh, outsideFrac] = edge_check(leg, s, k1, k2)
% with the Edge background line of Method A: CH4 minus the line within 1 km of
% each edge (ppb), and the share of the positive plume signal of the whole
% transect that lies outside the city segment
n = numel(leg.ch); ch = leg.ch(:); xx = leg_x(leg);
m10 = round(n/10);
if k1 > 1, bg1 = nmedian(ch(1:k1-1)); else, bg1 = nmedian(ch(1:m10)); end
if k2 < n, bg2 = nmedian(ch(k2+1:n)); else, bg2 = nmedian(ch(n-m10:n)); end
c = polyfit([xx(1) xx(n)], [bg1 bg2], 1); dev = ch - (c(1)*xx + c(2));
edgeEnh = 1000*[nmean(dev(abs(s - s(k1)) <= 1)), nmean(dev(abs(s - s(k2)) <= 1))];
pos = max(dev, 0); pos(isnan(pos)) = 0;
inside = false(n,1); inside(k1:k2) = true;
tot = sum(pos); if tot > 0, outsideFrac = sum(pos(~inside))/tot; else, outsideFrac = NaN; end
end

function [Enh, Dn, Up] = segment_enhancement(xx, ch, bgLine, k1, k2)
% mean enhancement over the segment k1:k2 on the 0.001 degree grid of Method A
mid = k1:k2; xm = xx(mid); cm = ch(mid); bm = bgLine(mid); good = ~isnan(cm) & ~isnan(xm);
grid = min(xx([k1 k2])):0.001:max(xx([k1 k2]));
[xu, iu] = unique(xm(good), 'first'); cu = cm(good); cu = cu(iu); bu = bm(good); bu = bu(iu);
if numel(xu) < 2 || numel(grid) < 2, Enh = NaN; Dn = NaN; Up = NaN; return; end
Dn = nmean(interp1(xu, cu, grid)); Up = nmean(interp1(xu, bu, grid));
Enh = Dn - Up;
end

function xx = leg_x(leg)
% Method A integrates along longitude for east-west transects, latitude otherwise
dLat = nmax(leg.lat) - nmin(leg.lat); dLon = (nmax(leg.lon) - nmin(leg.lon))*cosd(nmean(leg.lat));
if dLon > dLat, xx = leg.lon(:); else, xx = leg.lat(:); end
end

function y = fill_gaps(y)
for w = [120 240 360]
    T = mmean_omitnan(y, w); bad = isnan(y); y(bad) = T(bad);
end
end

function d = arcdeg(lat1, lon1, lat2, lon2)
a = sind((lat2 - lat1)/2).^2 + cosd(lat1).*cosd(lat2).*sind((lon2 - lon1)/2).^2;
d = 2*atan2d(sqrt(a), sqrt(1 - a));
end


%% ============================================================ modelled flow around the city (R2.8)
function H = wrf_homogeneity(W, D, P)
H.exists = true;
tmid = nmean(D.LegTime); [~, q] = min(abs(W.tHours - tmid)); H.hour = W.hours(q);
r = hypot(W.xs, W.ys); sel = r <= P.HOMOG_KM;
Ua = W.Uabl(:,:,q); Va = W.Vabl(:,:,q);
[H.dirMean, H.dirSD, H.spdMean, H.spdSD] = circ_stats(Ua(sel), Va(sel));
[uc, vc] = wrf_wind(0, 0, tmid, W); [H.dirCentre, H.spdCentre] = dir_spd(uc, vc);
H.pblhCentre = interp2(W.PBLH(:,:,q), W.aj(1), W.ai(1)); H.ziAircraft = W.zi(q);
% hour by hour at the county centre over the flight
hf = unique(round(D.LegTime/3600))'; hf = hf(hf >= W.hours(1) & hf <= W.hours(end));
H.hoursFlight = hf; H.dirCentreHours = NaN(size(hf)); H.spdCentreHours = NaN(size(hf));
for i = 1:numel(hf)
    [uc, vc] = wrf_wind(0, 0, hf(i)*3600, W); [H.dirCentreHours(i), H.spdCentreHours(i)] = dir_spd(uc, vc);
end
end

function [dirMean, dirSD, spdMean, spdSD] = circ_stats(u, v)
d = mod(atan2d(u, v) + 180, 360); s = hypot(u, v);
C = nmean(cosd(d)); S = nmean(sind(d)); Rbar = hypot(C, S);
dirMean = mod(atan2d(S, C), 360); dirSD = sqrt(-2*log(max(Rbar, 1e-9)))*180/pi;
spdMean = nmean(s); spdSD = std(s(~isnan(s)));
end

function [d, s] = dir_spd(u, v)
d = mod(atan2d(u, v) + 180, 360); s = hypot(u, v);
end

function a = wrap180(a)
a = mod(a + 180, 360) - 180;
end


%% ============================================================ small helpers (toolbox-free)
function m = nmean(v),   v = v(~isnan(v)); if isempty(v), m = NaN; else, m = mean(v); end, end
function m = nmedian(v), v = v(~isnan(v)); if isempty(v), m = NaN; else, m = median(v); end, end
function m = nmin(v),    v = v(~isnan(v)); if isempty(v), m = NaN; else, m = min(v); end, end
function [m, i] = nmax(v)
ok = find(~isnan(v)); if isempty(ok), m = NaN; i = NaN; else, [m, j] = max(v(ok)); i = ok(j); end
end
function y = mmean_omitnan(x, k)
% same as movmean(x, k, 'omitnan'): centred window of k points (for even k, k/2
% points before and k/2-1 after), shrinking at the ends, NaN ignored
sz = size(x); x = x(:); n = numel(x); ok = ~isnan(x); xv = x; xv(~ok) = 0;
cs = [0; cumsum(xv)]; cn = [0; cumsum(double(ok))];
if mod(k,2) == 0, b = k/2; a = k/2 - 1; else, b = (k-1)/2; a = b; end
i = (1:n)'; lo = max(1, i - b); hi = min(n, i + a);
sumv = cs(hi+1) - cs(lo); cntv = cn(hi+1) - cn(lo);
y = sumv./cntv; y(cntv == 0) = NaN; y = reshape(y, sz);
end
function p = pct(v, q)
v = sort(v(~isnan(v))); if isempty(v), p = NaN; return; end
if numel(v) == 1, p = v; return; end
p = interp1(linspace(0, 100, numel(v)), v, q, 'linear');
end


%% ============================================================ collect, print, plot
function collect_and_report(OUT_DIR, DIAG_FILE, P, FN)
files = dir(fullfile(OUT_DIR, '2*.mat')); names = sort({files.name});
if isempty(names), fprintf('No day results in %s\n', OUT_DIR); return; end
nM = numel(P.MODELS); allLeg = struct([]); dayRows = struct([]); EX = []; nUpDays = 0;
for i = 1:numel(names)
    R = load(fullfile(OUT_DIR, names{i}));
    for L = 1:numel(R.leg)
        Q = R.leg(L);
        if isfield(Q, 'example') && ~isempty(Q.example), EX = Q; end
        if isfield(Q, 'example'), Q = rmfield(Q, 'example'); end
        if isempty(allLeg), allLeg = Q; else, allLeg(end+1) = Q; end %#ok<AGROW>
    end
    dr.Date = R.Date; dr.n = numel(R.leg);
    for m = 1:nM
        dr.FluxDay(m) = nmean(arrayfun(@(q) q.Flux(m), R.leg));
        dr.FluxDDay(m) = nmean(arrayfun(@(q) q.FluxD(m), R.leg));
    end
    dr.FluxDayFile = nmean([R.leg.FluxA_file]); dr.FluxDDayFile = nmean([R.leg.FluxD_file]);
    dr.hasD = R.hasD; dr.homog = R.homog; dr.up = R.up; dr.wrfNote = R.wrfNote;
    if isempty(dayRows), dayRows = dr; else, dayRows(end+1) = dr; end %#ok<AGROW>
    nUpDays = nUpDays + R.hasD;
end
nL = numel(allLeg); nD = numel(dayRows);
fprintf('\n%d days, %d transects collected (%d days with Method D)\n', nD, nL, nUpDays);

% diagnostic ratio tags (optional)
ratio = NaN(nD,1);
if isfile(DIAG_FILE)
    try
        Dg = load(DIAG_FILE); Dg = Dg.Diag_All2;
        for i = 1:nD, j = find(Dg(:,5) == str2double(dayRows(i).Date), 1); if ~isempty(j), ratio(i) = Dg(j,3); end, end
    catch, end
end

% ---- per-transect arrays
getm = @(f, m) arrayfun(@(q) q.(f)(m), allLeg)';
FA_file = [allLeg.FluxA_file]'; FD_file = [allLeg.FluxD_file]';
Fl = NaN(nL,nM); FD = Fl; sE1 = Fl; sE2 = Fl; wid = Fl; tNear = Fl; tFar = Fl; lost = Fl; upMed = Fl; upAbs = Fl; contig = Fl;
for m = 1:nM
    Fl(:,m) = getm('Flux', m); FD(:,m) = getm('FluxD', m); contig(:,m) = getm('contig', m);
    sE1(:,m) = arrayfun(@(q) q.sEdge(m,1), allLeg)'; sE2(:,m) = arrayfun(@(q) q.sEdge(m,2), allLeg)';
    wid(:,m) = getm('width_km', m); tNear(:,m) = getm('transitNear_min', m); tFar(:,m) = getm('transitFar_min', m);
    lost(:,m) = getm('lostFrac', m); upMed(:,m) = getm('upShiftMed_km', m); upAbs(:,m) = getm('upShiftAbsMed_km', m);
end
sAn1 = arrayfun(@(q) q.sEdgeAnalytic(1), allLeg)'; sAn2 = arrayfun(@(q) q.sEdgeAnalytic(2), allLeg)';
dDir = [allLeg.dDirWrf]'; dDirLev = [allLeg.dDirWrfLev]'; dLegDay = [allLeg.dDirLegDay]';
spdRatio = [allLeg.spdWrf]'./[allLeg.spdAir]'; spdRatioLev = [allLeg.spdWrfLev]'./[allLeg.spdAir]';
FshS = [allLeg.FluxShapeSeg]'; FshF = [allLeg.FluxShapeFull]'; Fsh25 = [allLeg.FluxShape25Seg]';
edgeEnh = reshape([allLeg.edgeEnh_ppb], 2, nL)'; outFrac = [allLeg.outsideFrac]';
pairF = NaN(nL,nM); for m = 1:nM, pairF(:,m) = getm('pairFrac', m); end
DistFile = [allLeg.DistA_file]';
sLo = [allLeg.sMin]'; sHi = [allLeg.sMax]';
sAn1c = min(max(sAn1, sLo), sHi); sAn2c = min(max(sAn2, sLo), sHi);        % analytic edges clipped to the transect

% ---- CSV, one row per transect
csv = fullfile(OUT_DIR, 'Trajectory_check_legs.csv');
fid = fopen(csv, 'w');
fprintf(fid, 'Date,Leg,HourUTC,nPts,Height_m,ABL_m,cosTheta,dirDay,spdDay,dirLeg,spdLeg,dDirLegDay,dirWrfABL,spdWrfABL,dDirWrfABL,dirWrfLev,spdWrfLev,dDirWrfLev,pblhWrf');
for m = 1:nM, fprintf(fid, ',sEdge1_%s,sEdge2_%s,width_%s,Flux_%s,FluxD_%s,transitNear_%s,transitFar_%s,contig_%s,lost_%s,upShiftMed_%s,upShiftAbs_%s,pairFrac_%s', ...
        P.MODELS{m}, P.MODELS{m}, P.MODELS{m}, P.MODELS{m}, P.MODELS{m}, P.MODELS{m}, P.MODELS{m}, P.MODELS{m}, P.MODELS{m}, P.MODELS{m}, P.MODELS{m}, P.MODELS{m}); end
fprintf(fid, ',sEdge1_analytic,sEdge2_analytic,sMin,sMax,FluxA_file,DistA_file,FluxD_file,FluxShape10Seg,FluxShape25Seg,FluxShape10Full,edgeEnh1_ppb,edgeEnh2_ppb,outsideFrac\n');
for i = 1:nL
    q = allLeg(i);
    fprintf(fid, '%s,%d,%.3f,%d,%.0f,%.0f,%.3f,%.1f,%.2f,%.1f,%.2f,%.1f,%.1f,%.2f,%.1f,%.1f,%.2f,%.1f,%.0f', q.Date, q.Leg, q.hourUTC, q.nPts, q.height_m, q.ABL_m, q.cosTheta, ...
        q.dirDay, q.spdDay, q.dirAir, q.spdAir, q.dDirLegDay, q.dirWrf, q.spdWrf, q.dDirWrf, q.dirWrfLev, q.spdWrfLev, q.dDirWrfLev, q.pblhWrf);
    for m = 1:nM
        fprintf(fid, ',%.2f,%.2f,%.2f,%.2f,%.2f,%.1f,%.1f,%.3f,%.3f,%.2f,%.2f,%.3f', sE1(i,m), sE2(i,m), wid(i,m), Fl(i,m), FD(i,m), tNear(i,m), tFar(i,m), contig(i,m), lost(i,m), upMed(i,m), upAbs(i,m), pairF(i,m));
    end
    fprintf(fid, ',%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.3f\n', sAn1(i), sAn2(i), sLo(i), sHi(i), FA_file(i), DistFile(i), FD_file(i), FshS(i), Fsh25(i), FshF(i), edgeEnh(i,1), edgeEnh(i,2), outFrac(i));
end
fclose(fid);
fprintf('Wrote %s\n', csv);

% ---- printed summary
fprintf('\n================ CHECKS ================\n');
d = Fl(:,1) - FA_file;
fprintf('Edge flux, straight-line replica vs Solution-MethodA: max |diff| %.2f mol/s, %d transects differ by > 1 mol/s\n', nmax(abs(d)), nnz(abs(d) > 1));
fprintf('City-segment width, replica vs Solution-MethodA Distance: max |diff| %.2f km\n', nmax(abs(wid(:,1) - DistFile)));
fprintf('Edges, parcel-based vs analytic projection clipped to the transect (straight line): max |diff| %.2f km; %d edges beyond the end of the transect (Method A uses the end point there)\n', ...
    nmax(abs([sE1(:,1) - sAn1c; sE2(:,1) - sAn2c])), nnz(sAn1 < sLo) + nnz(sAn2 > sHi));
d = FD(:,1) - FD_file;
fprintf('Upwind Edge flux, straight-line replica vs Solution-MethodD: max |diff| %.2f mol/s, %d transects differ by > 1 mol/s, %d by > 5 (n = %d of %d with a file value)\n', ...
    nmax(abs(d)), nnz(abs(d) > 1), nnz(abs(d) > 5), nnz(~isnan(d)), nnz(~isnan(FD_file)));
fprintf('Upwind pairing coverage (share of city-segment points whose back trajectory crosses the upwind transect): straight median %.2f, obs %.2f, wrf %.2f; transects with no pairing: straight %d, obs %d, wrf %d\n', ...
    nmedian(pairF(:,1)), nmedian(pairF(:,3)), nmedian(pairF(:,4)), nnz(~isnan(FD_file) & isnan(FD(:,1))), nnz(~isnan(FD_file) & isnan(FD(:,3))), nnz(~isnan(FD_file) & isnan(FD(:,4))));
fprintf('Transit time, transect to the near county boundary: median %.0f min (IQR %.0f-%.0f); to the far boundary: median %.0f min (IQR %.0f-%.0f)\n', ...
    pct(tNear(:,1),50), pct(tNear(:,1),25), pct(tNear(:,1),75), pct(tFar(:,1),50), pct(tFar(:,1),25), pct(tFar(:,1),75));
fprintf('Transects whose WRF trajectories never reach the county: %d of %d\n', nnz(isnan(Fl(:,4)) & ~isnan(Fl(:,1))), nL);
fprintf('Parcels lost outside the WRF sub-domain: max fraction %.2f\n', nmax(lost(:,4)));

fprintf('\n================ EDGE POSITIONS (R2.7) ================\n');
fprintf('Shift of the city-segment edges relative to the straight-line projection, all transects, both edges (km):\n');
for m = 2:nM
    sh = [sE1(:,m) - sE1(:,1); sE2(:,m) - sE2(:,1)]; a = abs(sh); a = a(~isnan(a));
    fprintf('  %-9s n = %3d  |shift| median %.1f, IQR %.1f-%.1f, 90th pct %.1f, max %.1f km;  > 2 km: %d (%.0f%%), > 5 km: %d (%.0f%%);  signed: median %+.1f km\n', ...
        P.MODELS{m}, numel(a), pct(a,50), pct(a,25), pct(a,75), pct(a,90), nmax(a), nnz(a > 2), 100*mean(a > 2), nnz(a > 5), 100*mean(a > 5), nmedian(sh));
end
fprintf('City-segment width: straight-line median %.1f km; ', pct(wid(:,1),50));
for m = 2:nM, fprintf('%s median %.1f km (ratio to straight: median %.3f); ', P.MODELS{m}, pct(wid(:,m),50), pct(wid(:,m)./wid(:,1),50)); end
fprintf('\n');
fprintf('Trajectory-based segments that are not contiguous (points inside the envelope that miss the county): obs %d, wrf %d transects\n', ...
    nnz(contig(:,3) < 0.999), nnz(contig(:,4) < 0.999));

fprintf('\n================ EDGE FLUX (R2.7) ================\n');
for m = 2:nM
    r = Fl(:,m)./Fl(:,1); ok = ~isnan(r) & abs(Fl(:,1)) > 5; dF = Fl(:,m) - Fl(:,1);
    fprintf('  %-9s per transect: F/F_straight median %.3f, IQR %.3f-%.3f, 5-95%% %.3f-%.3f (|F| > 5, n = %d); F - F_straight: median %+.1f, IQR %+.1f to %+.1f, mean |dF| %.1f mol/s; |dF| > 10: %d\n', ...
        P.MODELS{m}, pct(r(ok),50), pct(r(ok),25), pct(r(ok),75), pct(r(ok),5), pct(r(ok),95), nnz(ok), ...
        nmedian(dF), pct(dF,25), pct(dF,75), nmean(abs(dF)), nnz(abs(dF) > 10));
end
fprintf('Daily means (mean +- SD of the %d daily values, SE in brackets):\n', nD);
FDay = reshape([dayRows.FluxDay], nM, nD)'; FDDay = reshape([dayRows.FluxDDay], nM, nD)';
fprintf('  %-22s %5.1f +- %4.1f (%3.1f)\n', 'Edge, Solution files', nmean([dayRows.FluxDayFile]), std([dayRows.FluxDayFile]), std([dayRows.FluxDayFile])/sqrt(nD));
for m = 1:nM
    v = FDay(:,m); v = v(~isnan(v));
    fprintf('  %-22s %5.1f +- %4.1f (%3.1f)   n = %d;  mean |day - straight| %.1f mol/s, largest %.1f\n', ['Edge, ' P.MODELS{m}], mean(v), std(v), std(v)/sqrt(numel(v)), numel(v), ...
        nmean(abs(FDay(:,m) - FDay(:,1))), nmax(abs(FDay(:,m) - FDay(:,1))));
end
if ~all(isnan(ratio))
    g = ratio > 2;
    for m = [3 4]
        fprintf('  Edge %s, ratio > 2 days: %.1f +- %.1f (n = %d); ratio <= 2 days: %.1f +- %.1f (n = %d)\n', P.MODELS{m}, ...
            nmean(FDay(g,m)), std(FDay(g & ~isnan(FDay(:,m)),m)), nnz(g & ~isnan(FDay(:,m))), nmean(FDay(~g,m)), std(FDay(~g & ~isnan(FDay(:,m)),m)), nnz(~g & ~isnan(FDay(:,m))));
    end
end
fprintf('Days where the daily Edge mean changes by more than 10 mol/s:\n');
for i = 1:nD
    for m = [3 4]
        if abs(FDay(i,m) - FDay(i,1)) > 10
            fprintf('  %s %-8s %6.1f -> %6.1f mol/s (ratio %.2f)\n', dayRows(i).Date, P.MODELS{m}, FDay(i,1), FDay(i,m), ratio(i));
        end
    end
end

fprintf('\n================ WIND: AIRCRAFT vs WRF, HOMOGENEITY (R2.8) ================\n');
a = abs(dLegDay); a = a(~isnan(a));
fprintf('Transect wind direction minus daily mean (aircraft): |diff| median %.0f, 90th pct %.0f deg; within 20 deg: %.0f%% of %d transects\n', pct(a,50), pct(a,90), 100*mean(a <= 20), numel(a));
a = abs(dDir); a = a(~isnan(a));
fprintf('WRF ABL-mean wind minus aircraft wind at the transects: |dir diff| median %.0f, IQR %.0f-%.0f, 90th pct %.0f deg (n = %d); signed median %+.0f; speed ratio WRF/aircraft median %.2f, IQR %.2f-%.2f\n', ...
    pct(a,50), pct(a,25), pct(a,75), pct(a,90), numel(a), nmedian(dDir), pct(spdRatio,50), pct(spdRatio,25), pct(spdRatio,75));
a = abs(dDirLev); a = a(~isnan(a));
fprintf('WRF wind at flight level minus aircraft wind:          |dir diff| median %.0f, IQR %.0f-%.0f, 90th pct %.0f deg; speed ratio median %.2f\n', ...
    pct(a,50), pct(a,25), pct(a,75), pct(a,90), pct(spdRatioLev,50));
fprintf('Transects with |WRF - aircraft direction| > 30 deg: %d of %d; > 45 deg: %d\n', nnz(abs(dDir) > 30), nnz(~isnan(dDir)), nnz(abs(dDir) > 45));
sdS = NaN(nD,1); rngT = NaN(nD,1); dUD = NaN(nD,1); pbl = NaN(nD,2);
for i = 1:nD
    H = dayRows(i).homog;
    if isfield(H,'exists') && H.exists
        sdS(i) = H.dirSD; rngT(i) = nmax(abs(wrap180(H.dirCentreHours - H.dirCentreHours(1))));
        pbl(i,:) = [H.pblhCentre H.ziAircraft];
    end
    u = dayRows(i).up;
    if u.exists && isfield(u,'dirWrf') && ~isnan(u.dirWrf)
        j = strcmp({allLeg.Date}, dayRows(i).Date);
        dUD(i) = nmedian(wrap180([allLeg(j).dirWrf] - u.dirWrf));
    end
end
fprintf('WRF spatial homogeneity (ABL-mean wind within %d km of the county centre, mid-flight hour): circular SD of direction median %.0f deg, IQR %.0f-%.0f, max %.0f (n = %d days)\n', ...
    P.HOMOG_KM, pct(sdS,50), pct(sdS,25), pct(sdS,75), nmax(sdS), nnz(~isnan(sdS)));
fprintf('WRF temporal change at the county centre over the flight hours: |change| median %.0f deg, 90th pct %.0f, max %.0f\n', pct(rngT,50), pct(rngT,90), nmax(rngT));
fprintf('WRF direction at the downwind transects minus at the upwind transect: median %+.0f deg, |diff| median %.0f, max %.0f (n = %d days)\n', nmedian(dUD), pct(abs(dUD),50), nmax(abs(dUD)), nnz(~isnan(dUD)));
fprintf('WRF PBLH at the county centre vs aircraft ABL height (mid-flight): ratio median %.2f, IQR %.2f-%.2f\n', pct(pbl(:,1)./pbl(:,2),50), pct(pbl(:,1)./pbl(:,2),25), pct(pbl(:,1)./pbl(:,2),75));

fprintf('\n================ UPWIND PAIRING (R2.11) ================\n');
for m = 2:nM
    a = upAbs(:,m); a = a(~isnan(a)); sgn = upMed(:,m); sgn = sgn(~isnan(sgn));
    fprintf('  %-9s landing point on the upwind transect vs straight line: per-transect median |shift| median %.1f km, IQR %.1f-%.1f, 90th pct %.1f, max %.1f (n = %d); signed median %+.1f km\n', ...
        P.MODELS{m}, pct(a,50), pct(a,25), pct(a,75), pct(a,90), nmax(a), numel(a), nmedian(sgn));
end
okD = ~isnan(FD(:,1));
fprintf('Upwind Edge flux with trajectory-based edges and pairing (n = %d transects):\n', nnz(okD));
for m = 2:nM
    dF = FD(:,m) - FD(:,1); r = FD(:,m)./FD(:,1); ok = okD & ~isnan(dF) & abs(FD(:,1)) > 5;
    fprintf('  %-9s F - F_straight: median %+.1f, IQR %+.1f to %+.1f, mean |dF| %.1f mol/s; ratio median %.3f (|F| > 5, n = %d)\n', ...
        P.MODELS{m}, nmedian(dF(okD)), pct(dF(okD),25), pct(dF(okD),75), nmean(abs(dF(okD))), pct(r(ok),50), nnz(ok));
end
hasD = [dayRows.hasD];
fprintf('Daily means, Upwind Edge (n = %d days): files %.1f +- %.1f;', nnz(hasD), nmean([dayRows(hasD).FluxDDayFile]), std([dayRows(hasD).FluxDDayFile]));
for m = 1:nM, v = FDDay(hasD,m); v = v(~isnan(v)); fprintf(' %s %.1f +- %.1f (n = %d);', P.MODELS{m}, mean(v), std(v), numel(v)); end
fprintf('\n');

fprintf('\n================ PROFILE SHAPE vs GEOMETRY (R2.7, second part) ================\n');
a = edgeEnh(:); a = a(~isnan(a));
fprintf('CH4 minus the Edge background line within 1 km of the geometric edges: median %+.1f ppb, IQR %+.1f to %+.1f, |dev| 90th pct %.1f ppb; |dev| <= 2 ppb on %.0f%%, <= 5 ppb on %.0f%% of %d edges\n', ...
    nmedian(a), pct(a,25), pct(a,75), pct(abs(a),90), 100*mean(abs(a) <= 2), 100*mean(abs(a) <= 5), numel(a));
fprintf('Share of the positive plume signal (above the Edge line) lying outside the city segment: median %.2f, IQR %.2f-%.2f, 90th pct %.2f; > 0.25 on %d of %d transects\n', ...
    pct(outFrac,50), pct(outFrac,25), pct(outFrac,75), pct(outFrac,90), nnz(outFrac > 0.25), nnz(~isnan(outFrac)));
bgE = cell2mat(arrayfun(@(q) q.bg(1,:), allLeg, 'UniformOutput', false)');
for v = 1:2
    if v == 1, Fs = FshS; bgS = reshape([allLeg.shapeBg], 2, nL)'; nm = 'lowest tenth'; else, Fs = Fsh25; bgS = reshape([allLeg.shapeBg25], 2, nL)'; nm = 'lowest quarter'; end
    ok = ~isnan(Fs) & ~isnan(Fl(:,1));
    dF = Fs - Fl(:,1); cc = corrcoef(Fs(ok), Fl(ok,1)); db = 1000*(bgS - bgE);
    dayShape = NaN(nD,1); for i = 1:nD, dayShape(i) = nmean(Fs(strcmp({allLeg.Date}, dayRows(i).Date))); end
    fprintf('Trough-fit background (%s on either side of the maximum), same edges: background minus edge median %+.1f ppb (IQR %+.1f to %+.1f); F - F_edge median %+.1f, IQR %+.1f to %+.1f mol/s, r = %.2f (n = %d); daily mean %.1f +- %.1f vs Edge %.1f +- %.1f\n', ...
        nm, nmedian(db(:)), pct(db(:),25), pct(db(:),75), nmedian(dF), pct(dF,25), pct(dF,75), cc(1,2), nnz(ok), nmean(dayShape), std(dayShape(~isnan(dayShape))), nmean(FDay(:,1)), std(FDay(:,1)));
end
dF = FshF - Fl(:,1);
fprintf('Trough-fit background over the WHOLE transect (no geometry): F - F_edge median %+.1f, IQR %+.1f to %+.1f mol/s; daily mean %.1f\n', nmedian(dF), pct(dF,25), pct(dF,75), ...
    nmean(arrayfun(@(i) nmean(FshF(strcmp({allLeg.Date}, dayRows(i).Date))), 1:nD)));

% ---- per-day table
fprintf('\n================ PER DAY ================\n');
fprintf('%-9s %5s %2s | Edge daily mean: file straight  obs   wrf | med |edge shift| obs  wrf (km) | WRF-aircraft dir (deg) | WRF\n', 'Date', 'ratio', 'n');
for i = 1:nD
    j = strcmp({allLeg.Date}, dayRows(i).Date);
    shO = abs([sE1(j,3) - sE1(j,1); sE2(j,3) - sE2(j,1)]); shW = abs([sE1(j,4) - sE1(j,1); sE2(j,4) - sE2(j,1)]);
    fprintf('%-9s %5.2f %2d | %6.1f %6.1f %6.1f %6.1f | %5.1f %5.1f | %+5.0f | %s\n', dayRows(i).Date, ratio(i), dayRows(i).n, ...
        dayRows(i).FluxDayFile, FDay(i,1), FDay(i,3), FDay(i,4), nmedian(shO), nmedian(shW), nmedian(dDir(j)), dayRows(i).wrfNote);
end

% ---- figure
try
    make_figure(OUT_DIR, EX, sE1, sE2, Fl, dDir, dDirLev, upMed, P, FN);
catch err
    fprintf('Figure not made: %s\n', err.message);
end
end


%% ============================================================ the SI figure
function make_figure(OUT_DIR, EX, sE1, sE2, Fl, dDir, dDirLev, upMed, P, FN)
cO = [0.10 0.40 0.85]; cW = [0.85 0.20 0.15]; cS = [0 0 0]; cL = [0.45 0.45 0.45];
fig = figure('Color','w','Units','pixels','Position',[60 60 1500 900]);

% (a) map of the example transect
ax = subplot(2,3,1); hold(ax,'on'); box(ax,'on');
if ~isempty(EX)
    E = EX.example;
    patch(ax, E.cty.x, E.cty.y, [0.93 0.93 0.93], 'EdgeColor','k', 'LineWidth',1.2);
    if ~isempty(E.upx), plot(ax, E.upx, E.upy, '-', 'Color', cL, 'LineWidth', 1.5); end
    for m = [3 4]
        if numel(E.traj) >= m && ~isempty(E.traj(m)) && isfield(E.traj(m),'X') && ~isempty(E.traj(m).X)
            cm = cO; if m == 4, cm = cW; end
            nP = min(size(E.traj(m).X, 2), round(150*60/P.DT) + 1);          % first 150 min of each trajectory
            plot(ax, E.traj(m).X(:,1:nP)', E.traj(m).Y(:,1:nP)', '-', 'Color', 0.45*cm + 0.55, 'LineWidth', 0.6);
        end
    end
    % straight-line projection lines from the two envelope corners
    for c = 1:4
        pc = [E.cty.x(c) E.cty.y(c)];
        plot(ax, pc(1) + [0 60]*E.uhat(1), pc(2) + [0 60]*E.uhat(2), '--', 'Color', cS, 'LineWidth', 0.9);
    end
    scatter(ax, E.x, E.y, 14, 1000*(E.ch - nmedian(E.ch)), 'filled');
    colormap(ax, parula); cb = colorbar(ax); cb.Label.String = 'CH_4 minus transect median (ppb)'; cb.FontSize = 9;
    mk = {'s','d','o','^'}; cc = {cS, cL, cO, cW}; nm = {'straight line','transect wind','observed wind (time varying)','WRF trajectories'};
    hh = []; lab = {};
    for m = 1:4
        if ~isnan(EX.k1(m))
            h = plot(ax, E.x([EX.k1(m) EX.k2(m)]), E.y([EX.k1(m) EX.k2(m)]), mk{m}, 'MarkerSize', 9, 'LineWidth', 1.4, 'Color', cc{m}, 'MarkerFaceColor', 'w');
            hh(end+1) = h; lab{end+1} = nm{m}; %#ok<AGROW>
        end
    end
    axis(ax, 'equal'); xlim(ax, [-60 60]); ylim(ax, [-60 60]);
    legend(ax, hh, lab, 'Location','southoutside', 'FontSize', 8, 'Box','off');
    title(ax, sprintf('%s/%s/%s  T%d  %02d:%02d UTC', EX.Date(5:6), EX.Date(7:8), EX.Date(1:4), EX.Leg, floor(EX.hourUTC), round(mod(EX.hourUTC,1)*60)), 'FontSize', 10);
end
xlabel(ax, 'km east of the county centre'); ylabel(ax, 'km north of the county centre');
set(ax, 'FontName', FN, 'FontSize', 10, 'TickDir', 'out');
text(ax, 0.02, 0.97, '(a)', 'Units','normalized', 'FontSize', 12, 'FontWeight','bold', 'VerticalAlignment','top');

% (b) the example transect with the edges
ax = subplot(2,3,2); hold(ax,'on'); box(ax,'on');
if ~isempty(EX)
    E = EX.example;
    plot(ax, E.s, 1000*E.ch, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    yl = [nmin(1000*E.ch) - 5, nmax(1000*E.ch) + 5];
    for m = 1:4
        if ~isnan(EX.k1(m))
            for e2 = 1:2
                plot(ax, EX.sEdge(m,e2)*[1 1], yl, '-', 'Color', cc{m}, 'LineWidth', 1.3 + 0.5*(m==1));
            end
        end
    end
    ylim(ax, yl); xlim(ax, [nmin(E.s) nmax(E.s)]);
end
xlabel(ax, 'Distance along the transect (km)'); ylabel(ax, 'CH_4 (ppb)');
set(ax, 'FontName', FN, 'FontSize', 10, 'TickDir', 'out');
text(ax, 0.02, 0.97, '(b)', 'Units','normalized', 'FontSize', 12, 'FontWeight','bold', 'VerticalAlignment','top');
title(ax, 'edges: black straight line, grey transect wind, blue observed wind, red WRF', 'FontSize', 9, 'FontWeight','normal');

% (c) edge shifts
ax = subplot(2,3,3); hold(ax,'on'); box(ax,'on');
shO = [sE1(:,3) - sE1(:,1); sE2(:,3) - sE2(:,1)]; shW = [sE1(:,4) - sE1(:,1); sE2(:,4) - sE2(:,1)];
edges = -15:1:15;
hist_bar(ax, shO, edges, cO, 0.55); hist_bar(ax, shW, edges, cW, 0.55);
xlabel(ax, 'Edge shift relative to the straight line (km)'); ylabel(ax, 'Number of edges');
legend(ax, {sprintf('observed wind (median |shift| %.1f km)', nmedian(abs(shO))), sprintf('WRF (median |shift| %.1f km)', nmedian(abs(shW)))}, 'Location','northeast', 'FontSize', 8, 'Box','off');
set(ax, 'FontName', FN, 'FontSize', 10, 'TickDir', 'out'); xlim(ax, [-15 15]);
text(ax, 0.02, 0.97, '(c)', 'Units','normalized', 'FontSize', 12, 'FontWeight','bold', 'VerticalAlignment','top');

% (d) fluxes
ax = subplot(2,3,4); hold(ax,'on'); box(ax,'on');
lim = [nmin([Fl(:); 0]) - 10, nmax(Fl(:)) + 10];
plot(ax, lim, lim, '-', 'Color', [0.6 0.6 0.6]);
plot(ax, Fl(:,1), Fl(:,3), 'o', 'Color', cO, 'MarkerSize', 5, 'MarkerFaceColor', cO);
plot(ax, Fl(:,1), Fl(:,4), '^', 'Color', cW, 'MarkerSize', 5);
xlim(ax, lim); ylim(ax, lim); axis(ax, 'square');
xlabel(ax, 'Downwind Edge flux, straight-line edges (mol s^{-1})'); ylabel(ax, 'Downwind Edge flux, trajectory edges (mol s^{-1})');
legend(ax, {'1:1', 'observed wind', 'WRF'}, 'Location','northwest', 'FontSize', 8, 'Box','off');
set(ax, 'FontName', FN, 'FontSize', 10, 'TickDir', 'out');
text(ax, 0.02, 0.97, '(d)', 'Units','normalized', 'FontSize', 12, 'FontWeight','bold', 'VerticalAlignment','top');

% (e) WRF vs aircraft wind direction
ax = subplot(2,3,5); hold(ax,'on'); box(ax,'on');
edges = -90:10:90;
hist_bar(ax, dDir, edges, cW, 0.55); hist_bar(ax, dDirLev, edges, [0.95 0.65 0.2], 0.45);
xlabel(ax, 'WRF minus aircraft wind direction at the transect (deg)'); ylabel(ax, 'Number of transects');
a1 = abs(dDir); a1 = a1(~isnan(a1)); a2 = abs(dDirLev); a2 = a2(~isnan(a2));
legend(ax, {sprintf('ABL mean (median |diff| %.0f deg)', nmedian(a1)), sprintf('at flight level (median |diff| %.0f deg)', nmedian(a2))}, 'Location','northeast', 'FontSize', 8, 'Box','off');
set(ax, 'FontName', FN, 'FontSize', 10, 'TickDir', 'out'); xlim(ax, [-90 90]);
text(ax, 0.02, 0.97, '(e)', 'Units','normalized', 'FontSize', 12, 'FontWeight','bold', 'VerticalAlignment','top');

% (f) landing point on the upwind transect (R2.11)
ax = subplot(2,3,6); hold(ax,'on'); box(ax,'on');
edges = -30:2:30;
hist_bar(ax, upMed(:,3), edges, cO, 0.55); hist_bar(ax, upMed(:,4), edges, cW, 0.55);
xlabel(ax, 'Shift of the upwind landing point vs straight line (km)'); ylabel(ax, 'Number of transects');
legend(ax, {sprintf('observed wind (median |shift| %.1f km)', nmedian(abs(upMed(:,3)))), sprintf('WRF (median |shift| %.1f km)', nmedian(abs(upMed(:,4))))}, 'Location','northeast', 'FontSize', 8, 'Box','off');
set(ax, 'FontName', FN, 'FontSize', 10, 'TickDir', 'out'); xlim(ax, [-30 30]);
text(ax, 0.02, 0.97, '(f)', 'Units','normalized', 'FontSize', 12, 'FontWeight','bold', 'VerticalAlignment','top');

set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile(OUT_DIR, 'FigS6_trajectory_check.png'), '-dpng', '-r300');
fprintf('Saved %s\n', fullfile(OUT_DIR, 'FigS6_trajectory_check.png'));
end

function hist_bar(ax, v, edges, col, alpha)
v = v(~isnan(v)); if isempty(v), return; end
v = min(max(v, edges(1)), edges(end) - 1e-9);
bw = edges(2) - edges(1); nb = numel(edges) - 1;
idx = min(nb, max(1, floor((v - edges(1))/bw) + 1));
cnt = accumarray(idx(:), 1, [nb 1]);
ctr = edges(1:end-1) + bw/2;
bar(ax, ctr, cnt, 1, 'FaceColor', col, 'FaceAlpha', alpha, 'EdgeColor', col);
end

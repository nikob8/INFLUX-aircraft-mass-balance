%% ========================================================================
%  Figure7_v4.m
%  ------------------------------------------------------------------------
%  Revised Figure 7 for the reviewer's request at line 370.
%
%  The original figure is KEPT: three maps (regional background, Marion
%  County, combined) with the eight idealized transects, and three line
%  panels showing the model along the DOWNWIND theoretical transect.
%
%  What is added:
%    * the aircraft flight leg drawn on all three maps
%    * modelled and aircraft mean wind arrows on the maps
%    * the observed aircraft CH4 overlaid on the three line panels, as a
%      reference, so the reader can see the idealized transect and the real
%      flight leg sample the same plume
%    * the aircraft's edge-line background drawn on panels (d) and (f)
%
%  The model still follows the theoretical transect; the aircraft is there
%  for comparison only.  The point of the figure is that the two agree.
%
%  The x-axis is distance in km along the transect (Reviewer 1, R1.6), with
%  the origin where the aircraft mass-balance mean wind through the centre of
%  Marion County crosses the transect and distance increasing eastward (or
%  northward for a north-south transect), the same convention as Figure 2.
%  Model transect and flight leg get the same construction, so they share
%  one axis.
%
%  Requires  <COLO_DIR>/Colocated_YYYYMMDD.mat  from Compare_WRF_Obs_v5
%  plus one wrfout for the maps and the theoretical transect.
%
%  A ready-to-paste caption is printed at the end.
%
%  Nikolay Balashov -- INFLUX aircraft paper, revision
%  ========================================================================

clc; clear; close all;

%% ========================= SETTINGS =====================================

CASE_DATE = '20110301';
LEG       = 'T3';        % aircraft leg used for the overlay
MAP_HOUR  = 17;          % UTC for the maps and the theoretical transect
                         % ('auto' = hour nearest the leg midpoint)
MAP_LEVEL = 8;           % model level, as in the original figure
                         % ('auto' = nearest the leg mean altitude)

DIRECTION = 'auto';      % downwind transect to plot: 'auto' picks it from the
                         % modelled wind exactly as Generate_Upwind_v4 does,
                         % or force one of N NE E SE S SW W NW

COLO_DIR  = fullfile(pwd,'figs_v4');
WRF_ROOT  = '/storage/group/zrb5027/default/INFLUX/nvb5011/';
OBS_DIR   = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/Observations/INFLUX';
SHAPE_PATH= '/storage/work/nvb5011/INFLUX_Paper_2/Thomas_2016/Shapefiles';

DOMAIN = 2;  KMAX = 25;
JJJ    = 10;                       % subset half-width, as in the original
LON_E  = 153-JJJ : 167+JJJ;
LAT_E  = 142-JJJ : 157+JJJ;

LAT_LIM = [39.25 40.2912];
LON_LIM = [-86.80 -85.45];
QUIV_SKIP = 10;

SHOW_OBS_DOTS = true;    % aircraft as dots on the line panels
SHOW_FIG = true;
PNG_DPI  = 400;
OUT_PNG  = sprintf('Figure7_v4_%s.png', CASE_DATE);
LATC = 39.7806;  LONC = -86.1401;   % centre of Marion County (mean of the four corners used in the mass balance)

%% ========================= LOAD =========================================

addpath(SHAPE_PATH);
yyyy = CASE_DATE(1:4); mm = CASE_DATE(5:6); dd = CASE_DATE(7:8);

S  = load(fullfile(COLO_DIR, ['Colocated_' CASE_DATE '.mat']));
CO = S.CO;  MB = S.MB;
k  = find(strcmp({CO.label}, LEG), 1);
if isempty(k), error('Leg %s not in Colocated_%s.mat', LEG, CASE_DATE); end
C  = CO(k);  st = C.stats;

A = readmatrix(fullfile(OBS_DIR, ['INFLUX_' CASE_DATE]), 'FileType','text');
ok = A(:,2) > -90 & A(:,3) > -360 & A(:,6) > -99;
trkLat = A(ok,2);  trkLon = A(ok,3);

if ischar(MAP_HOUR) && strcmpi(MAP_HOUR,'auto')
    MAP_HOUR = round(mean(C.secs)/3600);
end
fn = find_wrfout(WRF_ROOT, DOMAIN, yyyy, mm, dd, MAP_HOUR);
if isempty(fn), error('No wrfout for %s %02d UTC', CASE_DATE, MAP_HOUR); end
H = load_wrf_hour(fn, MAP_HOUR, KMAX);
[~,wdirName] = fileparts(fileparts(fn));

XLAT = double(ncread(fn,'XLAT'));  XLAT = XLAT(:,:,1);
XLON = double(ncread(fn,'XLONG')); XLON = XLON(:,:,1);
[nx,ny] = size(XLAT);

if ischar(MAP_LEVEL) && strcmpi(MAP_LEVEL,'auto')
    zc = squeeze(mean(mean(H.zasl(LON_E,LAT_E,:),1,'omitnan'),2,'omitnan'));
    [~,MAP_LEVEL] = min(abs(zc - st.altASL));
end
zLev = mean(H.zasl(LON_E,LAT_E,MAP_LEVEL),'all','omitnan');
tLev = mean(H.hgt(LON_E,LAT_E),'all','omitnan');

fprintf('Maps: %02d UTC from %s, level %d = %.0f m ASL (%.0f m AGL)\n', ...
    MAP_HOUR, wdirName, MAP_LEVEL, zLev, zLev-tLev);

Fbg  = H.bg(:,:,MAP_LEVEL);   Fcty = H.city(:,:,MAP_LEVEL);
Ftot = Fbg + Fcty;
Uk = H.U(:,:,MAP_LEVEL);      Vk = H.V(:,:,MAP_LEVEL);

LatS = XLAT(LON_E,LAT_E);  LonS = XLON(LON_E,LAT_E);
BgS  = Fbg(LON_E,LAT_E);   CtyS = Fcty(LON_E,LAT_E);

%% =============== THEORETICAL DOWNWIND TRANSECT ==========================

wdModel = mod(atan2d(mean(Uk(LON_E,LAT_E),'all','omitnan'), ...
                     mean(Vk(LON_E,LAT_E),'all','omitnan')) + 180, 360);
wsModel = mean(sqrt(Uk(LON_E,LAT_E).^2 + Vk(LON_E,LAT_E).^2),'all','omitnan');

if strcmpi(DIRECTION,'auto')
    DIRECTION = downwind_name(wdModel);
end
[ti, tj] = downwind_indices(DIRECTION);
lin  = sub2ind(size(BgS), ti, tj);
tLat = LatS(lin);  tLon = LonS(lin);
tBg  = BgS(lin);   tCty = CtyS(lin);
tTot = tBg + tCty;

% common axis: km along the transect (Figure 2 convention).  The origin of
% each transect is where the aircraft mass-balance mean wind, drawn through
% the centre of Marion County, crosses it.
if isfield(MB,'MeanWDIR'), wdRef = MB.MeanWDIR(1); else, wdRef = st.wd_obs; end
[xT, orientT] = km_axis(tLat(:), tLon(:), LATC, LONC, wdRef);
xO            = km_axis(C.lat(:), C.lon(:), LATC, LONC, wdRef, orientT);
xlab = 'Distance along transect (km)';
XV = orientT;
[xT, srt] = sort(xT);
tBg = tBg(srt);  tCty = tCty(srt);  tTot = tTot(srt);

fprintf('Wind %.0f deg -> downwind transect "%s", %d points, km axis increasing %s\n', ...
    wdModel, DIRECTION, numel(xT), XV);

%% =============== AIRCRAFT, PUT ON THE SAME BASELINE =====================

obsPPB  = C.obs_ppm(:)*1000;
edgeObs = obsPPB - C.obs_enh_ppb(:);          % aircraft edge-line background

% one constant: match the aircraft edge background to the modelled regional
% background, in the mean, over the range the two transects share
loX = max(min(xT), min(xO));  hiX = min(max(xT), max(xO));
sel = xO >= loX & xO <= hiX;
if nnz(sel) < 10
    warning('Flight leg and %s transect barely overlap; using full-range means.', DIRECTION);
    sel = true(size(xO));
end
bgAtObs = interp1(xT, tBg, xO(sel), 'linear', 'extrap');
shift   = mean(bgAtObs,'omitnan') - mean(edgeObs(sel),'omitnan');
obsS    = obsPPB  + shift;
edgeS   = edgeObs + shift;
refPPM  = -shift/1000;

fprintf('Aircraft referenced to %.3f ppm; overlap %.2f to %.2f\n', refPPM, loX, hiX);

%% ========================= FIGURE =======================================

COL.obs  = [0.10 0.45 0.85];
COL.bg   = [0.00 0.00 0.00];
COL.city = [0.85 0.15 0.10];
COL.tot  = [0.10 0.25 0.80];
COL.acw  = [1.00 0.80 0.10];

fig = figure('Color','w','Units','pixels','Position',[60 60 1420 760], ...
             'Visible', onoff(SHOW_FIG));

xL = [0.055 0.365 0.675];  wP = 0.265;
yMap = 0.545;  hMap = 0.40;
yLin = 0.095;  hLin = 0.34;

cmax = max([BgS(:); CtyS(:); BgS(:)+CtyS(:)],[],'omitnan');
mapTitle = {sprintf('%s/%s/%s Background %02d UTC CH_4', mm, dd, yyyy, MAP_HOUR), ...
            sprintf('%s/%s/%s City %02d UTC CH_4',       mm, dd, yyyy, MAP_HOUR), ...
            sprintf('%s/%s/%s Combined %02d UTC CH_4',   mm, dd, yyyy, MAP_HOUR)};
mapField = {Fbg, Fcty, Ftot};
axM = gobjects(1,3);

for p = 1:3
    axM(p) = axes('Parent',fig,'Position',[xL(p) yMap wP hMap]);
    worldmap(LAT_LIM, LON_LIM); gridm('off');
    hS = surfm(LatS, LonS, mapField{p}(LON_E,LAT_E));
    zdatam(hS,-100);
    colormap(axM(p), parula(256)); caxis([0 cmax]);
    hold on
    draw_basemap();

    z1 = 1:QUIV_SKIP:nx;  z2 = 1:QUIV_SKIP:ny;
    q = quiverm(XLAT(z1,z2), XLON(z1,z2), Vk(z1,z2), Uk(z1,z2));
    set(q,'Color',[0.92 0.92 0.92],'LineWidth',0.8);

    draw_ideal_transects(LatS, LonS, p == 1);      % labels on the first map only

    % aircraft: whole flight faint, the leg used for the overlay bold
    plotm(trkLat, trkLon, '-', 'Color',[0.75 0.75 0.75], 'LineWidth',0.6);
    plotm(C.lat, C.lon, '-', 'Color','w', 'LineWidth',3.4);
    plotm(C.lat, C.lon, '-', 'Color',COL.obs, 'LineWidth',1.9);

    % aircraft mean wind for the leg, black-outlined so it reads anywhere
    la0 = mean(C.lat); lo0 = mean(C.lon);  L = 0.20;
    for lw = [3.4 1.8]
        qa = quiverm(la0, lo0, -L*cosd(st.wd_obs), -L*sind(st.wd_obs), 0);
        if lw > 2, set(qa,'Color','k','LineWidth',lw);
        else,      set(qa,'Color',COL.acw,'LineWidth',lw); end
    end

    title(mapTitle{p}, 'FontSize',10.5, 'FontWeight','bold');
end

cb = colorbar(axM(3));
axM(3).Position = [xL(3) yMap wP hMap];
cb.Position = [xL(3)+wP+0.015 yMap 0.013 hMap];
cb.Label.String = 'CH_4 enhancement (ppb)';
cb.FontSize = 9.5;  cb.TickDirection = 'out';

% ---------------- (d) background --------------------------------------
ax4 = axes('Parent',fig,'Position',[xL(1) yLin wP hLin]); hold on; box off
plot(xO, edgeS, '--', 'Color',COL.obs, 'LineWidth',1.5);
plot(xT, tBg, '-', 'Color',COL.bg, 'LineWidth',2.2);
title(sprintf('Background of %s Transect', DIRECTION), 'FontSize',11,'FontWeight','bold');
ylabel('CH_4 enhancement (ppb)','FontSize',10.5);

% ---------------- (e) city plume --------------------------------------
ax5 = axes('Parent',fig,'Position',[xL(2) yLin wP hLin]); hold on; box off
if SHOW_OBS_DOTS
    plot(xO, C.obs_enh_ppb, '.', 'MarkerSize',5, 'Color',COL.obs);
end
plot(xT, tCty, '-', 'Color',COL.city, 'LineWidth',2.2);
yline(0,'-','Color',[0.85 0.85 0.85]);
title(sprintf('City Plume of %s Transect', DIRECTION), 'FontSize',11,'FontWeight','bold');

% ---------------- (f) combined ----------------------------------------
ax6 = axes('Parent',fig,'Position',[xL(3) yLin wP hLin]); hold on; box off
if SHOW_OBS_DOTS
    plot(xO, obsS, '.', 'MarkerSize',5, 'Color',COL.obs);
end
plot(xO, edgeS, '--', 'Color',COL.obs, 'LineWidth',1.2);
plot(xT, tTot, '-', 'Color',COL.tot, 'LineWidth',2.2);
title('Combined CH_4', 'FontSize',11,'FontWeight','bold');

% ---- shared cosmetics -------------------------------------------------
yAll = [tBg(:); tCty(:); tTot(:); obsS(:); edgeS(:); C.obs_enh_ppb(:)];
yl = [floor(min(yAll,[],'omitnan')/10)*10, ceil(max(yAll,[],'omitnan')/10)*10];
xl = [min([xT(:); xO(:)]) max([xT(:); xO(:)])];
for ax = [ax4 ax5 ax6]
    set(ax,'FontSize',10.5,'TickDir','out','LineWidth',1.0, ...
           'YLim',yl,'XLim',xl,'Box','off');
    grid(ax,'on'); ax.GridAlpha = 0.12;
    xlabel(ax, xlab, 'FontSize',10.5);
end
set([ax5 ax6],'YTickLabel',[]);

exportgraphics(fig, OUT_PNG, 'Resolution', PNG_DPI);
fprintf('\nWrote %s\n', OUT_PNG);

%% ========================= CAPTION ======================================

if isfield(MB,'MeanWDIR'), mbW = MB.MeanWDIR(1); else, mbW = NaN; end
cap = sprintf([ ...
'Figure 7. WRF-Chem CH4 for %s/%s/%s at %02d UTC, model level %d (%.0f m above ground), ' ...
'broken into three components: regional background only (left), Marion County only ' ...
'(centre), and the two combined (right). White lines on the maps mark the eight ' ...
'idealized upwind/downwind transects, labelled by the direction the wind comes from; ' ...
'the yellow outline is Marion County. Grey arrows are the modelled wind field and the ' ...
'gold arrow is the aircraft mean wind for the leg (%.0f deg at %.1f m/s; the modelled ' ...
'mean is %.0f deg at %.1f m/s%s). The blue line is aircraft leg %s (%s-%s UTC, %.0f m ' ...
'above ground); the faint grey line is the rest of the flight. Bottom panels follow the ' ...
'modelled %s downwind transect (solid) with the aircraft measurements overlaid for ' ...
'reference (blue dots). The horizontal axis is distance along the transect, with the ' ...
'origin where the aircraft mean wind through the center of Marion County crosses each ' ...
'transect and distance increasing toward the %s, as in Figure 2. The ' ...
'dashed blue line is the straight edge-line background fitted to the flight data by the ' ...
'mass balance method. Aircraft CH4 is shown relative to %.3f ppm, chosen so its edge ' ...
'background matches the modelled regional background in the mean over the range the two ' ...
'transects share. The diagnostic ratio for this flight day is %.1f.'], ...
mm, dd, yyyy, MAP_HOUR, MAP_LEVEL, zLev-tLev, ...
st.wd_obs, st.ws_obs, wdModel, wsModel, ...
iff(isnan(mbW), '', sprintf(', mass balance %.0f deg', mbW)), ...
LEG, st.t0, st.t1, st.altAGL, DIRECTION, XV, refPPM, 4.6);

fprintf('\n---------------- caption ----------------\n%s\n', cap);
fid = fopen(sprintf('Figure7_v4_%s_caption.txt', CASE_DATE), 'w');
fprintf(fid, '%s\n', cap);  fclose(fid);


%% ========================= FUNCTIONS ====================================

function nm = downwind_name(wd)
% Downwind transect label for a wind coming FROM wd, matching the
% eight-sector logic in Generate_Upwind_v4.
    if     wd >= 157.5 && wd < 202.5, nm = 'N';
    elseif wd >= 202.5 && wd < 247.5, nm = 'NE';
    elseif wd >= 247.5 && wd < 292.5, nm = 'E';
    elseif wd >= 292.5 && wd < 337.5, nm = 'SE';
    elseif wd >= 337.5 || wd <  22.5, nm = 'S';
    elseif wd >=  22.5 && wd <  67.5, nm = 'SW';
    elseif wd >=  67.5 && wd < 112.5, nm = 'W';
    else,                             nm = 'NW';
    end
end

function [ii, jj] = downwind_indices(nm)
% Index pairs into the LON_E x LAT_E subset for each downwind transect,
% exactly the lines Generate_Upwind_v4 samples.
    switch upper(nm)
        case 'N',  ii = (6:28)';            jj = repmat(28, 23, 1);
        case 'S',  ii = (6:28)';            jj = repmat(8, 23, 1);
        case 'E',  ii = repmat(28, 27, 1);  jj = (6:32)';
        case 'W',  ii = repmat(7, 27, 1);   jj = (6:32)';
        case 'NE', i0 = (0:15)';  ii = 33-i0;  jj = 18+i0;
        case 'SW', i0 = (0:13)';  ii = 16-i0;  jj = 3+i0;
        case 'NW', i0 = (0:14)';  ii = 19+i0;  jj = 2+i0;
        case 'SE', i0 = (0:15)';  ii = 2+i0;   jj = 18+i0;
        otherwise, error('Unknown transect direction %s', nm);
    end
end

function draw_ideal_transects(L, N, withLabels)
    w = {'Color','w','LineWidth',1.0};
    plotm(L(6:28,8),  N(6:28,8),  w{:});
    plotm(L(6:28,28), N(6:28,28), w{:});
    plotm([L(16,3) L(3,16)],   [N(16,3) N(3,16)],   w{:});
    plotm([L(33,18) L(18,33)], [N(33,18) N(18,33)], w{:});
    plotm(L(7,6:32),  N(7,6:32),  w{:});
    plotm(L(28,6:32), N(28,6:32), w{:});
    plotm([L(2,18) L(17,33)],  [N(2,18) N(17,33)],  w{:});
    plotm([L(19,2) L(33,16)],  [N(19,2) N(33,16)],  w{:});
    if ~withLabels, return; end
    t = {'Color','w','FontSize',9.5,'FontWeight','bold','HorizontalAlignment','center'};
    textm(39.98,-86.10,'N',t{:});  textm(39.44,-86.10,'S',t{:});
    textm(39.75,-85.72,'E',t{:});  textm(39.75,-86.60,'W',t{:});
    textm(39.85,-85.60,'NE',t{:}); textm(39.39,-86.40,'SW',t{:});
    textm(40.15,-86.35,'NW',t{:}); textm(39.60,-85.60,'SE',t{:});
end

function [s, orient] = km_axis(lat, lon, latC, lonC, wdFrom, orient)
% Signed distance (km) along the principal axis of a transect.  s = 0 where
% the line through (latC, lonC) in the direction the wind blows (wind FROM
% wdFrom) crosses that axis; s increases eastward ('east') or northward
% ('north'); the orientation is chosen from the transect unless given.
    x = (lon - lonC)*111.32*cosd(latC);  y = (lat - latC)*110.57;
    ok = ~isnan(x) & ~isnan(y);
    p0 = [mean(x(ok)) mean(y(ok))];
    [V, D] = eig(cov([x(ok)-p0(1), y(ok)-p0(2)]));  [~, m] = max(diag(D));  e = V(:,m)';
    if nargin < 6 || isempty(orient)
        if abs(e(1)) >= abs(e(2)), orient = 'east'; else, orient = 'north'; end
    end
    if strcmpi(orient,'east')  && e(1) < 0, e = -e; end
    if strcmpi(orient,'north') && e(2) < 0, e = -e; end
    u  = [-sind(wdFrom), -cosd(wdFrom)];        % direction the wind blows toward
    ab = [e', -u'] \ (-p0');                     % p0 + a*e = t*u  ->  a, t
    pO = p0 + ab(1)*e;                          % crossing point (km from the county centre)
    s  = (x - pO(1))*e(1) + (y - pO(2))*e(2);
end

function v = iff(c, a, b)
    if c, v = a; else, v = b; end
end

function v = onoff(tf)
    if tf, v = 'on'; else, v = 'off'; end
end

function fn = find_wrfout(wrfRoot, dom, yyyy, mm, dd, hh)
    key = sprintf('%s-%s-%s_%02d', yyyy, mm, dd, hh);
    pat = sprintf('wrfout_d0%d_(\\d{4}-\\d{2}-\\d{2})_(\\d{2})', dom);
    dirs = dir(fullfile(wrfRoot,'WRF_*')); dirs = dirs([dirs.isdir]);
    fn = '';  cands = {};
    for i = 1:numel(dirs)
        L = dir(fullfile(wrfRoot, dirs(i).name, sprintf('wrfout_d0%d_*', dom)));
        for j = 1:numel(L)
            tok = regexp(L(j).name, pat, 'tokens','once');
            if ~isempty(tok) && strcmp([tok{1} '_' tok{2}], key)
                cands{end+1} = fullfile(wrfRoot, dirs(i).name, L(j).name); %#ok<AGROW>
            end
        end
    end
    if isempty(cands), return; end
    isdi = endsWith(cands,'_di');
    if any(isdi), fn = cands{find(isdi,1)}; else, fn = cands{1}; end
end

function H = load_wrf_hour(fn, hh, KMAX)
    rd3 = @(v,nk) double(ncread(fn, v, [1 1 1 1], [Inf Inf nk 1]));
    bg = 0;  for t = 1:5, bg = bg + rd3(sprintf('tracer_%d',t), KMAX); end
    bg = (bg - 5*300)/1e6;
    ct = 0;  for t = 6:8, ct = ct + rd3(sprintf('tracer_%d',t), KMAX); end
    ct = (ct - 3*300)/1e6;
    U = rd3('U',KMAX);  V = rd3('V',KMAX);
    U = 0.5*(U(1:end-1,:,:)+U(2:end,:,:));  V = 0.5*(V(:,1:end-1,:)+V(:,2:end,:));
    zw = (rd3('PH',KMAX+1)+rd3('PHB',KMAX+1))/9.81;
    zasl = 0.5*(zw(:,:,1:end-1)+zw(:,:,2:end));
    hgt = double(ncread(fn,'HGT')); hgt = hgt(:,:,1);
    H = struct('hour',hh,'bg',bg,'city',ct,'U',U,'V',V,'zasl',zasl,'hgt',hgt,'file',fn);
end

function draw_basemap()
    wsv = warning('off','all');
    try
        S  = shaperead('Indianapolis_Interstate_Mileposts.shp');
        ms = defaultm('tranmerc');
        ms.geoid = almanac('earth','grs80','survey foot');
        ms.origin = [37.5 -85.66666666666667 0];
        ms.falseeasting = 328083.3333333333;  ms.falsenorthing = 820208.3333333331;
        ms.scalefactor = 0.9999666666666667;  ms = defaultm(ms);
        [LATr,LONr] = minvtran(ms,[S.X],[S.Y]);
        scatterm(LATr,LONr,12,'.','MarkerEdgeColor',[0.5 0.5 0.5]);
    catch
    end
    try
        MCB = shaperead('Marion_County_Boundary.shp');
        plotm(MCB.Y, MCB.X, 'y', 'LineWidth', 1.8);
    catch
    end
    warning(wsv);
end
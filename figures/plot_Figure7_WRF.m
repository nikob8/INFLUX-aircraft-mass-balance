%% ========================================================================
%  plot_Figure7_WRF.m   (replaces Figure7_v4.m)
%  ------------------------------------------------------------------------
%  Figure 7: WRF-Chem regional background, Marion County and total CH4 on
%  03/01/2011 (maps, a to c) and along the downwind idealized transect selected
%  by the simulated wind, with the aircraft observations of transect T3 (d to f).
%
%  Revision of 2026-11-07: same data and computation as Figure7_v4.m, redrawn in
%  the style of Figures 2, 3 and 5: panel letters and short titles, a key to the
%  map symbols between the two rows, a short legend in every line panel, the
%  diagnostic ratio in panel (f), one sans-serif font, and a shorter caption.
%  Note: the transects are labeled by their side of the county; the N transect is
%  the downwind transect for southerly winds (downwind_name below).
%
%  Requires  <COLO_DIR>/Colocated_YYYYMMDD.mat  from Compare_WRF_Obs_v5
%  plus one wrfout for the maps and the theoretical transect (Mapping Toolbox).
%  Output: Figure7_YYYYMMDD.png (300 dpi); the caption and the numbers quoted in
%  the text are printed before the export.
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
DIAG_RATIO = 4.6;        % diagnostic ratio of 03/01/2011 (Data Set S2), printed in panel (f)
SHOW_FIG = true;
PNG_DPI  = 300;
VECTOR_PDF = false;      % also try a vector PDF (on the PSU server the vector export
                         % of the map axes fails and closes the figure, hence PNG only)
OUT_NAME = sprintf('Figure7_%s', CASE_DATE);   % .png (300 dpi) and .pdf
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
xlab = sprintf('Distance along %s transect (km)', DIRECTION);
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
% Style of Figures 2, 3 and 5: one sans-serif font, panel letters, a key for the map
% symbols between the two rows, and a legend in every line panel.

FS = 10.5;
[FONT, UNI] = pick_font();
set(0, 'DefaultAxesFontName', FONT, 'DefaultTextFontName', FONT);
if UNI
    INT = 'none';  sCH4 = ['CH' char(8324)];
else
    INT = 'tex';   sCH4 = 'CH_4';
end
fprintf('Font used for all text: %s (%s)\n', FONT, ternary(UNI, 'Unicode sub/superscripts', 'TeX sub/superscripts'));

COL.obs  = [0.10 0.45 0.85];     % aircraft
COL.bg   = [0.00 0.00 0.00];     % simulated regional background
COL.city = [0.85 0.15 0.10];     % simulated Marion County enhancement
COL.tot  = [0.00 0.00 0.00];     % simulated total (black: the aircraft data are blue)
COL.acw  = [1.00 0.80 0.10];     % measured mean wind
COL.sim  = [0.92 0.92 0.92];     % simulated wind
COL.trk  = [0.75 0.75 0.75];     % rest of the flight
pal = parula(256);  COL.swatch = pal(1,:);   % dark background for the key symbols

fig = figure('Color','w','Units','pixels','Position',[60 60 1420 880], 'Visible', onoff(SHOW_FIG));
xL = [0.055 0.365 0.675];  wP = 0.265;
yMap = 0.575;  hMap = 0.37;
yLin = 0.075;  hLin = 0.30;

cmax = max([BgS(:); CtyS(:); BgS(:)+CtyS(:)],[],'omitnan');
mapTitle = {'(a) Regional background', '(b) Marion County', '(c) Total'};
mapField = {Fbg, Fcty, Ftot};
axM = gobjects(1,3);
for p = 1:3
    axes('Parent',fig,'Position',[xL(p) yMap wP hMap]);
    axM(p) = worldmap(LAT_LIM, LON_LIM);  gridm('off');
    setm(axM(p), 'FontName', FONT, 'FontSize', FS-1, 'MLabelLocation', 0.5, 'PLabelLocation', 0.5);
    if p > 1, setm(axM(p), 'ParallelLabel', 'off'); end
    hS = surfm(LatS, LonS, mapField{p}(LON_E,LAT_E));
    zdatam(hS,-100);
    colormap(axM(p), parula(256)); caxis(axM(p), [0 cmax]);
    hold on
    draw_basemap();
    z1 = 1:QUIV_SKIP:nx;  z2 = 1:QUIV_SKIP:ny;
    q = quiverm(XLAT(z1,z2), XLON(z1,z2), Vk(z1,z2), Uk(z1,z2));
    set(q,'Color',COL.sim,'LineWidth',0.8);
    draw_ideal_transects(LatS, LonS, p == 1, DIRECTION, FONT);
    plotm(trkLat, trkLon, '-', 'Color',COL.trk, 'LineWidth',0.6);
    plotm(C.lat, C.lon, '-', 'Color','w', 'LineWidth',3.4);
    plotm(C.lat, C.lon, '-', 'Color',COL.obs, 'LineWidth',1.9);
    la0 = mean(C.lat); lo0 = mean(C.lon);  L = 0.20;
    for lw = [3.4 1.8]
        qa = quiverm(la0, lo0, -L*cosd(st.wd_obs), -L*sind(st.wd_obs), 0);
        if lw > 2, set(qa,'Color','k','LineWidth',lw);
        else,      set(qa,'Color',COL.acw,'LineWidth',lw); end
    end
    text(axM(p), 0, 1.05, mapTitle{p}, 'Units','normalized', 'FontSize',FS+1, 'FontName',FONT, 'Interpreter','none');
    fix_map_labels(axM(p), FONT, FS-1);
end
cb = colorbar(axM(3));
set(axM(3), 'Position', [xL(3) yMap wP hMap]);
cbPos = [xL(3)+wP+0.012 yMap 0.012 hMap];
set(cb, 'Position', cbPos, 'FontSize', FS-1, 'FontName', FONT);
try, set(cb, 'TickDirection', 'out'); catch, end
% units above the bar: with the bar placed by hand MATLAB draws the tick labels on its
% left, where a rotated label collided with them; the quantity is named in the caption
annotation(fig, 'textbox', [cbPos(1)-0.03 cbPos(2)+cbPos(4)+0.004 cbPos(3)+0.06 0.03], 'String', 'ppb', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', FS-1, 'FontName', FONT, 'Interpreter', 'none', 'Margin', 1, 'FitBoxToText', 'off');

% ---------------- key to the map symbols --------------------------------
axK = axes('Parent',fig,'Position',[xL(1) 0.445 xL(3)+wP-xL(1) 0.075], 'XLim',[0 1], 'YLim',[0 1], 'Visible','off');
hold(axK,'on');
key = {
  'line',  {'Color','w','LineWidth',0.9},             'Idealized downwind transects'
  'line',  {'Color','w','LineWidth',2.8},             sprintf('%s transect (panels d to f)', DIRECTION)
  'line',  {'Color',[1 1 0],'LineWidth',1.8},         'Marion County'
  'line',  {'Color',[0.5 0.5 0.5],'LineWidth',2.6},   'Interstates'
  'track', {},                                        sprintf('Aircraft transect %s', LEG(2:end))
  'line',  {'Color',COL.trk,'LineWidth',0.8},         'Rest of flight'
  'arrow', {COL.sim, 1.0},                            'Simulated wind'
  'arrow', {COL.acw, 2.0},                            'Measured mean wind'
};
xs = [0.00 0.26 0.52 0.76  0.00 0.26 0.52 0.76];
ys = [0.74 0.74 0.74 0.74  0.24 0.24 0.24 0.24];
wS = 0.032;  hSw = 0.40;
for k = 1:size(key,1)
    x0 = xs(k);  y0 = ys(k);
    patch(axK, [x0 x0+wS x0+wS x0], [y0-hSw/2 y0-hSw/2 y0+hSw/2 y0+hSw/2], COL.swatch, 'EdgeColor','none');
    a = x0 + 0.004;  b = x0 + wS - 0.004;
    switch key{k,1}
        case 'line'
            plot(axK, [a b], [y0 y0], '-', key{k,2}{:});
        case 'track'
            plot(axK, [a b], [y0 y0], '-', 'Color','w', 'LineWidth',3.4);
            plot(axK, [a b], [y0 y0], '-', 'Color',COL.obs, 'LineWidth',1.9);
        case 'arrow'
            c = key{k,2}{1};  lw = key{k,2}{2};
            if lw > 1.5
                plot(axK, [a b-0.004], [y0 y0], '-', 'Color','k', 'LineWidth',lw+1.6);
                plot(axK, b-0.002, y0, '>', 'MarkerSize',7.5, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
            end
            plot(axK, [a b-0.004], [y0 y0], '-', 'Color',c, 'LineWidth',lw);
            plot(axK, b-0.003, y0, '>', 'MarkerSize',5.5, 'MarkerFaceColor',c, 'MarkerEdgeColor',c);
    end
    text(axK, x0+wS+0.006, y0, key{k,3}, 'FontSize',FS-0.5, 'FontName',FONT, 'Interpreter','none', 'VerticalAlignment','middle');
end

% ---------------- line panels ---------------------------------------------
yAll = [tBg(:); tCty(:); tTot(:); obsS(:); edgeS(:); C.obs_enh_ppb(:)];
yl = [floor(min(yAll,[],'omitnan')/10)*10, ceil((max(yAll,[],'omitnan') + 0.25*(max(yAll,[],'omitnan') - min(yAll,[],'omitnan')))/10)*10];   % head room for the legends
xl = [min([xT(:); xO(:)]) max([xT(:); xO(:)])];
ax4 = axes('Parent',fig,'Position',[xL(1) yLin wP hLin]); hold(ax4,'on');
h4b = plot(ax4, xO, edgeS, '--', 'Color',COL.obs, 'LineWidth',1.5);
h4a = plot(ax4, xT, tBg, '-', 'Color',COL.bg, 'LineWidth',2.2);
lg4 = legend(ax4, [h4a h4b], {'WRF-Chem', 'Downwind Edge background'}, 'Location','northeast', 'FontSize',FS-1, 'Interpreter','none');
ax5 = axes('Parent',fig,'Position',[xL(2) yLin wP hLin]); hold(ax5,'on');
plot(ax5, xl, [0 0], '-', 'Color',[0.85 0.85 0.85]);
plot(ax5, xO, C.obs_enh_ppb, '.', 'MarkerSize',5, 'Color',COL.obs);
h5a = plot(ax5, xT, tCty, '-', 'Color',COL.city, 'LineWidth',2.2);
h5p = plot(ax5, NaN, NaN, '.', 'MarkerSize',14, 'Color',COL.obs);       % legend proxy: larger dot
lg5 = legend(ax5, [h5a h5p], {'WRF-Chem', 'Aircraft'}, 'Location','northeast', 'FontSize',FS-1, 'Interpreter','none');
ax6 = axes('Parent',fig,'Position',[xL(3) yLin wP hLin]); hold(ax6,'on');
plot(ax6, xO, obsS, '.', 'MarkerSize',5, 'Color',COL.obs);
h6c = plot(ax6, xO, edgeS, '--', 'Color',COL.obs, 'LineWidth',1.2);
h6a = plot(ax6, xT, tTot, '-', 'Color',COL.tot, 'LineWidth',2.2);
h6p = plot(ax6, NaN, NaN, '.', 'MarkerSize',14, 'Color',COL.obs);       % legend proxy: larger dot
lg6 = legend(ax6, [h6a h6p h6c], {'WRF-Chem', 'Aircraft', 'Downwind Edge background'}, ...
    'Location','northeast', 'FontSize',FS-1, 'Interpreter','none');
for lg = [lg4 lg5 lg6]
    set(lg, 'FontName', FONT);
    try, set(lg, 'ItemTokenSize', [18 10]); catch, end
end
lineTitle = {'(d) Regional background', '(e) Marion County', '(f) Total'};
axL = [ax4 ax5 ax6];
for k = 1:3
    ax = axL(k);
    set(ax,'FontSize',FS,'FontName',FONT,'TickDir','out','LineWidth',1.0,'YLim',yl,'XLim',xl,'Box','off');
    grid(ax,'on');  try, set(ax, 'GridAlpha', 0.12); catch, end
    xlabel(ax, xlab, 'FontSize',FS);
    text(ax, 0, 1.06, lineTitle{k}, 'Units','normalized', 'FontSize',FS+1, 'FontName',FONT, 'Interpreter','none');
end
set([ax5 ax6],'YTickLabel',[]);
ylabel(ax4, [sCH4 ' enhancement (ppb)'], 'Interpreter',INT, 'FontSize',FS);
text(ax6, 0.03, 0.06, sprintf('Diagnostic ratio: %.1f', DIAG_RATIO), 'Units','normalized', ...
    'FontSize',FS-0.5, 'FontName',FONT, 'Interpreter','none');


%% ========================= CAPTION ======================================

cap = sprintf([ ...
'Figure 7. Example of the WRF-Chem CH4 simulations used to characterize background conditions, for %s/%s/%s at ' ...
'%02d UTC and model level %d (approximately %.0f m above ground). (a to c) Simulated CH4 enhancement from emissions ' ...
'outside Marion County (regional background), from Marion County, and their sum, with the eight idealized downwind ' ...
'transects, labeled in (a) by their position relative to the county (the %s transect is downwind for %serly winds). ' ...
'(d to f) The same quantities sampled from the 3 km model grid along the %s transect, which the simulated wind selects ' ...
'on this day, with the aircraft observations of transect %s (%s to %s UTC, approximately %.0f m above ground), as a function of ' ...
'distance along the transect with the origin and sign convention as in Figure 2. Aircraft CH4 is shown relative to ' ...
'%.3f ppm, chosen so that its Downwind Edge background matches the simulated regional background in the mean over ' ...
'the range the two transects share.'], ...
mm, dd, yyyy, MAP_HOUR, MAP_LEVEL, zLev-tLev, DIRECTION, opposite_name(DIRECTION), DIRECTION, LEG(2:end), ...
hhmm(st.t0), hhmm(st.t1), round(st.altAGL/10)*10, refPPM);
fprintf('\n---------------- caption ----------------\n%s\n', cap);
fprintf('measured mean wind %.0f deg at %.1f m/s, simulated mean %.0f deg at %.1f m/s, aircraft transect %s at %.0f m above ground\n', ...
    st.wd_obs, st.ws_obs, wdModel, wsModel, LEG, st.altAGL);

% numbers quoted in Section 2.6.1 and in the response letter, for checking against this run
eA  = C.obs_enh_ppb(:);
eAs = movmean(eA, 15, 'omitnan');                    % aircraft, 15-point running mean
[pS, iS] = max(tCty);  [pA, iA] = max(eA);  [pAs, iAs] = max(eAs);
fprintf(['city plume: simulated peak %.0f ppb at %.1f km; aircraft peak %.0f ppb at %.1f km ' ...
         '(15-point running mean %.0f ppb at %.1f km); simulated minus aircraft position %.1f km\n'], ...
    pS, xT(iS), pA, xO(iA), pAs, xO(iAs), xT(iS) - xO(iAs));
inO = xT >= loX & xT <= hiX;
if ~any(inO), inO = true(size(xT)); end
fprintf('simulated regional background: range %.0f ppb along the whole transect, %.0f ppb over the shared range\n', ...
    max(tBg) - min(tBg), max(tBg(inO)) - min(tBg(inO)));


%% ========================= EXPORT =======================================
% PNG first; then the PDF, as vector graphics if possible (falls back to an image PDF)
drawnow;
exportgraphics(fig, [OUT_NAME '.png'], 'Resolution', PNG_DPI);
fprintf('\nWrote %s.png (%d dpi)\n', OUT_NAME, PNG_DPI);
if VECTOR_PDF
    fprintf('Writing the vector PDF (can take a minute; keep the figure open) ...\n');
    try
        exportgraphics(fig, [OUT_NAME '.pdf'], 'ContentType', 'vector');
        fprintf('Wrote %s.pdf (vector)\n', OUT_NAME);
    catch err
        fprintf('Vector PDF not written (%s)\n', err.message);
        if isgraphics(fig)
            exportgraphics(fig, [OUT_NAME '.pdf'], 'ContentType', 'image', 'Resolution', PNG_DPI);
            fprintf('Wrote %s.pdf as an image (%d dpi)\n', OUT_NAME, PNG_DPI);
        end
    end
end

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

function draw_ideal_transects(L, N, withLabels, selected, font)
% the eight idealized downwind transects (white), the one sampled in panels d to f
% thicker; labels (position relative to the county) on the first map only
    names = {'N','NE','E','SE','S','SW','W','NW'};
    for k = 1:numel(names)
        [ii, jj] = downwind_indices(names{k});
        lin = sub2ind(size(L), ii, jj);
        lw = 0.9;  if strcmpi(names{k}, selected), lw = 2.8; end
        plotm(L(lin), N(lin), '-', 'Color', 'w', 'LineWidth', lw);
    end
    if ~withLabels, return; end
    t = {'Color','w','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','FontName',font};
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
        plotm(LATr, LONr, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 6);   % one object (scatterm makes one per point)
    catch
    end
    try
        MCB = shaperead('Marion_County_Boundary.shp');
        plotm(MCB.Y, MCB.X, 'y', 'LineWidth', 1.8);
    catch
    end
    warning(wsv);
end

function [font, uni] = pick_font()
% One sans-serif font for all text. On Linux the installed families are read from
% fontconfig (fc-list): MATLAB draws TeX-formatted text (m s^{-1}, CH_4) in a serif
% fallback when FontName is only an alias such as Helvetica or Arial. A font that has
% the superscript and subscript characters (FreeSans, DejaVu Sans, Noto Sans) lets the
% units be written as plain Unicode text (uni = true), with no TeX involved.
    font = 'Helvetica';  uni = false;
    if ispc || ismac
        if any(strcmpi(listfonts, 'Arial')), font = 'Arial'; end
        return
    end
    [st, out] = system('fc-list : family');
    if st ~= 0, return; end
    fams = strtrim(strsplit(out, {newline, ','}));
    [st2, out2] = system('fc-list '':charset=b0 b9 207b 2084 2212'' family');
    glyph = {};
    if st2 == 0, glyph = strtrim(strsplit(out2, {newline, ','})); end
    for c = {'FreeSans', 'DejaVu Sans', 'Noto Sans'}
        if any(strcmp(glyph, c{1})), font = c{1}; uni = true; return; end
    end
    for c = {'Liberation Sans', 'Nimbus Sans', 'Nimbus Sans L', 'Arial', 'Helvetica'}
        if any(strcmp(fams, c{1})), font = c{1}; return; end
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function fix_map_labels(ax, font, fs)
% Mapping Toolbox writes the grid labels with TeX (86.25^{\circ} W); rewrite them as
% plain text with the degree sign so that they use the same font as the rest
    h = findall(ax, 'Type', 'text');
    for k = 1:numel(h)
        str = get(h(k), 'String');
        if ischar(str), str = {str}; end
        if ~iscell(str) || ~any(contains(str, '\circ')), continue; end
        str = strrep(str, '^{\circ}', char(176));
        str = strrep(str, '\circ', char(176));
        if numel(str) == 1, str = str{1}; end
        set(h(k), 'String', str, 'Interpreter', 'none', 'FontName', font, 'FontSize', fs);
    end
end

function out = hhmm(t)
% 'HH:MM:SS' (or 'HH:MM') rounded to the minute
    v = sscanf(strrep(char(t), ':', ' '), '%f');
    if numel(v) < 3, v(3) = 0; end
    m = round(v(1)*60 + v(2) + v(3)/60);
    out = sprintf('%02d:%02d', floor(m/60), mod(m, 60));
end

function w = opposite_name(nm)
% wind direction (from) for which the transect on side NM of the county is downwind
    switch upper(nm)
        case 'N',  w = 'south';      case 'S',  w = 'north';
        case 'E',  w = 'west';       case 'W',  w = 'east';
        case 'NE', w = 'southwest';  case 'SW', w = 'northeast';
        case 'NW', w = 'southeast';  case 'SE', w = 'northwest';
        otherwise, w = '?';
    end
end

%% plot_Figure3_5_Upwind.m
%  Figure 3 (Upwind approach: map and downwind transect) and Figure 5 (Upwind,
%  Adjusted ABL Upwind and Upwind Edge backgrounds side by side) for the first
%  downwind transect of 10/04/2012.
%
%  Revision of 2026-11-06: redrawn in the style of the new Figure 2 (map zoomed to the
%  county and both transects, band swept by Marion County along the mean wind shaded,
%  upwind transect outlined in black with dashed lines showing how points of the
%  downwind transect are traced back to it, one sans-serif font with
%  units written without TeX when the font allows, ends outside the county projection
%  shaded gray, dashed backgrounds, softer fills, legends). The computation is unchanged.
%
%  The upwind background is rebuilt EXACTLY as in MassBalanceMethodB.m: every
%  point of the downwind transect is traced along the mean wind to the upwind
%  transect, the gaps are filled with moving means, the line is shifted by the
%  cross-correlation lag (+-50 points) and filled again.  Methods C and D use
%  the same line multiplied by one factor per transect:
%     Adjusted ABL Upwind (C):  x CorrectionFactor    (ABL dilution factor)
%     Upwind Edge         (D):  x CorrectionFactor x CorrectionFactor2
%  Those factors and the flux labels are read from <Date>/Solution-Method{B,C,D}.mat,
%  so the figure shows the Table S1 numbers.  (If CorrectionFactor2 is missing it
%  is recomputed here with the Method D formula.)
%
%  Needs on the path: INFLUX_<Date>, LegAndProfileList/<Date>.mat, InterX.m,
%  alongTrackKm.m, PlotINFLUX.m or the Indianapolis interstate shapefile,
%  <Date>/Solution-MethodB.mat, -MethodC.mat, -MethodD.mat (in one of SOL_ROOTS),
%  Indy_Mass_Balance_Table_v5.csv (optional; flux labels from Table S1).
%  Requires the Mapping Toolbox. Output: Fig3_Upwind_<Date>_T<k>.png/.pdf and
%  Fig5_UpwindVariants_<Date>_T<k>.png/.pdf

clear; close all; clc;

%% ------------------------------------------------------------ settings
Date    = '20121004';
LegLoop = 1;                                  % first downwind transect
MB_ROOT = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';
OBS_DIR = fullfile(MB_ROOT,'AircraftObs');    % folder with INFLUX_<Date>  (change if elsewhere)
addpath(MB_ROOT); if isfolder(OBS_DIR), addpath(OBS_DIR); end
ROADS_CODE = '/storage/work/zrb5027/WRFGUI_V3/Code';                         % has PlotINFLUX.m (Zach)
ROADS_SHP  = '/storage/work/nvb5011/INFLUX_Paper_2/Thomas_2016/Shapefiles';  % Indianapolis_Interstate_Mileposts.shp
if isfolder(ROADS_CODE), addpath(ROADS_CODE); end
if isfolder(ROADS_SHP),  addpath(ROADS_SHP);  end
outName3 = sprintf('Fig3_Upwind_%s_T%d', Date, LegLoop);
outName5 = sprintf('Fig5_UpwindVariants_%s_T%d', Date, LegLoop);

SSLF_LATLON = [39.714691162109375, -86.205429077148438];   % South Side Landfill
% Folders searched, in this order, for <Date>/Solution-MethodB/C/D.mat: your rerun (MassB,
% as for Figure 2), the current folder, and Zach's originals in MB_ROOT. The correction
% factors do not depend on cos theta, and the flux labels are taken from Table S1 when
% TABLE_S1 is found, so the originals give the same figure.
SOL_ROOTS = {'/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB', pwd, MB_ROOT};
TABLE_S1  = 'Indy_Mass_Balance_Table_v5.csv';   % data/ folder of the repository (optional)
SOL_ROOT = '';
for k = 1:numel(SOL_ROOTS)
    if isfile(fullfile(SOL_ROOTS{k}, Date, 'Solution-MethodC.mat')), SOL_ROOT = SOL_ROOTS{k}; break; end
end
assert(~isempty(SOL_ROOT), 'No %s/Solution-MethodC.mat in any folder of SOL_ROOTS; add the folder that holds it.', Date);
fprintf('Solution files read from %s\n', fullfile(SOL_ROOT, Date));

% ---- geometry, identical to the Method scripts ----------------------------
%BottomRight BottomLeft TopLeft TopRight
CountyCorners=[39.63760423971629, -85.95148383144084; 39.63337385377353, -86.32570560518701; 39.92446141560437, -86.32734193809195; 39.9271330073512, -85.95582253907553];

%% ------------------------------------------------------------ load (as Method B)
load(fullfile(MB_ROOT,'LegAndProfileList',Date));          % LegList, ProfileList, UpwindLeg, CH4Flagged?
assert(exist('UpwindLeg','var')==1, '%s has no upwind transect (UpwindLeg) in LegAndProfileList', Date);
UpwindLeg=UpwindLeg(1):UpwindLeg(2);
AllLeg=[UpwindLeg];
for k=1:size(LegList,1)
    AllLeg=[AllLeg,LegList(k,1):LegList(k,2)]; %#ok<AGROW>
end

x=load(['INFLUX_',Date]);
x(x<-998)=NaN;
Time=x(:,1); Lat=x(:,2); Lon=x(:,3); Alt=x(:,4);
WSPD=x(:,11); WDIR=x(:,12);
U=-1*WSPD.*sind(WDIR);
V=-1*WSPD.*cosd(WDIR);
CH4=movmean(x(:,6),30,'omitnan');
if exist('CH4Flagged','var')==1
    CH4=CH4Flagged;
end
MeanU=mean(U(AllLeg),'omitnan');
MeanV=mean(V(AllLeg),'omitnan');
MeanWSPD=(MeanU.^2+MeanV.^2).^.5;
MeanWDIR=atan2d(MeanU,MeanV)+180;

ppbFac = 1; if median(CH4,'omitnan') < 100, ppbFac = 1000; end   % data in ppm -> plot in ppb
CH4ppb = CH4*ppbFac;

%% ------------------------------------------------------------ segmentation (verbatim Method B)
if (MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0)
    PointB1=CountyCorners(3,:);
    PointB2=CountyCorners(1,:);
elseif (MeanU>0 & MeanV<0) | (MeanU<0 & MeanV>0)
    PointB1=CountyCorners(2,:);
    PointB2=CountyCorners(4,:);
end
PointC=[PointB1(1)+MeanV, PointB1(2)+MeanU*(111/86)];
PointA=[PointB1(1)-MeanV, PointB1(2)-MeanU*(111/86)];
Line1a=[PointA(2) PointB1(2) PointC(2);PointA(1) PointB1(1) PointC(1)];
PointC=[PointB2(1)+MeanV, PointB2(2)+MeanU*(111/86)];
PointA=[PointB2(1)-MeanV, PointB2(2)-MeanU*(111/86)];
Line1b=[PointA(2) PointB2(2) PointC(2);PointA(1) PointB2(1) PointC(1)];

Leg=LegList(LegLoop,1):LegList(LegLoop,2);
Line2=[Lon(Leg)';Lat(Leg)'];
Intersection1=flipud(InterX(Line1a,Line2))';
Intersection2=flipud(InterX(Line1b,Line2))';
if ((MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0))
    if isempty(Intersection1)==1, Intersection1=[max(Lat(Leg)) min(Lon(Leg))]; end
    if isempty(Intersection2)==1, Intersection2=[min(Lat(Leg)) max(Lon(Leg))]; end
    if abs(Intersection1(1)-Intersection2(1)) < 0.01
        Upwind1=find(Lon(Leg)<Intersection1(2));
    elseif abs(Intersection1(2)-Intersection2(2)) < 0.01
        Upwind1=find(Lat(Leg)>Intersection1(1));
    else
        Upwind1=find(Lat(Leg)>Intersection1(1) | Lon(Leg)<Intersection1(2));
    end
    if abs(Intersection1(1)-Intersection2(1)) < 0.01
        Upwind2=find(Lon(Leg)>Intersection2(2));
    elseif abs(Intersection2(2)-Intersection2(2)) < 0.01
        Upwind2=find(Lat(Leg)<Intersection2(1));
    else
        Upwind2=find(Lat(Leg)<Intersection2(1) | Lon(Leg)>Intersection2(2));
    end
elseif((MeanU>0 & MeanV<0) | (MeanU<0 & MeanV>0))
    if isempty(Intersection1)==1, Intersection1=[min(Lat(Leg)) min(Lon(Leg))]; end
    if isempty(Intersection2)==1, Intersection2=[max(Lat(Leg)) max(Lon(Leg))]; end
    if abs(Intersection1(1)-Intersection2(1)) < 0.01
        Upwind1=find(Lon(Leg)<Intersection1(2));
    elseif abs(Intersection1(2)-Intersection2(2)) < 0.01
        Upwind1=find(Lat(Leg)<Intersection1(1));
    else
        Upwind1=find(Lat(Leg)<Intersection1(1) | Lon(Leg)<Intersection1(2));
    end
    if abs(Intersection1(1)-Intersection2(1)) < 0.01
        Upwind2=find(Lon(Leg)>Intersection2(2));
    elseif abs(Intersection2(2)-Intersection2(2)) < 0.01
        Upwind2=find(Lat(Leg)>Intersection2(1));
    else
        Upwind2=find(Lat(Leg)>Intersection2(1) | Lon(Leg)>Intersection2(2));
    end
end
Intersection=[Intersection1;Intersection2];
Upwind1=Leg(Upwind1);
Upwind2=Leg(Upwind2);
Temp=ismember(Leg,[Upwind1 Upwind2]);
Downwind=Leg(Temp==0);
isDown = ismember(Leg,Downwind);

%% ------------------------------------------------------------ upwind background (verbatim Method B)
UpwindLine=[Lon(UpwindLeg), Lat(UpwindLeg)];
BackgroundLine=[];
UpX = NaN(length(Leg), 2);
for UpwindLoop=1:length(Leg)
    PointB=[Lat(Leg(UpwindLoop)), Lon(Leg(UpwindLoop))];
    PointC=[PointB(1)+MeanV, PointB(2)+MeanU*(111/86)];
    PointA=[PointB(1)-MeanV, PointB(2)-MeanU*(111/86)];
    Line1A=[PointA(2) PointB(2) PointC(2);PointA(1) PointB(1) PointC(1)];
    UpwindIntersection=flipud(InterX(Line1A,UpwindLine'))';
    if isempty(UpwindIntersection)==1
        BackgroundLine(UpwindLoop)=NaN; %#ok<SAGROW>
        UpX(UpwindLoop,:) = [NaN NaN]; %#ok<SAGROW>
    else
        UpX(UpwindLoop,:) = UpwindIntersection(1,1:2); %#ok<SAGROW>   % [lat lon] where this point traces to the upwind transect
        TempDistance=distance(UpwindIntersection(1),UpwindIntersection(2),UpwindLine(:,2),UpwindLine(:,1));
        GoodSpot=find(TempDistance==min(TempDistance));
        UpwindAligned=UpwindLeg(GoodSpot);
        BackgroundLine(UpwindLoop)=CH4(UpwindAligned); %#ok<SAGROW>
    end
end
nTraced = nnz(~isnan(BackgroundLine));
BackgroundLine = fill_gaps(BackgroundLine);

% cross-correlation shift (+-50 points), as in Method B
DownwindLine=CH4(Leg);
R = NaN(1,101); ShiftCounter=0;
for Shift=-50:50
    ShiftCounter=ShiftCounter+1;
    TempBackground=NaN(length(BackgroundLine)+100,1);
    TempDownwind=NaN(length(BackgroundLine)+100,1);
    TempDownwind(51:end-50)=DownwindLine;
    TempBackground(51+Shift:end-50+Shift)=BackgroundLine;
    TempCorr=corrcoef(TempBackground,TempDownwind,'Rows','complete');
    R(ShiftCounter)=TempCorr(2);
end
Shift=find(R==max(R),1)-51;
BackgroundLineDoppel=NaN(length(BackgroundLine)+100,1);
BackgroundLineDoppel(51+Shift:end-50+Shift)=BackgroundLine;
BackgroundLine=BackgroundLineDoppel(51:end-50)';
BackgroundLine = fill_gaps(BackgroundLine);
BackgroundLine(isnan(BackgroundLine))=median(BackgroundLine,'omitnan');
bgB = BackgroundLine(:)*ppbFac;                 % ppb, one value per Leg point

%% ------------------------------------------------------------ factors and fluxes (saved solutions)
[FluxB, ~ , ~ ] = read_solution(SOL_ROOT, Date, 'B', LegLoop);
[FluxC, CF, ~ ] = read_solution(SOL_ROOT, Date, 'C', LegLoop);
[FluxD, ~ , CF2] = read_solution(SOL_ROOT, Date, 'D', LegLoop);
% flux labels from Table S1 (the published values) when the csv is available
if isfile(TABLE_S1) && ~isempty(which('readtable'))
    T1 = readtable(TABLE_S1);
    if isdatetime(T1.Date), dd = T1.Date == datetime(Date, 'InputFormat', 'yyyyMMdd');
    else, dd = strcmp(cellstr(string(T1.Date)), datestr(datenum(Date,'yyyymmdd'),'yyyy-mm-dd')); end
    r = find(dd & T1.Transect == LegLoop, 1);
    if ~isempty(r)
        fprintf('Solution files: Upwind %.1f, Adjusted ABL Upwind %.1f, Upwind Edge %.1f mol/s; labels from %s\n', FluxB, FluxC, FluxD, TABLE_S1);
        FluxB = T1.Flux_Upwind_mols(r);  FluxC = T1.Flux_ABLAdjust_mols(r);  FluxD = T1.Flux_UpwindEdge_mols(r);
    end
end
if isnan(CF)
    error('CorrectionFactor for transect %d not found in %s/Solution-MethodC.mat; run MassBalanceMethodC.m first.', LegLoop, Date);
end
bgC = bgB*CF;
% Method D factor: 10th percentile of the downwind trace over the 10th percentile of
% the (ABL adjusted) background, both over the points where neither is NaN
BadSpot=find(isnan(DownwindLine)==0 & isnan(bgC)==0);
TempSpot=round(length(BadSpot)*0.1);
TempSort1=sort(bgC(BadSpot)); TempSort2=sort(DownwindLine(BadSpot)*ppbFac);
CF2here=TempSort2(TempSpot)/TempSort1(TempSpot);
if isnan(CF2), CF2 = CF2here; fprintf('CorrectionFactor2 recomputed here: %.4f\n', CF2); end
bgD = bgC*CF2;

%% ------------------------------------------------------------ km axis (as Figure 2)
latO = mean(Intersection(:,1));  lonO = mean(Intersection(:,2));
LatDistance=distance(min(Lat(Leg)),mean(Lon(Leg),'omitnan'),max(Lat(Leg)),mean(Lon(Leg),'omitnan'));
LonDistance=distance(mean(Lat(Leg),'omitnan'),min(Lon(Leg)),mean(Lat(Leg),'omitnan'),max(Lon(Leg)));
if LonDistance>LatDistance, orient='east'; else, orient='north'; end
s = alongTrackKm(Lat(Leg), Lon(Leg), latO, lonO, orient);

%% ------------------------------------------------------------ report
tD = mean(Time(Leg),'omitnan'); tU = mean(Time(UpwindLeg),'omitnan');
fprintf('%s T%d: downwind %02d:%02d UTC (mean MSL altitude %.0f m), upwind %02d:%02d UTC (mean MSL altitude %.0f m)\n', ...
    Date, LegLoop, floor(tD/3600), round(mod(tD,3600)/60), mean(Alt(Leg),'omitnan'), ...
    floor(tU/3600), round(mod(tU,3600)/60), mean(Alt(UpwindLeg),'omitnan'));
fprintf('wind %.0f deg / %.1f m/s; %d of %d downwind points traced to the upwind transect; shift %d points\n', ...
    MeanWDIR, MeanWSPD, nTraced, numel(Leg), Shift);
fprintf('CorrectionFactor (C) = %.4f, CorrectionFactor2 (D) = %.4f (recomputed here: %.4f)\n', CF, CF2, CF2here);
fprintf('Flux labels: Upwind %.0f, Adjusted ABL Upwind %.0f, Upwind Edge %.0f mol/s  (Table S1 for this transect: -40, 4, 34)\n', FluxB, FluxC, FluxD);

%% ------------------------------------------------------------ style (as Figure 2)
C_BAND  = [0.90 0.92 0.95];     % band swept by the county along the wind (map)
C_OUT   = [0.91 0.91 0.91];     % parts of the transect outside the county projection
C_POS   = [0.55 0.78 0.47];     % positive enhancement
C_NEG   = [0.91 0.50 0.45];     % negative enhancement
C_ROAD  = [0.62 0.62 0.62];
C_TRACE = [0.45 0.45 0.45];     % tracing from the downwind to the upwind transect
FS = 11;
[FONT, UNI] = pick_font();
set(0, 'DefaultAxesFontName', FONT, 'DefaultTextFontName', FONT);
if UNI   % plain text with Unicode sub- and superscripts
    INT = 'none';  sCH4 = ['CH' char(8324)];  sPerS = [' s' char(8315) char(185)];  numstr = @(v) strrep(sprintf('%.0f', v), '-', char(8722));
else     % TeX markup
    INT = 'tex';   sCH4 = 'CH_4';               sPerS = ' s^{-1}';                     numstr = @(v) sprintf('%.0f', v);
end
fprintf('Font used for all text: %s (%s)\n', FONT, ternary(UNI, 'Unicode sub/superscripts', 'TeX sub/superscripts'));
txt = struct('FONT', FONT, 'INT', INT, 'sCH4', sCH4, 'sPerS', sPerS, 'FS', FS, 'numstr', numstr, ...
             'C_OUT', C_OUT, 'C_POS', C_POS, 'C_NEG', C_NEG);
trace = CH4ppb(Leg);  trace = trace(:);  s = s(:);
[hhD, mmD] = hhmm(tD);  [hhU, mmU] = hhmm(tU);

%% ------------------------------------------------------------ Figure 3
lim = [min([CH4ppb(Leg); CH4ppb(UpwindLeg)]) max([CH4ppb(Leg); CH4ppb(UpwindLeg)])];
fig3 = figure('Color','w','Units','pixels','Position',[80 80 1250 520]);

% ---- (a) map --------------------------------------------------------------
axes('Parent', fig3, 'Position', [0.085 0.11 0.33 0.79]);
allLat = [CountyCorners(:,1); Lat(Leg); Lat(UpwindLeg)];  allLon = [CountyCorners(:,2); Lon(Leg); Lon(UpwindLeg)];
latlim = [min(allLat) max(allLat)] + [-0.07 0.07];
lonlim = [min(allLon) max(allLon)] + [-0.09 0.09];
ax1 = worldmap(latlim, lonlim);
set(ax1, 'SortMethod', 'childorder');
setm(ax1, 'FontName', FONT, 'FontSize', FS-1, 'FLineWidth', 1, 'GColor', [0.82 0.82 0.82], ...
    'GLineStyle', '-', 'GLineWidth', 0.4, 'MLabelLocation', 0.25, 'PLabelLocation', 0.2, ...
    'MLabelRound', -2, 'PLabelRound', -1);
hold on
la0 = mean(Lat(Leg), 'omitnan');
dwScore = (CountyCorners(:,2) - mean(CountyCorners(:,2)))*cosd(la0)*MeanU + (CountyCorners(:,1) - mean(CountyCorners(:,1)))*MeanV;
[~, iDown] = max(dwScore);                       % most downwind corner of the county
bandLat = [PointB1(1); Intersection1(1,1); Intersection2(1,1); PointB2(1); CountyCorners(iDown,1)];
bandLon = [PointB1(2); Intersection1(1,2); Intersection2(1,2); PointB2(2); CountyCorners(iDown,2)];
hBand = patchm(bandLat, bandLon, C_BAND);
set(hBand, 'EdgeColor', 'none');
draw_roads(ax1, C_ROAD, hBand);
cc = [CountyCorners; CountyCorners(1,:)];
plotm(cc(:,1), cc(:,2), '-', 'Color', [0 0 0], 'LineWidth', 1.3);

% how the Upwind method pairs the transects: points of the downwind transect traced
% back along the mean wind to the upwind transect (every ~5 km within the band)
iD = find(isDown(:) & ~isnan(UpX(:,1)));
if ~isempty(iD)
    sel = iD(unique(round(linspace(1, numel(iD), 8))));
    for k = sel(:)'
        plotm([UpX(k,1) Lat(Leg(k))], [UpX(k,2) Lon(Leg(k))], '--', 'Color', C_TRACE, 'LineWidth', 0.9);
    end
end
plotm(Lat(UpwindLeg), Lon(UpwindLeg), '-', 'Color', [0 0 0], 'LineWidth', 8);   % upwind transect outlined in black
scatterm(Lat(UpwindLeg), Lon(UpwindLeg), 22, CH4ppb(UpwindLeg), 'filled');
scatterm(Lat(Leg), Lon(Leg), 22, trace, 'filled');
caxis(ax1, lim);  colormap(ax1, parula(256));
cb = colorbar(ax1, 'eastoutside');
set(cb, 'FontName', FONT, 'FontSize', FS-1);
try
    set(get(cb, 'Label'), 'String', [sCH4 ' (ppb)'], 'Interpreter', INT, 'FontName', FONT, 'FontSize', FS);
catch
    ylabel(cb, [sCH4 ' (ppb)'], 'Interpreter', INT, 'FontName', FONT, 'FontSize', FS);
end
plotm(SSLF_LATLON(1), SSLF_LATLON(2), '^', 'MarkerSize', 9, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1], 'LineWidth', 1);
textm(SSLF_LATLON(1) - 0.035, SSLF_LATLON(2) + 0.015, 'SSLF', 'FontName', FONT, 'FontSize', FS-1);
occ = [Lat(Leg) Lon(Leg); Lat(UpwindLeg) Lon(UpwindLeg); CountyCorners; SSLF_LATLON];
[aLat, aLon, lLat, lLon] = wind_arrow_place(latlim, lonlim, occ, MeanU, MeanV, MeanWSPD);
Lw = 0.11*diff(latlim);
hq = quiverm(aLat, aLon, Lw*MeanV/MeanWSPD, Lw*MeanU/MeanWSPD/cosd(mean(latlim)), 'k', 0);
set(hq, 'LineWidth', 1.6);
textm(lLat, lLon, [sprintf('%.1f m', MeanWSPD) sPerS], 'Interpreter', INT, 'FontName', FONT, 'FontSize', FS-1, 'HorizontalAlignment', 'center');
text(ax1, 0, 1.06, sprintf('(a) %s/%s/%s: upwind %02d:%02d UTC, downwind transect %d %02d:%02d UTC', ...
    Date(5:6), Date(7:8), Date(1:4), hhU, mmU, LegLoop, hhD, mmD), 'Units', 'normalized', 'FontSize', FS+1, 'FontName', FONT);
fix_map_labels(ax1, FONT, FS-1);

% ---- (b) downwind transect with the Upwind background ----------------------
ax2 = axes('Parent', fig3, 'Position', [0.53 0.14 0.45 0.76]);
ylB = [floor(min([trace; bgB])) - 3, ceil(max([trace; bgB]) + 0.40*(max([trace; bgB]) - min([trace; bgB])))];
hs = draw_profile(ax2, s, trace, bgB, isDown, ylB, FluxB, txt);
ylabel(ax2, [sCH4 ' mole fraction (ppb)'], 'Interpreter', INT, 'FontSize', FS+1);
text(ax2, 0, 1.06, '(b) Upwind background and enhancement', 'Units', 'normalized', 'FontSize', FS+1, 'FontName', FONT);
labs = {sCH4, 'Upwind background', 'Positive enhancement', 'Negative enhancement', 'Outside the county projection'};
lg = legend(ax2, hs, labs(1:numel(hs)), ...
    'Location', 'northwest', 'FontSize', FS-1, 'Interpreter', INT);
set(lg, 'FontName', FONT, 'Box', 'on');
try, set(lg, 'NumColumns', 2, 'ItemTokenSize', [20 12]); catch, end
save_fig(fig3, outName3);

%% ------------------------------------------------------------ Figure 5
fig5 = figure('Color','w','Units','pixels','Position',[80 80 1500 500]);
yAll = [trace; bgB; bgC; bgD];  yl = [floor(min(yAll)) - 2, ceil(max(yAll)) + 2];   % identical limits on all panels
pos = [0.065 0.13 0.275 0.70; 0.385 0.13 0.275 0.70; 0.705 0.13 0.275 0.70];
bgs = {bgB, bgC, bgD};  flx = [FluxB, FluxC, FluxD];
ttl = {'(a) Upwind', '(b) Adjusted ABL Upwind', '(c) Upwind Edge'};
ax = gobjects(1,3);
for k = 1:3
    ax(k) = axes('Parent', fig5, 'Position', pos(k,:));
    hs = draw_profile(ax(k), s, trace, bgs{k}, isDown, yl, flx(k), txt);
    text(ax(k), 0, 1.05, ttl{k}, 'Units', 'normalized', 'FontSize', FS+1, 'FontName', FONT);
end
ylabel(ax(1), [sCH4 ' mole fraction (ppb)'], 'Interpreter', INT, 'FontSize', FS+1);
labs = {sCH4, 'Background of the method', 'Positive enhancement', 'Negative enhancement', 'Outside the county projection'};
lg = legend(ax(1), hs, labs(1:numel(hs)), ...
    'FontSize', FS-1, 'Interpreter', INT, 'Orientation', 'horizontal');
set(lg, 'FontName', FONT, 'Box', 'off', 'Units', 'normalized', 'Position', [0.12 0.93 0.76 0.05]);
try, set(lg, 'NumColumns', 5, 'ItemTokenSize', [20 12]); catch, end
save_fig(fig5, outName5);


%% ============================================================ local functions
function y = fill_gaps(y)
% the three-pass moving-mean gap filling of MassBalanceMethodB.m
    for w = [120 240 360]
        Temp = movmean(y, w, 'omitnan');
        BadSpot = find(isnan(y)==1);
        y(BadSpot) = Temp(BadSpot);
    end
end

function [Flux, CF, CF2] = read_solution(MB_ROOT, Date, Method, LegLoop)
% flux and correction factors for one transect from <Date>/Solution-Method<Method>.mat
    Flux = NaN; CF = NaN; CF2 = NaN;
    f = fullfile(MB_ROOT, Date, ['Solution-Method' Method '.mat']);   % MB_ROOT here is the folder passed in (SOL_ROOT)
    if ~isfile(f), warning('%s not found', f); return; end
    S = load(f);
    if isfield(S,'Flux') && numel(S.Flux) >= LegLoop, Flux = S.Flux(LegLoop); end
    if isfield(S,'CorrectionFactor')  && numel(S.CorrectionFactor)  >= LegLoop, CF  = S.CorrectionFactor(LegLoop);  end
    if isfield(S,'CorrectionFactor2') && numel(S.CorrectionFactor2) >= LegLoop, CF2 = S.CorrectionFactor2(LegLoop); end
end

function hs = draw_profile(ax, s, trace, bg, isDown, yl, fluxVal, txt)
% CH4 along the downwind transect (black), background of the method (dashed), positive
% (green) and negative (red) enhancement over the part directly downwind of Marion
% County, the ends outside the county projection shaded gray, and the flux label
    s = s(:);  trace = trace(:);  bg = bg(:);  isDown = logical(isDown(:));
    hold(ax, 'on');  box(ax, 'on');
    set(ax, 'XLim', [min(s) max(s)], 'YLim', yl);
    hOut = gobjects(0);
    d = diff([0; ~isDown; 0]);  r0 = find(d == 1);  r1 = find(d == -1) - 1;
    for k = 1:numel(r0)
        xa = min(s(r0(k):r1(k)));  xb = max(s(r0(k):r1(k)));
        hOut(end+1) = patch(ax, [xa xb xb xa], [yl(1) yl(1) yl(2) yl(2)], txt.C_OUT, 'EdgeColor', 'none'); %#ok<AGROW>
    end
    xD = s(isDown);  yD = trace(isDown);  bD = bg(isDown);
    y1 = yD;  bad = (y1 - bD) < 0;  y1(bad) = bD(bad);  g = ~isnan(y1) & ~isnan(bD);
    hPos = patch(ax, [xD(g); flipud(xD(g))], [bD(g); flipud(y1(g))], txt.C_POS, 'EdgeColor', 'none');
    y1 = yD;  bad = (y1 - bD) > 0;  y1(bad) = bD(bad);  g = ~isnan(y1) & ~isnan(bD);
    hNeg = patch(ax, [xD(g); flipud(xD(g))], [bD(g); flipud(y1(g))], txt.C_NEG, 'EdgeColor', 'none');
    hBg = plot(ax, s, bg, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.4);
    hTr = plot(ax, s, trace, '-', 'Color', [0 0 0], 'LineWidth', 1.3);
    set(ax, 'FontSize', txt.FS, 'FontName', txt.FONT, 'TickDir', 'out', 'Layer', 'top', 'XLim', [min(s) max(s)], 'YLim', yl);
    xlabel(ax, 'Distance along flight track (km)', 'FontSize', txt.FS+1);
    text(ax, 0.025, 0.06, ['Marion County: ' txt.numstr(fluxVal) ' mol' txt.sPerS], 'Units', 'normalized', ...
        'Interpreter', txt.INT, 'FontSize', txt.FS, 'FontName', txt.FONT);
    hs = [hTr, hBg, hPos, hNeg];
    if ~isempty(hOut), hs = [hs, hOut(1)]; end
end

function save_fig(fig, name)
    if ~isempty(which('exportgraphics'))
        exportgraphics(fig, [name '.png'], 'Resolution', 300);
        exportgraphics(fig, [name '.pdf'], 'ContentType', 'vector');
        fprintf('Saved %s.png / .pdf\n', name);
    else
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, [name '.png'], '-dpng', '-r300');
        fprintf('Saved %s.png\n', name);
    end
end

function [hh, mm] = hhmm(t)
    hh = floor(t/3600);  mm = round((t - 3600*hh)/60);
    if mm == 60, hh = hh + 1; mm = 0; end
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

function draw_roads(ax, col, keep)
% Major Indianapolis highways. Uses Zach's PlotINFLUX.m (the routine
% behind the submitted figures) if it is on the path; otherwise draws the
% interstates from Indianapolis_Interstate_Mileposts.shp (Indiana State Plane
% East, US survey feet). Line and patch objects drawn so far (apart from the map frame,
% the grid, and KEEP) are set to COL; PlotINFLUX draws its roads with markers, which
% keep their colour (black in the published figure).
    if ~isempty(which('PlotINFLUX'))
        PlotINFLUX;
    elseif ~isempty(which('Indianapolis_Interstate_Mileposts.shp'))
        wsv = warning('off','all');
        S  = shaperead('Indianapolis_Interstate_Mileposts.shp');
        ms = defaultm('tranmerc');
        ms.geoid         = almanac('earth','grs80','survey foot');
        ms.origin        = [37.5 -85.66666666666667 0];
        ms.falseeasting  = 328083.3333333333;
        ms.falsenorthing = 820208.3333333331;
        ms.scalefactor   = 0.9999666666666667;
        ms = defaultm(ms);
        [LATr, LONr] = minvtran(ms, [S.X], [S.Y]);
        scatterm(LATr, LONr, 4, col, 'filled');
        warning(wsv);
        return
    else
        warning('Neither PlotINFLUX.m nor Indianapolis_Interstate_Mileposts.shp is on the path: no roads drawn.');
        return
    end
    skipTags = {'Frame', 'Parallel', 'Meridian', 'MLabel', 'PLabel'};
    objs = findall(ax);
    for k = 1:numel(objs)
        h = objs(k);
        if h == ax || any(h == keep(:)) || any(strcmp(get(h, 'Tag'), skipTags)), continue; end
        switch get(h, 'Type')
            case 'line'
                set(h, 'Color', col, 'LineWidth', 0.9);
            case 'patch'
                set(h, 'EdgeColor', col, 'LineWidth', 0.9);
        end
    end
end

function [aLat, aLon, lLat, lLon] = wind_arrow_place(latlim, lonlim, occ, MeanU, MeanV, MeanWSPD)
% arrow centred in the corner of the map farthest from the occupied points OCC
% (n x 2, lat lon); returns the arrow tail and the label position
    la0 = mean(latlim);  dLat = diff(latlim);  dLon = diff(lonlim);
    cand = [latlim(1)+0.16*dLat, lonlim(1)+0.17*dLon; latlim(1)+0.16*dLat, lonlim(2)-0.17*dLon; ...
            latlim(2)-0.14*dLat, lonlim(1)+0.17*dLon; latlim(2)-0.14*dLat, lonlim(2)-0.17*dLon];
    score = zeros(4,1);
    for k = 1:4
        score(k) = min(hypot(occ(:,1) - cand(k,1), (occ(:,2) - cand(k,2))*cosd(la0)));
    end
    [~, kBest] = max(score);
    Lw = 0.11*dLat;
    aLat = cand(kBest,1) - 0.5*Lw*MeanV/MeanWSPD;
    aLon = cand(kBest,2) - 0.5*Lw*MeanU/MeanWSPD/cosd(la0);
    lLat = cand(kBest,1) - 0.10*dLat;  lLon = cand(kBest,2);
end

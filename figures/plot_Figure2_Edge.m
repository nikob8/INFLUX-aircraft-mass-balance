%% plot_Figure2_Edge.m
%  Figure 2: the Downwind Edge approach for the first downwind transect of 10/04/2012.
%
%  Revision of 2026-11-06 (second pass): latitude and longitude labels fully visible and
%  in the same font as the rest; units and CH4 written without the TeX interpreter when
%  an installed font has the superscript and subscript characters (MATLAB on Linux
%  otherwise draws TeX text in a serif fallback font).
%  Revision of 2026-11-06: redesigned for readability. The computation is the
%  same as in the previous version (and as in MassBalanceMethodA.m): the edge and
%  downwind segments, the background line, and the green and red patches use the
%  code of Method A verbatim, so the city flux is the Table S1 value. What changed
%  is the presentation:
%     (a) map zoomed to the county and the transect; the band swept by Marion
%         County along the mean wind is shaded, so the integrated part of the
%         transect is visible directly; the background segments are outlined in
%         black; wind arrow placed in an empty corner
%     (b) background segments shaded gray, Downwind Edge background dashed,
%         softer fill colors, plume edges marked, legend, flux box
%     one sans-serif font everywhere (the old figure mixed serif and sans), no
%     overall title (date and time are in the panel title and the caption)
%
%  Needs on the path: INFLUX_<Date>, LegAndProfileList/<Date>.mat, InterX.m,
%  alongTrackKm.m, PlotINFLUX.m or the interstate shapefile (optional, roads),
%  Landfill/<Date>_Loop<k>.mat (optional), <Date>/Solution-MethodA.mat (optional),
%  Indy_Mass_Balance_Table_v5.csv (optional; data/ folder of the repository).
%  Requires the Mapping Toolbox (worldmap) like the other map figures.
%  Output: Fig2_Edge_<Date>_T<k>.png (300 dpi) and .pdf

clear; close all; clc;

%% ------------------------------------------------------------ settings
Date    = '20121004';
LegLoop = 1;                                  % first downwind transect
MB_ROOT = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';
OBS_DIR = fullfile(MB_ROOT,'AircraftObs');    % folder with INFLUX_<Date>  (change if elsewhere)
ROADS_CODE = '/storage/work/zrb5027/WRFGUI_V3/Code';                         % has PlotINFLUX.m (Zach)
ROADS_SHP  = '/storage/work/nvb5011/INFLUX_Paper_2/Thomas_2016/Shapefiles';  % Indianapolis_Interstate_Mileposts.shp
SOL_ROOT   = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB';   % <Date>/Solution-MethodA.mat
TABLE_S1   = 'Indy_Mass_Balance_Table_v5.csv';   % flux labels are taken from Table S1 when this file is found
if isfolder(MB_ROOT), addpath(MB_ROOT); end
if isfolder(OBS_DIR), addpath(OBS_DIR); end
if isfolder(ROADS_CODE), addpath(ROADS_CODE); end
if isfolder(ROADS_SHP),  addpath(ROADS_SHP);  end
outName = sprintf('Fig2_Edge_%s_T%d', Date, LegLoop);

SSLF_LATLON = [39.714691162109375, -86.205429077148438];   % South Side Landfill
PEAK_WIN_KM = 6;        % search for the landfill peak within this distance of the SSLF crossing
SMOOTH_PTS  = 15;       % extra smoothing used ONLY to find the plume edge minima
FLUX_CITY_TABLE = 50.23;   % Table S1 v5, 10/04/2012 transect 1 (used if the csv is not found)
FLUX_LF_TABLE   = 31.67;

% colours and font
C_BAND = [0.90 0.92 0.95];     % band swept by the county along the wind (map)
C_BGSEG = [0.91 0.91 0.91];    % background segments (panel b)
C_POS  = [0.55 0.78 0.47];     % positive enhancement
C_NEG  = [0.91 0.50 0.45];     % negative enhancement
C_LF   = [0.16 0.42 0.80];     % SSLF plume
C_ROAD = [0.62 0.62 0.62];
FS = 11;                        % base font size
[FONT, UNI] = pick_font();
set(0, 'DefaultAxesFontName', FONT, 'DefaultTextFontName', FONT);
if UNI   % plain text with Unicode sub- and superscripts
    INT = 'none';  sCH4 = ['CH' char(8324)];  sPerS = [' s' char(8315) char(185)];  numstr = @(v) strrep(sprintf('%.0f', v), '-', char(8722));
else     % TeX markup
    INT = 'tex';   sCH4 = 'CH_4';               sPerS = ' s^{-1}';                     numstr = @(v) sprintf('%.0f', v);
end
fprintf('Font used for all text: %s (%s)\n', FONT, ternary(UNI, 'Unicode sub/superscripts', 'TeX sub/superscripts'));

% ---- geometry, identical to MassBalanceMethodA.m --------------------------
%BottomRight BottomLeft TopLeft TopRight
CountyCorners=[39.63760423971629, -85.95148383144084; 39.63337385377353, -86.32570560518701; 39.92446141560437, -86.32734193809195; 39.9271330073512, -85.95582253907553];

%% ------------------------------------------------------------ load (as Method A)
load(fullfile(MB_ROOT,'LegAndProfileList',Date));          % LegList, ProfileList, UpwindLeg?, CH4Flagged?
if exist('UpwindLeg','var')==1
    UpwindLeg=UpwindLeg(1):UpwindLeg(2);
    AllLeg=[UpwindLeg];
else
    AllLeg=[];
end
for k=1:size(LegList,1)
    AllLeg=[AllLeg,LegList(k,1):LegList(k,2)]; %#ok<AGROW>
end

x=load(['INFLUX_',Date]);
x(x<-998)=NaN;
Time=x(:,1); Lat=x(:,2); Lon=x(:,3); Alt=x(:,4); %#ok<NASGU>
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

if median(CH4,'omitnan') < 100, CH4ppb = CH4*1000; else, CH4ppb = CH4; end   % plot in ppb

%% ------------------------------------------------------------ edge segmentation (verbatim Method A)
if (MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0) %#ok<*OR2,*AND2>
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
    if abs(Intersection1(1)-Intersection2(1)) < 0.05
        Upwind1=find(Lon(Leg)<Intersection1(2));
    elseif abs(Intersection1(2)-Intersection2(2)) < 0.05
        Upwind1=find(Lat(Leg)>Intersection1(1));
    else
        Upwind1=find(Lat(Leg)>Intersection1(1) | Lon(Leg)<Intersection1(2));
    end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05
        Upwind2=find(Lon(Leg)>Intersection2(2));
    elseif abs(Intersection2(2)-Intersection2(2)) < 0.05
        Upwind2=find(Lat(Leg)<Intersection2(1));
    else
        Upwind2=find(Lat(Leg)<Intersection2(1) | Lon(Leg)>Intersection2(2));
    end
elseif((MeanU>0 & MeanV<0) | (MeanU<0 & MeanV>0))
    if isempty(Intersection1)==1, Intersection1=[min(Lat(Leg)) min(Lon(Leg))]; end
    if isempty(Intersection2)==1, Intersection2=[max(Lat(Leg)) max(Lon(Leg))]; end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05
        Upwind1=find(Lon(Leg)<Intersection1(2));
    elseif abs(Intersection1(2)-Intersection2(2)) < 0.05
        Upwind1=find(Lat(Leg)<Intersection1(1));
    else
        Upwind1=find(Lat(Leg)<Intersection1(1) | Lon(Leg)<Intersection1(2));
    end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05
        Upwind2=find(Lon(Leg)>Intersection2(2));
    elseif abs(Intersection2(2)-Intersection2(2)) < 0.05
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

% ---- edge backgrounds (verbatim Method A; bgIdx1/bgIdx2 record which points were used)
Upwind1Flag = 1; if isempty(Upwind1), Upwind1Flag = -1; end
Upwind2Flag = 1; if isempty(Upwind2), Upwind2Flag = -1; end
Upwind1Background=median(CH4(Upwind1),'omitnan');
Upwind2Background=median(CH4(Upwind2),'omitnan');
bgIdx1 = Upwind1;  bgIdx2 = Upwind2;
mTemp=round(length(Leg)/10);
if Upwind1Flag==-1
    CheckSpot1=Leg(1)-Upwind2(1);  CheckSpot2=Leg(end)-Upwind2(end);
    if CheckSpot1==0,     Upwind1Background=median(CH4(Leg(end-mTemp:end)),'omitnan'); bgIdx1 = Leg(end-mTemp:end);
    elseif CheckSpot2==0, Upwind1Background=median(CH4(Leg(1:mTemp)),'omitnan');       bgIdx1 = Leg(1:mTemp);
    else, error('Check upwind edge methodology'); end
elseif Upwind2Flag==-1
    CheckSpot1=Leg(1)-Upwind1(1);  CheckSpot2=Leg(end)-Upwind1(end);
    if CheckSpot1==0,     Upwind2Background=median(CH4(Leg(end-mTemp:end)),'omitnan'); bgIdx2 = Leg(end-mTemp:end);
    elseif CheckSpot2==0, Upwind2Background=median(CH4(Leg(1:mTemp)),'omitnan');       bgIdx2 = Leg(1:mTemp);
    else, error('Check upwind edge methodology'); end
end

LatDistance=distance(min(Lat(Leg)),mean(Lon(Leg),'omitnan'),max(Lat(Leg)),mean(Lon(Leg),'omitnan'));
LonDistance=distance(mean(Lat(Leg),'omitnan'),min(Lon(Leg)),mean(Lat(Leg),'omitnan'),max(Lon(Leg)));
if LonDistance>LatDistance
    if Lon(Leg(1))>Lon(Leg(end)), [Upwind1Background,Upwind2Background]=deal(Upwind2Background,Upwind1Background); end
    LineCoeff=polyfit([Lon(Leg(1)),Lon(Leg(end))],[Upwind1Background,Upwind2Background],1);
    BackgroundLine=LineCoeff(1).*Lon(Leg)+LineCoeff(2);
    orient='east';
else
    if (MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0), [Upwind1Background,Upwind2Background]=deal(Upwind2Background,Upwind1Background); end
    if Lat(Leg(1))>Lat(Leg(end)), [Upwind1Background,Upwind2Background]=deal(Upwind2Background,Upwind1Background); end
    LineCoeff=polyfit([Lat(Leg(1)),Lat(Leg(end))],[Upwind1Background,Upwind2Background],1);
    BackgroundLine=LineCoeff(1).*Lat(Leg)+LineCoeff(2);
    orient='north';
end
if median(CH4,'omitnan') < 100, BackgroundLine = BackgroundLine*1000; end    % ppb
isDown = ismember(Leg,Downwind);
BackgroundLineLeg = BackgroundLine(isDown);

%% ------------------------------------------------------------ km axis
% origin = midpoint of the two county-corner crossings = where the mean-wind
% line through the centre of Marion County crosses the transect
latO = mean(Intersection(:,1));  lonO = mean(Intersection(:,2));
s = alongTrackKm(Lat(Leg), Lon(Leg), latO, lonO, orient);     % km, one per Leg point
s = s(:);
sDown = s(isDown);

%% ------------------------------------------------------------ flux labels
FluxCity = FLUX_CITY_TABLE;  FluxLFshow = FLUX_LF_TABLE;  src = 'constants in the settings block';
if isfile(TABLE_S1) && ~isempty(which('readtable'))
    T = readtable(TABLE_S1);
    if isdatetime(T.Date), dd = T.Date == datetime(Date, 'InputFormat', 'yyyyMMdd');
    else, dd = strcmp(cellstr(string(T.Date)), datestr(datenum(Date,'yyyymmdd'),'yyyy-mm-dd')); end
    r = find(dd & T.Transect == LegLoop, 1);
    if ~isempty(r)
        FluxCity = T.Flux_DownwindEdge_mols(r);  FluxLFshow = T.Flux_Landfill_mols(r);  src = TABLE_S1;
    end
end
solFile = fullfile(SOL_ROOT, Date, 'Solution-MethodA.mat');
if isfile(solFile)
    S = load(solFile, 'Flux');
    fprintf('Solution-MethodA.mat flux of this transect: %.1f mol/s\n', S.Flux(LegLoop));
end
fprintf('Flux labels from %s: Marion County %.1f mol/s, SSLF %.1f mol/s\n', src, FluxCity, FluxLFshow);

legTime = mean(Time(Leg),'omitnan');
hh = floor(legTime/3600);  mm = round((legTime - 3600*hh)/60);
if mm == 60, hh = hh + 1; mm = 0; end
fprintf('%s T%d: %02d:%02d UTC, wind %.0f deg / %.1f m/s, background points %d + %d, downwind %d\n', ...
    Date, LegLoop, hh, mm, MeanWDIR, MeanWSPD, numel(bgIdx1), numel(bgIdx2), numel(Downwind));

%% ------------------------------------------------------------ landfill plume (as before)
% (1) where does the SSLF project onto this transect along the mean wind?
PC = [SSLF_LATLON(1)+MeanV, SSLF_LATLON(2)+MeanU*(111/86)];
PA = [SSLF_LATLON(1)-MeanV, SSLF_LATLON(2)-MeanU*(111/86)];
LineLF = [PA(2) SSLF_LATLON(2) PC(2); PA(1) SSLF_LATLON(1) PC(1)];
XLF = flipud(InterX(LineLF, Line2))';                  % [lat lon]
lfI = [];  lfSource = 'none';

% (2) edge indices saved by the landfill script, if the Loop file has them
lfFile = fullfile(MB_ROOT, 'Landfill', sprintf('%s_Loop%d.mat', Date, LegLoop));
if isfile(lfFile)
    L = load(lfFile);  fn = fieldnames(L);
    for q = 1:numel(fn)                                 % a 2-element index pair inside this leg?
        v = double(L.(fn{q}));
        if isnumeric(v) && numel(v) == 2 && all(v == round(v))
            if all(v >= Leg(1) & v <= Leg(end))          % absolute indices
                lfI = sort(v(:)' - Leg(1) + 1);  lfSource = ['file: ' fn{q}]; break
            elseif all(v >= 1 & v <= numel(Leg)) && diff(sort(v)) > 3   % indices within the leg
                lfI = sort(v(:)');               lfSource = ['file: ' fn{q}]; break
            end
        end
    end
end

% (3) otherwise detect the plume around the SSLF crossing
if isempty(lfI) && ~isempty(XLF)
    [~, eAx] = alongTrackKm(Lat(Leg), Lon(Leg), latO, lonO, orient);
    la0 = mean(Lat(Leg),'omitnan');
    sX  = (XLF(1,2)-lonO)*111.32*cosd(la0)*eAx(1) + (XLF(1,1)-latO)*110.57*eAx(2);
    enh = CH4ppb(Leg) - BackgroundLine;
    win = find(abs(s - sX) <= PEAK_WIN_KM & isDown(:));
    if ~isempty(win)
        [~, k] = max(enh(win));  iPk = win(k);
        es = movmean(CH4ppb(Leg), SMOOTH_PTS, 'omitnan') - BackgroundLine;
        while iPk > 1          && es(iPk-1) > es(iPk), iPk = iPk - 1; end   % to the smoothed peak
        while iPk < numel(Leg) && es(iPk+1) > es(iPk), iPk = iPk + 1; end
        iL = iPk; while iL > 1          && es(iL-1) <= es(iL), iL = iL - 1; end
        iR = iPk; while iR < numel(Leg) && es(iR+1) <= es(iR), iR = iR + 1; end
        if iR - iL >= 3, lfI = [iL iR]; lfSource = 'detected around the SSLF crossing'; end
    end
end
if ~isempty(lfI)
    fprintf('Landfill plume (%s): s = %.1f to %.1f km\n', lfSource, s(lfI(1)), s(lfI(2)));
else
    fprintf('No landfill plume identified on this transect.\n');
end

%% ------------------------------------------------------------ figure
fig = figure('Color','w','Units','pixels','Position',[80 80 1250 520]);
yTr = CH4ppb(Leg);  yTr = yTr(:);
lim = [min(yTr) max(yTr)];

% ---- (a) map --------------------------------------------------------------
ax1 = axes('Parent', fig, 'Position', [0.085 0.11 0.33 0.79]);
allLat = [CountyCorners(:,1); Lat(Leg)];  allLon = [CountyCorners(:,2); Lon(Leg)];
latlim = [min(allLat) max(allLat)] + [-0.07 0.07];
lonlim = [min(allLon) max(allLon)] + [-0.09 0.09];
ax1 = worldmap(latlim, lonlim);
set(ax1, 'SortMethod', 'childorder');
setm(ax1, 'FontName', FONT, 'FontSize', FS-1, 'FLineWidth', 1, 'GColor', [0.82 0.82 0.82], ...
    'GLineStyle', '-', 'GLineWidth', 0.4, 'MLabelLocation', 0.25, 'PLabelLocation', 0.2, ...
    'MLabelRound', -2, 'PLabelRound', -1);
hold on

% band swept by the county along the mean wind, from the county to the transect
la0 = mean(Lat(Leg), 'omitnan');
dwScore = (CountyCorners(:,2) - mean(CountyCorners(:,2)))*cosd(la0)*MeanU + (CountyCorners(:,1) - mean(CountyCorners(:,1)))*MeanV;
[~, iDown] = max(dwScore);                       % most downwind corner of the county
bandLat = [PointB1(1); Intersection1(1,1); Intersection2(1,1); PointB2(1); CountyCorners(iDown,1)];
bandLon = [PointB1(2); Intersection1(1,2); Intersection2(1,2); PointB2(2); CountyCorners(iDown,2)];
hBand = patchm(bandLat, bandLon, C_BAND);
set(hBand, 'EdgeColor', 'none');

draw_roads(ax1, C_ROAD, hBand);                                   % Indianapolis highways
cc = [CountyCorners; CountyCorners(1,:)];
plotm(cc(:,1), cc(:,2), '-', 'Color', [0 0 0], 'LineWidth', 1.3); % study area (Marion County)

% transect: background segments outlined in black, then the points coloured by CH4
for seg = {bgIdx1, bgIdx2}
    ii = seg{1};
    if ~isempty(ii), plotm(Lat(ii), Lon(ii), '-', 'Color', [0 0 0], 'LineWidth', 8); end
end
scatterm(Lat(Leg), Lon(Leg), 22, yTr, 'filled');
caxis(ax1, lim);  colormap(ax1, parula(256));
cb = colorbar(ax1, 'eastoutside');
set(cb, 'FontName', FONT, 'FontSize', FS-1);
try
    set(get(cb, 'Label'), 'String', [sCH4 ' (ppb)'], 'Interpreter', INT, 'FontName', FONT, 'FontSize', FS);
catch
    ylabel(cb, [sCH4 ' (ppb)'], 'Interpreter', INT, 'FontName', FONT, 'FontSize', FS);
end

% SSLF
plotm(SSLF_LATLON(1), SSLF_LATLON(2), '^', 'MarkerSize', 9, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [1 1 1], 'LineWidth', 1);
textm(SSLF_LATLON(1) - 0.035, SSLF_LATLON(2) + 0.015, 'SSLF', 'FontName', FONT, 'FontSize', FS-1);

% mean wind: arrow in the corner of the map farthest from the transect and the county
dLat = diff(latlim);  dLon = diff(lonlim);
cand = [latlim(1)+0.16*dLat, lonlim(1)+0.17*dLon; latlim(1)+0.16*dLat, lonlim(2)-0.17*dLon; ...
        latlim(2)-0.14*dLat, lonlim(1)+0.17*dLon; latlim(2)-0.14*dLat, lonlim(2)-0.17*dLon];
occ = [Lat(Leg) Lon(Leg); CountyCorners; SSLF_LATLON];
score = zeros(4,1);
for k = 1:4
    score(k) = min(hypot(occ(:,1) - cand(k,1), (occ(:,2) - cand(k,2))*cosd(la0)));
end
[~, kBest] = max(score);
Lw = 0.11*dLat;                                   % arrow length (degrees of latitude)
aLat = cand(kBest,1) - 0.5*Lw*MeanV/MeanWSPD;
aLon = cand(kBest,2) - 0.5*Lw*MeanU/MeanWSPD/cosd(la0);
hq = quiverm(aLat, aLon, Lw*MeanV/MeanWSPD, Lw*MeanU/MeanWSPD/cosd(la0), 'k', 0);
set(hq, 'LineWidth', 1.6);
textm(cand(kBest,1) - 0.10*dLat, cand(kBest,2), [sprintf('%.1f m', MeanWSPD) sPerS], ...
    'Interpreter', INT, 'FontName', FONT, 'FontSize', FS-1, 'HorizontalAlignment', 'center');
text(ax1, 0, 1.06, sprintf('(a) Transect %d, %s/%s/%s, %02d:%02d UTC', LegLoop, Date(5:6), Date(7:8), Date(1:4), hh, mm), 'Units', 'normalized', 'FontSize', FS+1, 'FontName', FONT);

% ---- (b) CH4 along the transect -------------------------------------------
ax2 = axes('Parent', fig, 'Position', [0.53 0.14 0.45 0.76]);
hold(ax2, 'on');  box(ax2, 'on');
yl = [floor(lim(1)) - 3, ceil(lim(2) + 0.40*diff(lim))];   % head room for the legend
xl = [min(s) max(s)];
set(ax2, 'XLim', xl, 'YLim', yl);
% background segments
hBg = gobjects(0);
for seg = {bgIdx1, bgIdx2}
    ii = seg{1};
    if isempty(ii), continue; end
    sb = s(ismember(Leg, ii));
    hBg(end+1) = patch(ax2, [min(sb) max(sb) max(sb) min(sb)], [yl(1) yl(1) yl(2) yl(2)], C_BGSEG, 'EdgeColor', 'none'); %#ok<SAGROW>
end
% positive and negative enhancement over the downwind segment (as Method A)
xD = sDown(:);  yD = CH4ppb(Downwind);  yD = yD(:);  bD = BackgroundLineLeg(:);
y1 = max(yD, bD);  g = ~isnan(y1);
hPos = patch(ax2, [xD(g); flipud(xD(g))], [bD(g); flipud(y1(g))], C_POS, 'EdgeColor', 'none');
y1 = min(yD, bD);  g = ~isnan(y1);
hNeg = patch(ax2, [xD(g); flipud(xD(g))], [bD(g); flipud(y1(g))], C_NEG, 'EdgeColor', 'none');
% SSLF plume above the straight line joining its two edges
hLF = gobjects(0);
if ~isempty(lfI)
    ii = (lfI(1):lfI(2))';
    lfLine = interp1([s(ii(1)) s(ii(end))], [yTr(ii(1)) yTr(ii(end))], s(ii));
    yL = max(yTr(ii), lfLine);
    hLF = patch(ax2, [s(ii); flipud(s(ii))], [lfLine; flipud(yL)], C_LF, 'EdgeColor', 'none');
    plot(ax2, s(ii([1 end])), yTr(ii([1 end])), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 1 1], 'LineWidth', 1);
    [~, ip] = max(yL - lfLine);
    xPk = s(ii(ip));  yPk = yTr(ii(ip));
    xT = xPk + 0.07*diff(xl);  yT = yPk + 0.02*diff(yl);
    plot(ax2, [xPk + 0.008*diff(xl), xT - 0.005*diff(xl)], [yPk, yT], 'k-', 'LineWidth', 0.8);
    text(ax2, xT, yT, 'South Side Landfill', 'FontSize', FS-1, 'FontName', FONT, 'VerticalAlignment', 'middle');
end
hLine = plot(ax2, s, BackgroundLine(:), '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.4);
hTr = plot(ax2, s, yTr, '-', 'Color', [0 0 0], 'LineWidth', 1.3);
set(ax2, 'FontSize', FS, 'FontName', FONT, 'TickDir', 'out', 'Layer', 'top', 'XLim', xl, 'YLim', yl);
xlabel(ax2, 'Distance along flight track (km)', 'FontSize', FS+1);
ylabel(ax2, [sCH4 ' mole fraction (ppb)'], 'Interpreter', INT, 'FontSize', FS+1);
text(ax2, 0, 1.06, '(b) Downwind Edge background and enhancement', 'Units', 'normalized', 'FontSize', FS+1, 'FontName', FONT);

hs = [hTr, hLine, hPos, hNeg];  labs = {sCH4, 'Downwind Edge background', 'Positive enhancement', 'Negative enhancement'};
if ~isempty(hLF), hs = [hs, hLF]; labs{end+1} = 'SSLF plume'; end
if ~isempty(hBg), hs = [hs, hBg(1)]; labs{end+1} = 'Background segments'; end
lg = legend(ax2, hs, labs, 'Location', 'northwest', 'FontSize', FS-1, 'Interpreter', INT);
set(lg, 'FontName', FONT, 'Box', 'on');
try, set(lg, 'NumColumns', 2, 'ItemTokenSize', [20 12]); catch, end   % two columns of three entries

% emission rates of this transect
text(ax2, 0.025, 0.11, ['Marion County: ' numstr(FluxCity) ' mol' sPerS], 'Units', 'normalized', ...
    'Interpreter', INT, 'FontSize', FS, 'FontName', FONT, 'Color', [0 0 0]);
if ~isempty(lfI) && ~isnan(FluxLFshow)
    text(ax2, 0.025, 0.045, ['SSLF: ' numstr(FluxLFshow) ' mol' sPerS], 'Units', 'normalized', ...
        'Interpreter', INT, 'FontSize', FS, 'FontName', FONT, 'Color', C_LF*0.85);
end

%% ------------------------------------------------------------ save
fix_map_labels(ax1, FONT, FS-1);          % degree signs without TeX, same font as the rest
if ~isempty(which('exportgraphics'))
    exportgraphics(fig, [outName '.png'], 'Resolution', 300);
    exportgraphics(fig, [outName '.pdf'], 'ContentType', 'vector');
    fprintf('Saved %s.png / .pdf\n', outName);
else
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, [outName '.png'], '-dpng', '-r300');
    fprintf('Saved %s.png\n', outName);
end


%% ============================================================ local functions
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

%% plot_Figure3_5_Upwind.m
%  Revised Figure 3 (Upwind approach: map + downwind transect) and revised
%  Figure 5 (Upwind, Adjusted ABL Upwind and Upwind Edge backgrounds side by
%  side) for the first downwind transect of 10/04/2012.
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
%  Presentation changes (Reviewer 1, R1.6 / R1.7): SSLF marked on the map,
%  x-axis in km along the flight track with the same origin and sign convention
%  as Figure 2 (alongTrackKm.m), upwind transect outlined in black, identical
%  axis limits on the three panels of Figure 5.
%
%  Needs on the path: INFLUX_<Date>, LegAndProfileList/<Date>.mat, InterX.m,
%  alongTrackKm.m, PlotINFLUX.m or the Indianapolis interstate shapefile,
%  <Date>/Solution-MethodB.mat, -MethodC.mat, -MethodD.mat.

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
SOL_ROOT = pwd;         % folder holding <Date>/Solution-MethodB/C/D.mat from YOUR rerun (cos theta); MB_ROOT has Zach's originals

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
for UpwindLoop=1:length(Leg)
    PointB=[Lat(Leg(UpwindLoop)), Lon(Leg(UpwindLoop))];
    PointC=[PointB(1)+MeanV, PointB(2)+MeanU*(111/86)];
    PointA=[PointB(1)-MeanV, PointB(2)-MeanU*(111/86)];
    Line1A=[PointA(2) PointB(2) PointC(2);PointA(1) PointB(1) PointC(1)];
    UpwindIntersection=flipud(InterX(Line1A,UpwindLine'))';
    if isempty(UpwindIntersection)==1
        BackgroundLine(UpwindLoop)=NaN; %#ok<SAGROW>
    else
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
fprintf('Flux labels: Upwind %.0f, Adjusted ABL Upwind %.0f, Upwind Edge %.0f mol/s  (Table S1 for this transect: -41, 4, 34)\n', FluxB, FluxC, FluxD);

%% ------------------------------------------------------------ Figure 3
lim = [min([CH4ppb(Leg); CH4ppb(UpwindLeg)]) max([CH4ppb(Leg); CH4ppb(UpwindLeg)])];
fig3 = figure('Color','w','Units','pixels','Position',[80 80 1250 560]);
sgtitle(sprintf('Upwind Approach: %s/%s/%s', Date(5:6), Date(7:8), Date(1:4)), 'FontSize', 15, 'FontWeight', 'bold');

% ---- (a) map --------------------------------------------------------------
ax1 = subplot(1,2,1);
worldmap([39.2 40.6],[-87 -85.2]); setm(gca,'FontSize',12); hold on
scatterm(Lat(Upwind1),Lon(Upwind1),200,[0 0 0],'filled')        % black outline: downwind transect ends
scatterm(Lat(Downwind),Lon(Downwind),200,[1 0 0],'filled')      % red outline: directly downwind of Marion County
scatterm(Lat(Upwind2),Lon(Upwind2),200,[0 0 0],'filled')
scatterm(Lat(UpwindLeg),Lon(UpwindLeg),200,[0 0 0],'filled')    % black outline: upwind transect
scatterm(Lat(UpwindLeg),Lon(UpwindLeg),50,CH4ppb(UpwindLeg),'filled')
scatterm(Lat(Leg),Lon(Leg),50,CH4ppb(Leg),'filled')
caxis(lim); colormap(ax1, parula);
c = colorbar; c.FontSize = 12; c.Label.String = 'Methane Mole Fraction (ppb)'; c.Label.FontSize = 12;
draw_roads();
for k = 1:4
    k2 = mod(k,4)+1;
    plotm([CountyCorners(k,1),CountyCorners(k2,1)],[CountyCorners(k,2),CountyCorners(k2,2)],'LineWidth',3,'Color',[0 0 0])
end
plotm(SSLF_LATLON(1), SSLF_LATLON(2), '^', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'LineWidth', 1);
aLat = 39.42; aLon = -85.55; L = 0.25;
quiverm(aLat, aLon, L*MeanV/MeanWSPD, L*MeanU/MeanWSPD*(111/86), 'k', 0);
textm(39.30, -85.75, sprintf('U = %.1f m s^{-1}', MeanWSPD), 'FontSize', 12, 'FontWeight', 'bold');
text(0.03, 0.95, 'a.)', 'Units','normalized', 'FontSize', 13, 'FontWeight','bold');

% ---- (b) downwind transect with the Upwind background ----------------------
ax2 = subplot(1,2,2);
draw_panel(ax2, s, CH4ppb(Leg), bgB, isDown, 'b.)', FluxB, '');
ylabel('Methane Mole Fraction (ppb)', 'FontSize', 13)

set(fig3, 'PaperPositionMode','auto');
print(fig3, [outName3 '.png'], '-dpng', '-r300');
print(fig3, [outName3 '.eps'], '-depsc');
fprintf('Saved %s.png / .eps\n', outName3);

%% ------------------------------------------------------------ Figure 5
fig5 = figure('Color','w','Units','pixels','Position',[80 80 1500 470]);
yAll = [CH4ppb(Leg); bgB; bgC; bgD]; yl = [min(yAll)-2 max(yAll)+2];
ax = gobjects(1,3);
ax(1) = subplot(1,3,1); draw_panel(ax(1), s, CH4ppb(Leg), bgB, isDown, 'a.)', FluxB, 'Upwind');
ylabel('Methane Mole Fraction (ppb)', 'FontSize', 13)
ax(2) = subplot(1,3,2); draw_panel(ax(2), s, CH4ppb(Leg), bgC, isDown, 'b.)', FluxC, 'Adjusted ABL Upwind');
ax(3) = subplot(1,3,3); draw_panel(ax(3), s, CH4ppb(Leg), bgD, isDown, 'c.)', FluxD, 'Upwind Edge');
for k = 1:3, ylim(ax(k), yl); end                                % identical limits on all panels
set(fig5, 'PaperPositionMode','auto');
print(fig5, [outName5 '.png'], '-dpng', '-r300');
print(fig5, [outName5 '.eps'], '-depsc');
fprintf('Saved %s.png / .eps\n', outName5);


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

function draw_panel(ax, s, trace, bg, isDown, letter, fluxVal, ttl)
% downwind trace (black), background (blue), green/red enhancement patches over
% the part of the transect directly downwind of Marion County, flux label
    axes(ax); hold(ax,'on'); box(ax,'on');
    plot(ax, s, trace, 'Color',[0 0 0], 'LineWidth', 1.5)
    plot(ax, s, bg,    'Color', 'blue', 'LineWidth', 1.5)
    xD = s(isDown); yD = trace(isDown); bD = bg(isDown);
    y1 = yD; bad = (y1 - bD) < 0; y1(bad) = bD(bad); g = ~isnan(y1);
    patch(ax, [xD(g)' flip(xD(g))'], [bD(g)' flip(y1(g)')], 'g', 'EdgeColor','none')
    y1 = yD; bad = (y1 - bD) > 0; y1(bad) = bD(bad); g = ~isnan(y1);
    patch(ax, [xD(g)' flip(xD(g))'], [bD(g)' flip(y1(g)')], 'r', 'EdgeColor','none')
    plot(ax, s, trace, 'Color',[0 0 0], 'LineWidth', 1.5)          % trace on top
    plot(ax, s, bg,    'Color', 'blue', 'LineWidth', 1.5)
    xlabel(ax, 'Distance along flight track (km)', 'FontSize', 13)
    xlim(ax, [min(s) max(s)]); set(ax, 'FontSize', 12)
    text(ax, 0.03, 0.95, letter, 'Units','normalized', 'FontSize', 13, 'FontWeight','bold');
    text(ax, 0.05, 0.08, sprintf('Flux = %.0f mol s^{-1}', fluxVal), 'Units','normalized', 'FontSize', 12, 'FontWeight','bold');
    if ~isempty(ttl), title(ax, ttl, 'FontSize', 13, 'FontWeight','bold'); end
end

function draw_roads()
% Major Indianapolis highways.  Uses Zach's PlotINFLUX.m (the routine behind the
% submitted figures) if it is on the path; otherwise draws the interstates from
% Indianapolis_Interstate_Mileposts.shp (Indiana State Plane East, US survey feet),
% as in draw_basemap() of Compare_WRF_Obs_v5.m.
    if ~isempty(which('PlotINFLUX'))
        PlotINFLUX;
        return
    end
    if isempty(which('Indianapolis_Interstate_Mileposts.shp'))
        warning('Neither PlotINFLUX.m nor Indianapolis_Interstate_Mileposts.shp is on the path -- no roads drawn.');
        return
    end
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
    scatterm(LATr, LONr, 6, [0.25 0.25 0.25], 'filled');
    warning(wsv);
end

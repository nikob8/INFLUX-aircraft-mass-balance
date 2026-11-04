%% plot_Figure2_Edge.m
%  Revised Figure 2: Downwind Edge approach example, first downwind transect of 10/04/2012
%  (Reviewer 1: mark the SSLF; x-axis in km along the flight track from the
%   centre of the urban plume.  The red/black outlines are explained in the caption.)
%
%  The edge/downwind segmentation, the background line and the green/red patches
%  are computed EXACTLY as in MassBalanceMethodA.m (same code), so the city flux is
%  the Table S1 value.  Changes are presentation only:
%     * panel (b) x-axis = distance along the flight track (km); 0 where the mean
%       wind line through the centre of Marion County crosses the transect;
%       increasing eastward (alongTrackKm.m)
%     * SSLF marked on the map (triangle)
%     * the blue landfill polygon is now COMPUTED instead of hand drawn:
%         1. the SSLF is projected along the mean wind onto the transect (the
%            same InterX construction the code uses for the county corners);
%         2. the plume peak is the largest enhancement within PEAK_WIN_KM of
%            that crossing point;
%         3. each plume edge is the first local minimum of the (lightly
%            smoothed) enhancement above the Edge background when walking
%            outward from the peak -- i.e. the
%            "single point measurements on either side of the plume" of
%            Section 2.4;
%         4. the polygon is the trace above the straight line joining the two
%            edge points.
%       If Landfill/<Date>_Loop<k>.mat holds edge indices, those are used
%       instead (the script prints what it finds in that file).
%
%  Needs on the path: INFLUX_<Date>, LegAndProfileList/<Date>.mat, InterX.m,
%  PlotINFLUX.m (optional), alongTrackKm.m, <Date>/Solution-MethodA.mat (optional).

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
outName = sprintf('Fig2_Edge_%s_T%d', Date, LegLoop);

SSLF_LATLON   = [39.714691162109375, -86.205429077148438];   % South Side Landfill
PEAK_WIN_KM   = 6;      % search for the landfill peak within this distance of the SSLF crossing
SMOOTH_PTS    = 15;     % extra smoothing used ONLY to find the edge minima
FLUX_LF_TABLE = 32;     % Table S1 v3 (with cos theta), 10/04/2012 transect 1; NaN = use computed
SOL_ROOT = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB';   % <Date>/Solution-MethodA.mat from YOUR rerun (cos theta)

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

if median(CH4,'omitnan') < 100, CH4ppb = CH4*1000; else, CH4ppb = CH4; end   % plot in ppb

%% ------------------------------------------------------------ edge segmentation (verbatim Method A)
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

% ---- edge backgrounds (verbatim Method A) ----------------------------------
Upwind1Flag = 1; if isempty(Upwind1), Upwind1Flag = -1; end
Upwind2Flag = 1; if isempty(Upwind2), Upwind2Flag = -1; end
Upwind1Background=median(CH4(Upwind1),'omitnan');
Upwind2Background=median(CH4(Upwind2),'omitnan');
mTemp=round(length(Leg)/10);
if Upwind1Flag==-1
    CheckSpot1=Leg(1)-Upwind2(1);  CheckSpot2=Leg(end)-Upwind2(end);
    if CheckSpot1==0,     Upwind1Background=median(CH4(Leg(end-mTemp:end)),'omitnan');
    elseif CheckSpot2==0, Upwind1Background=median(CH4(Leg(1:mTemp)),'omitnan');
    else, error('Check upwind edge methodology'); end
elseif Upwind2Flag==-1
    CheckSpot1=Leg(1)-Upwind1(1);  CheckSpot2=Leg(end)-Upwind1(end);
    if CheckSpot1==0,     Upwind2Background=median(CH4(Leg(end-mTemp:end)),'omitnan');
    elseif CheckSpot2==0, Upwind2Background=median(CH4(Leg(1:mTemp)),'omitnan');
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

%% ------------------------------------------------------------ km axis (new)
% origin = midpoint of the two county-corner crossings = where the mean-wind
% line through the centre of Marion County crosses the transect
latO = mean(Intersection(:,1));  lonO = mean(Intersection(:,2));
s = alongTrackKm(Lat(Leg), Lon(Leg), latO, lonO, orient);     % km, one per Leg point
sDown = s(isDown);

% flux for the annotation (from the saved solution, if present)
FluxCity = NaN;
solFile = fullfile(SOL_ROOT, Date, 'Solution-MethodA.mat');
assert(isfile(solFile), 'Solution-MethodA.mat with cos(theta) not found at %s (fix SOL_ROOT)', solFile);
S = load(solFile, 'Flux'); FluxCity = S.Flux(LegLoop);
fprintf('Flux label from %s: %.1f mol/s (Table S1 v4: 50)\n', solFile, FluxCity);

legTime = mean(Time(Leg),'omitnan');
fprintf('%s T%d: %02d:%02d UTC, wind %.0f deg / %.1f m/s, edge points %d + %d, downwind %d, flux %.0f mol/s\n', ...
    Date, LegLoop, floor(legTime/3600), round(mod(legTime,3600)/60), MeanWDIR, MeanWSPD, numel(Upwind1), numel(Upwind2), numel(Downwind), FluxCity);

%% ------------------------------------------------------------ landfill plume
% (1) where does the SSLF project onto this transect along the mean wind?
PC = [SSLF_LATLON(1)+MeanV, SSLF_LATLON(2)+MeanU*(111/86)];
PA = [SSLF_LATLON(1)-MeanV, SSLF_LATLON(2)-MeanU*(111/86)];
LineLF = [PA(2) SSLF_LATLON(2) PC(2); PA(1) SSLF_LATLON(1) PC(1)];
XLF = flipud(InterX(LineLF, Line2))';                  % [lat lon]
lfI = [];  FluxLF = NaN;  lfSource = 'none';

% (2) edge indices saved by FindTheLandfill, if the Loop file has them
lfFile = fullfile(MB_ROOT, 'Landfill', sprintf('%s_Loop%d.mat', Date, LegLoop));
if isfile(lfFile)
    L = load(lfFile);  fn = fieldnames(L);
    fprintf('\nContents of %s:\n', lfFile);
    for q = 1:numel(fn)
        v = L.(fn{q});
        if isnumeric(v) && numel(v) <= 6
            fprintf('   %-24s %-10s %s\n', fn{q}, mat2str(size(v)), mat2str(v(:)', 6));
        else
            fprintf('   %-24s %-10s %s\n', fn{q}, mat2str(size(v)), class(v));
        end
    end
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
    for q = 1:numel(fn)                                  % a scalar flux?
        if ~isempty(regexpi(fn{q}, 'flux', 'once')) && isnumeric(L.(fn{q})) && isscalar(L.(fn{q}))
            FluxLF = double(L.(fn{q})); fprintf('   -> landfill flux in file: %s = %.1f mol/s\n', fn{q}, FluxLF);
        end
    end
end

% (3) otherwise detect the plume around the SSLF crossing
if isempty(lfI) && ~isempty(XLF)
    % position of the crossing on the same along-track axis as s
    [~, eAx] = alongTrackKm(Lat(Leg), Lon(Leg), latO, lonO, orient);
    la0 = mean(Lat(Leg),'omitnan');
    sX  = (XLF(1,2)-lonO)*111.32*cosd(la0)*eAx(1) + (XLF(1,1)-latO)*110.57*eAx(2);
    enh = CH4ppb(Leg) - BackgroundLine;
    win = find(abs(s - sX) <= PEAK_WIN_KM & isDown(:));
    if ~isempty(win)
        [~, k] = max(enh(win));  iPk = win(k);
        % walk on the smoothed enhancement ABOVE the Edge background, not on the raw
        % trace: with a sloping background the raw trace can keep falling past the
        % plume edge
        es = movmean(CH4ppb(Leg), SMOOTH_PTS, 'omitnan') - BackgroundLine;
        while iPk > 1          && es(iPk-1) > es(iPk), iPk = iPk - 1; end   % to the smoothed peak
        while iPk < numel(Leg) && es(iPk+1) > es(iPk), iPk = iPk + 1; end
        iL = iPk; while iL > 1          && es(iL-1) <= es(iL), iL = iL - 1; end
        iR = iPk; while iR < numel(Leg) && es(iR+1) <= es(iR), iR = iR + 1; end
        if iR - iL >= 3, lfI = [iL iR]; lfSource = 'detected around the SSLF crossing'; end
    end
end

if ~isempty(lfI)
    ii = lfI(1):lfI(2);
    lfLine = interp1([s(ii(1)) s(ii(end))], CH4ppb(Leg([ii(1) ii(end)])), s(ii));
    lfEnh  = max(CH4ppb(Leg(ii)) - lfLine, 0);                       % ppb
    widthKm = abs(s(ii(end)) - s(ii(1)));
    if isnan(FluxLF) && isfile(solFile)
        S2 = load(solFile, 'Density', 'DownwindABL');
        if isfield(S2,'Density') && isfield(S2,'DownwindABL')
            % same form as Method A: Density * WSPD * (enhancement integrated across the plume) * ABL
            FluxLF = S2.Density(LegLoop) * MeanWSPD * (abs(trapz(s(ii), lfEnh))*1000/1e9) * S2.DownwindABL(LegLoop);
        end
    end
    fprintf('\nLandfill plume (%s): s = %.1f to %.1f km (width %.1f km), peak +%.0f ppb, flux %.0f mol/s (Table S1: %g)\n', ...
        lfSource, s(ii(1)), s(ii(end)), widthKm, max(lfEnh), FluxLF, FLUX_LF_TABLE);
else
    fprintf('\nNo landfill plume identified on this transect.\n');
end

%% ------------------------------------------------------------ figure
fig = figure('Color','w','Units','pixels','Position',[80 80 1250 560]);
sgtitle(sprintf('Downwind Edge Approach: %s/%s/%s', Date(5:6), Date(7:8), Date(1:4)), 'FontSize', 15, 'FontWeight', 'bold');
lim = [min(CH4ppb(Leg)) max(CH4ppb(Leg))];

% ---- (a) map --------------------------------------------------------------
ax1 = subplot(1,2,1);
worldmap([39.2 40.6],[-87 -85.2]); setm(gca,'FontSize',12); hold on
scatterm(Lat(Upwind1),Lon(Upwind1),200,[0 0 0],'filled')        % black outline: edges
scatterm(Lat(Downwind),Lon(Downwind),200,[1 0 0],'filled')      % red outline: directly downwind of Marion County
scatterm(Lat(Upwind2),Lon(Upwind2),200,[0 0 0],'filled')
scatterm(Lat(Leg),Lon(Leg),50,CH4ppb(Leg),'filled')
caxis(lim); colormap(ax1, parula);
c = colorbar; c.FontSize = 12; c.Label.String = 'Methane Mole Fraction (ppb)'; c.Label.FontSize = 12;
draw_roads();                                                    % Indianapolis highways
for k = 1:4
    k2 = mod(k,4)+1;
    plotm([CountyCorners(k,1),CountyCorners(k2,1)],[CountyCorners(k,2),CountyCorners(k2,2)],'LineWidth',3,'Color',[0 0 0])
end
plotm(SSLF_LATLON(1), SSLF_LATLON(2), '^', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'LineWidth', 1);
% mean wind arrow (bottom right), length ~0.25 deg
aLat = 39.42; aLon = -85.55; L = 0.25;
quiverm(aLat, aLon, L*MeanV/MeanWSPD, L*MeanU/MeanWSPD*(111/86), 'k', 0);
textm(39.30, -85.75, sprintf('U = %.1f m s^{-1}', MeanWSPD), 'FontSize', 12, 'FontWeight', 'bold');
text(0.03, 0.95, 'a.)', 'Units','normalized', 'FontSize', 13, 'FontWeight','bold');

% ---- (b) along-track plot ---------------------------------------------------
ax2 = subplot(1,2,2); hold on; box on
plot(s, CH4ppb(Leg), 'Color',[0 0 0], 'LineWidth', 1.5)
plot(s, BackgroundLine, 'Color', 'blue', 'LineWidth', 1.5)
xD = sDown;  yD = CH4ppb(Downwind);
% green: above background
y1 = yD; bad = (y1 - BackgroundLineLeg) < 0; y1(bad) = BackgroundLineLeg(bad);
g = ~isnan(y1);
patch([xD(g)' flip(xD(g))'], [BackgroundLineLeg(g)' flip(y1(g)')], 'g', 'EdgeColor','none')
% red: below background
y1 = yD; bad = (y1 - BackgroundLineLeg) > 0; y1(bad) = BackgroundLineLeg(bad);
g = ~isnan(y1);
patch([xD(g)' flip(xD(g))'], [BackgroundLineLeg(g)' flip(y1(g)')], 'r', 'EdgeColor','none')
% blue: SSLF plume above the straight line between its two edges
if ~isempty(lfI)
    ii = lfI(1):lfI(2);
    lfLine = interp1([s(ii(1)) s(ii(end))], CH4ppb(Leg([ii(1) ii(end)])), s(ii));
    yL = CH4ppb(Leg(ii)); yL = max(yL, lfLine);
    patch([s(ii)' flip(s(ii))'], [lfLine' flip(yL')], [0.05 0.33 0.80], 'EdgeColor','none')
    [pk, ip] = max(yL - lfLine);
    xPk = s(ii(ip));  yPk = CH4ppb(Leg(ii(ip)));
    yr = ylim;  xr = xlim;
    xT = xPk + 0.12*diff(xr);  yT = yPk - 0.10*diff(yr);
    plot([xT xPk + 0.01*diff(xr)], [yT yPk - 0.02*diff(yr)], 'k-', 'LineWidth', 1);
    text(xT, yT, ' South Side Landfill', 'FontSize', 11, 'VerticalAlignment','top');
end
plot(s, CH4ppb(Leg), 'Color',[0 0 0], 'LineWidth', 1.5)          % redraw the trace on top
xlabel('Distance along flight track (km)', 'FontSize', 13)
ylabel('Methane Mole Fraction (ppb)', 'FontSize', 13)
xlim([min(s) max(s)]); set(gca, 'FontSize', 12)
text(0.03, 0.95, 'b.)', 'Units','normalized', 'FontSize', 13, 'FontWeight','bold');
txt = sprintf('Flux = %.0f mol s^{-1}', FluxCity);
text(0.05, 0.12, txt, 'Units','normalized', 'FontSize', 12, 'FontWeight','bold');
FluxLFshow = FLUX_LF_TABLE; if isnan(FluxLFshow), FluxLFshow = FluxLF; end
if ~isempty(lfI) && ~isnan(FluxLFshow)
    text(0.05, 0.05, sprintf('Flux_{LF} = %.0f mol s^{-1}', FluxLFshow), 'Units','normalized', 'FontSize', 12, 'FontWeight','bold', 'Color',[0.05 0.33 0.80]);
end

set(fig, 'PaperPositionMode','auto');
print(fig, [outName '.png'], '-dpng', '-r300');
print(fig, [outName '.eps'], '-depsc');
fprintf('Saved %s.png / .eps\n', outName);


%% ============================================================ local functions
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

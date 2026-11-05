%% Flight_Atlas.m
%  One page per flight day (Reviewer 2, per-flight documentation request):
%     (a) map of the upwind and all downwind transects, coloured by CH4, with
%         Marion County, the SSLF, the Indianapolis highways and the mean wind
%     (b) potential temperature profiles with the diagnosed ABL heights
%     (c) CH4 profiles with the same ABL heights
%     (d) table of per-transect fluxes for the four methods and the day means
%     (e...) one panel per downwind transect: CH4 along the transect (km along
%         the flight track, origin and sign convention as in Figure 2) with the
%         four backgrounds: Downwind Edge (blue), Upwind (orange), Adjusted ABL
%         Upwind (green), Upwind Edge (purple); the part of the transect directly
%         downwind of Marion County is shaded.
%
%  Backgrounds are rebuilt with the same code as the Method scripts (Downwind Edge as in
%  MassBalanceMethodA.m / plot_Figure2_Edge.m, Upwind as in MassBalanceMethodB.m /
%  plot_Figure3_5_Upwind.m).  The Adjusted ABL Upwind and Upwind Edge lines are
%  the Upwind line times the per-transect factors stored by Methods C and D.
%  Flux values and those factors are read from All_Solutions.mat
%  (made by Collect_Solutions.m), so every number on a page is the one in Table S1.
%
%  Output: Atlas/<Date>.png (one per day) and Atlas/Flight_Atlas.pdf (all days).
%
%  Needs on the path: INFLUX_<Date>, LegAndProfileList/<Date>.mat, FlightList.mat,
%  ElevationData.nc, InterX.m, alongTrackKm.m, PlotINFLUX.m (or the interstate
%  shapefile), and All_Solutions.mat in the current folder.

clear; close all force; clc;

%% ------------------------------------------------------------ settings
DATES    = {};               % {} = all days in FlightList; or e.g. {'20121004','20141119'}
MB_ROOT  = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';
LEG_DIR  = fullfile(MB_ROOT, 'LegAndProfileList');
OBS_DIR  = '/storage/work/zrb5027/WRFGUI_V3/Observations/INFLUX';      % INFLUX_<Date>
CODE_DIR = '/storage/work/zrb5027/WRFGUI_V3/Code';                     % InterX.m, PlotINFLUX.m
ROADS_SHP = '/storage/work/nvb5011/INFLUX_Paper_2/Thomas_2016/Shapefiles';
SOL_FILE = 'All_Solutions.mat';
OUT_DIR  = 'Atlas';
PDF_FILE = fullfile(OUT_DIR, 'Flight_Atlas.pdf');
SSLF_LATLON = [39.714691162109375, -86.205429077148438];
%BottomRight BottomLeft TopLeft TopRight
CountyCorners=[39.63760423971629, -85.95148383144084; 39.63337385377353, -86.32570560518701; 39.92446141560437, -86.32734193809195; 39.9271330073512, -85.95582253907553];

addpath(MB_ROOT); addpath(LEG_DIR); addpath(OBS_DIR); addpath(CODE_DIR);
if isfolder(ROADS_SHP), addpath(ROADS_SHP); end
for q = {'FlightList.mat','ElevationData.nc','InterX','alongTrackKm'}
    assert(~isempty(which(q{1})), 'Cannot find %s on the path.', q{1});
end
assert(isfile(SOL_FILE), '%s not found in %s (run Collect_Solutions.m first)', SOL_FILE, pwd);
if ~isfolder(OUT_DIR), mkdir(OUT_DIR); end
if isfile(PDF_FILE), delete(PDF_FILE); end

Sol = load(SOL_FILE);                          % SolA, SolB, SolC, SolD
FL = load('FlightList'); FlightList = FL.FlightList;
if isempty(DATES)
    DATES = arrayfun(@(f) f.name(1:8), FlightList, 'UniformOutput', false);
end
E.Elev = ncread('ElevationData.nc','data'); E.Elev(E.Elev<0) = 0;
E.Lat  = ncread('ElevationData.nc','lat');
E.Lon  = ncread('ElevationData.nc','lon'); E.Lon(E.Lon>180) = E.Lon(E.Lon>180) - 360;

col.A = [0.00 0.00 1.00];   col.B = [0.90 0.50 0.00];
col.C = [0.00 0.60 0.20];   col.D = [0.55 0.10 0.70];
names = {'Downwind Edge','Upwind','Adj. ABL Upwind','Upwind Edge'};

%% ------------------------------------------------------------ loop over days
for iDay = 1:numel(DATES)
    Date = DATES{iDay};
    try
        make_page(Date, Sol, E, CountyCorners, SSLF_LATLON, col, names, OUT_DIR, PDF_FILE);
        fprintf('%2d/%2d  %s  done\n', iDay, numel(DATES), Date);
    catch err
        fprintf('%2d/%2d  %s  ERROR: %s (line %d)\n', iDay, numel(DATES), Date, err.message, err.stack(1).line);
    end
    close all force;
end
fprintf('Atlas written to %s\n', PDF_FILE);


%% ============================================================ one page
function make_page(Date, Sol, E, CountyCorners, SSLF_LATLON, col, names, OUT_DIR, PDF_FILE)

%% ---- data (as in the Method scripts)
L = load(Date);                                % LegList, ProfileList, UpwindLeg?, CH4Flagged?, ManualABLHeights?
LegList = L.LegList; ProfileList = L.ProfileList;
hasUp = isfield(L,'UpwindLeg');
if hasUp, UpwindLeg = L.UpwindLeg(1):L.UpwindLeg(2); AllLeg = UpwindLeg; else, UpwindLeg = []; AllLeg = []; end
for k = 1:size(LegList,1), AllLeg = [AllLeg, LegList(k,1):LegList(k,2)]; end %#ok<AGROW>

x = load(['INFLUX_',Date]); x(x<-998) = NaN;
Time=x(:,1); Lat=x(:,2); Lon=x(:,3); Alt=x(:,4);
WSPD=x(:,11); WDIR=x(:,12);
U=-1*WSPD.*sind(WDIR); V=-1*WSPD.*cosd(WDIR);
CH4=movmean(x(:,6),30,'omitnan');
if isfield(L,'CH4Flagged'), CH4 = L.CH4Flagged; end
T=x(:,10); P=x(:,9);
MeanU=mean(U(AllLeg),'omitnan'); MeanV=mean(V(AllLeg),'omitnan');
MeanWSPD=(MeanU.^2+MeanV.^2).^.5; MeanWDIR=atan2d(MeanU,MeanV)+180;
ppbFac = 1; if median(CH4,'omitnan') < 100, ppbFac = 1000; end
CH4ppb = CH4*ppbFac;

Elevation = NaN(size(Alt));
for i = 1:length(Alt)
    Row = find(abs(E.Lon-Lon(i))<0.2, 1); Col = find(abs(E.Lat-Lat(i))<0.2, 1);
    if ~isempty(Row) && ~isempty(Col), Elevation(i) = E.Elev(Row,Col); end
end
AltAGL = Alt - Elevation;

% saved solutions for this day
SA = pick(Sol.SolA, Date); SB = pick(Sol.SolB, Date); SC = pick(Sol.SolC, Date); SD = pick(Sol.SolD, Date);
nT = size(LegList,1);

%% ---- profiles
nP = size(ProfileList,1);
Heights = 0:50:5000; Hedge = -25:50:5025;
prof = struct('t',{},'th',{},'ch4',{},'abl',{});
for p = 1:nP
    Profile = ProfileList(p,1):ProfileList(p,2);
    AltAGLX = AltAGL(Profile); Theta = T(Profile).*(1000./P(Profile)).^(2/7); CH4X = CH4ppb(Profile);
    th = NaN(1,numel(Heights)); ch = th;
    for h = 1:numel(Heights)
        g = AltAGLX>=Hedge(h) & AltAGLX<Hedge(h+1);
        th(h) = mean(Theta(g),'omitnan'); ch(h) = mean(CH4X(g),'omitnan');
    end
    if isfield(L,'ManualABLHeights')
        abl = L.ManualABLHeights(p);
    elseif ~isempty(SA) && isfield(SA,'ABLHeight') && numel(SA.ABLHeight) >= p
        abl = SA.ABLHeight(p);                 % what the Method scripts diagnosed
    else
        dTh = (th(1:end-2)-th(3:end))/100; dCH = (ch(1:end-2)-ch(3:end))/100/ppbFac;
        top = find(dCH>2e-5 & dTh<-5e-3, 1); if isempty(top), abl = NaN; else, abl = Heights(top+1); end
    end
    prof(p) = struct('t',mean(Time(Profile),'omitnan'),'th',th,'ch4',ch,'abl',abl);
end

%% ---- per-transect backgrounds
tr = struct('s',{},'ch4',{},'isDown',{},'bgA',{},'bgB',{},'bgC',{},'bgD',{},'t',{},'z',{},'zi',{});
for k = 1:nT
    Leg = LegList(k,1):LegList(k,2);
    [bgA, DownA, Inter, orient] = edge_background(Lat, Lon, CH4, Leg, MeanU, MeanV, CountyCorners);
    latO = mean(Inter(:,1)); lonO = mean(Inter(:,2));
    s = alongTrackKm(Lat(Leg), Lon(Leg), latO, lonO, orient);
    tr(k).s = s; tr(k).ch4 = CH4ppb(Leg); tr(k).isDown = ismember(Leg, DownA);
    tr(k).bgA = bgA(:)*ppbFac; tr(k).bgB = []; tr(k).bgC = []; tr(k).bgD = [];
    if hasUp
        bgB = upwind_background(Lat, Lon, CH4, Leg, UpwindLeg, MeanU, MeanV, CountyCorners);
        tr(k).bgB = bgB(:)*ppbFac;
        if ~isempty(SC) && isfield(SC,'CorrectionFactor') && numel(SC.CorrectionFactor) >= k
            tr(k).bgC = tr(k).bgB * SC.CorrectionFactor(k);
            if ~isempty(SD) && isfield(SD,'CorrectionFactor2') && numel(SD.CorrectionFactor2) >= k
                tr(k).bgD = tr(k).bgC * SD.CorrectionFactor2(k);
            end
        end
    end
    tr(k).t = mean(Time(Leg),'omitnan'); tr(k).z = mean(AltAGL(Leg),'omitnan');
    tr(k).zi = NaN; if ~isempty(SA) && numel(SA.DownwindABL) >= k, tr(k).zi = SA.DownwindABL(k); end
end

%% ---- figure
nRows = 1 + ceil(nT/4);
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 17 4.2*nRows],'Visible','off');
tl = tiledlayout(fig, nRows, 4, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf('%s/%s/%s    mean wind %.1f m s^{-1} from %.0f°    %d downwind transect(s)%s', ...
    Date(5:6), Date(7:8), Date(1:4), MeanWSPD, MeanWDIR, nT, ternary(hasUp, ', upwind transect flown', ', no upwind transect')), ...
    'FontSize', 14, 'FontWeight', 'bold');

% (a) map
ax = nexttile(tl); %#ok<NASGU>
worldmap([39.2 40.6],[-87 -85.2]); setm(gca,'FontSize',9); hold on
lim = [min(CH4ppb(AllLeg)) max(CH4ppb(AllLeg))];
if hasUp
    scatterm(Lat(UpwindLeg),Lon(UpwindLeg),60,[0 0 0],'filled');
    scatterm(Lat(UpwindLeg),Lon(UpwindLeg),14,CH4ppb(UpwindLeg),'filled');
end
for k = 1:nT
    Leg = LegList(k,1):LegList(k,2);
    scatterm(Lat(Leg),Lon(Leg),14,CH4ppb(Leg),'filled');
end
caxis(lim); colormap(gca, parula);
c = colorbar('eastoutside'); c.FontSize = 8; c.Label.String = 'CH_4 (ppb)';
draw_roads();
for k = 1:4, k2 = mod(k,4)+1; plotm(CountyCorners([k k2],1),CountyCorners([k k2],2),'LineWidth',2,'Color','k'); end
plotm(SSLF_LATLON(1), SSLF_LATLON(2), '^', 'MarkerSize', 8, 'MarkerFaceColor','k', 'MarkerEdgeColor','w');
quiverm(39.42, -85.55, 0.25*MeanV/MeanWSPD, 0.25*MeanU/MeanWSPD*(111/86), 'k', 0);
title('(a) Flight track', 'FontSize', 10);

% (b) theta profiles, (c) CH4 profiles
pc = lines(max(nP,1)); zTop = 0; for p = 1:nP, zTop = max(zTop, 1.3*max([prof(p).abl 500])); end
zTop = min(max(zTop, 1500), 5000);
for panel = 1:2
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for p = 1:nP
        if panel == 1, v = prof(p).th; else, v = prof(p).ch4; end
        plot(ax, v, Heights, '-', 'Color', pc(p,:), 'LineWidth', 1.3, 'DisplayName', sprintf('Profile %d, %s UTC', p, hhmm(prof(p).t)));
        if ~isnan(prof(p).abl), yline(ax, prof(p).abl, '--', 'Color', pc(p,:), 'LineWidth', 1, 'HandleVisibility','off'); end
    end
    ylim(ax, [0 zTop]); ylabel(ax, 'Height AGL (m)');
    if panel == 1
        xlabel(ax, 'Potential temperature (K)'); title(ax, '(b) \theta profiles (dashed: ABL height)', 'FontSize', 10);
        legend(ax, 'Location','southeast', 'FontSize', 8);
    else
        xlabel(ax, 'CH_4 (ppb)'); title(ax, '(c) CH_4 profiles', 'FontSize', 10);
    end
    set(ax, 'FontSize', 9);
end

% (d) flux table
ax = nexttile(tl); axis(ax,'off');
lines_ = {sprintf('%-3s %-6s %5s %5s %7s %7s %7s %7s', 'T', 'UTC', 'z m', 'zi m', 'DwEdge', 'Upwind', 'AdjABL', 'UpEdge')};
fl = NaN(nT,4);
for k = 1:nT
    fl(k,:) = [fluxOf(SA,k) fluxOf(SB,k) fluxOf(SC,k) fluxOf(SD,k)];
    lines_{end+1} = sprintf('%-3d %-6s %5.0f %5.0f %7s %7s %7s %7s', k, hhmm(tr(k).t), tr(k).z, tr(k).zi, n0(fl(k,1)), n0(fl(k,2)), n0(fl(k,3)), n0(fl(k,4))); %#ok<AGROW>
end
m = mean(fl,1,'omitnan');
lines_{end+1} = sprintf('%-3s %-6s %5s %5s %7s %7s %7s %7s', 'mean', '', '', '', n0(m(1)), n0(m(2)), n0(m(3)), n0(m(4)));
if ~isempty(SA) && isfield(SA,'CosTheta')
    lines_{end+1} = sprintf('cos(theta) %.2f to %.2f', min(SA.CosTheta), max(SA.CosTheta));
end
text(ax, 0, 1, strjoin(lines_, newline), 'FontName','FixedWidth', 'FontSize', 9, 'VerticalAlignment','top', 'Interpreter','none');
title(ax, '(d) Flux by method (mol s^{-1})', 'FontSize', 10);
% one legend per page for the transect panels, drawn here with a fixed location
% (a legend placed with 'best' inside an invisible tiled figure can end up outside its panel)
hold(ax, 'on');
hl = gobjects(1, 6);
hl(1) = plot(ax, NaN, NaN, 'k-', 'LineWidth', 1.4);
hl(2) = plot(ax, NaN, NaN, '-', 'Color', col.A, 'LineWidth', 1.3);
hl(3) = plot(ax, NaN, NaN, '-', 'Color', col.B, 'LineWidth', 1.1);
hl(4) = plot(ax, NaN, NaN, '-', 'Color', col.C, 'LineWidth', 1.1);
hl(5) = plot(ax, NaN, NaN, '-', 'Color', col.D, 'LineWidth', 1.1);
hl(6) = patch(ax, [NaN NaN NaN], [NaN NaN NaN], [0.93 0.93 0.93], 'EdgeColor', 'none');
lab = [{'aircraft CH_4'}, names, {'directly downwind of Marion County'}];
keepL = [true true hasUp hasUp hasUp true];
legend(ax, hl(keepL), lab(keepL), 'Location', 'south', 'FontSize', 8, 'NumColumns', 2, 'Box', 'on');
xlim(ax, [0 1]); ylim(ax, [0 1]);

% (e...) transects
letters = 'efghijklmnopqrstuvwxyz';
for k = 1:nT
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
    s = tr(k).s; y = tr(k).ch4; d = tr(k).isDown;
    yl = [min([y; tr(k).bgA; tr(k).bgB; tr(k).bgC; tr(k).bgD]) max([y; tr(k).bgA; tr(k).bgB; tr(k).bgC; tr(k).bgD])];
    yl = yl + [-0.05 0.05]*diff(yl);
    if any(d)
        patch(ax, [min(s(d)) max(s(d)) max(s(d)) min(s(d))], [yl(1) yl(1) yl(2) yl(2)], [0.93 0.93 0.93], 'EdgeColor','none');
    end
    h = gobjects(1,4);
    h(1) = plot(ax, s, tr(k).bgA, '-', 'Color', col.A, 'LineWidth', 1.3);
    if ~isempty(tr(k).bgB), h(2) = plot(ax, s, tr(k).bgB, '-', 'Color', col.B, 'LineWidth', 1.1); end
    if ~isempty(tr(k).bgC), h(3) = plot(ax, s, tr(k).bgC, '-', 'Color', col.C, 'LineWidth', 1.1); end
    if ~isempty(tr(k).bgD), h(4) = plot(ax, s, tr(k).bgD, '-', 'Color', col.D, 'LineWidth', 1.1); end
    plot(ax, s, y, 'k-', 'LineWidth', 1.4);
    xlim(ax, [min(s) max(s)]); ylim(ax, yl); set(ax, 'FontSize', 9);
    xlabel(ax, 'Distance along flight track (km)'); if mod(k-1,4) == 0, ylabel(ax, 'CH_4 (ppb)'); end
    if any(~isnan(fl(k,2:4)))          % abbreviations as in the column headers of table (d)
        txt = sprintf('DwEdge %s, Upwind %s, AdjABL %s, UpEdge %s mol s^{-1}', n0(fl(k,1)), n0(fl(k,2)), n0(fl(k,3)), n0(fl(k,4)));
    else
        txt = sprintf('%s %s mol s^{-1}', names{1}, n0(fl(k,1)));
    end
    title(ax, {sprintf('(%s) T%d  %s UTC  %.0f m AGL', letters(k), k, hhmm(tr(k).t), tr(k).z), txt}, 'FontSize', 9);
end

%% ---- save
exportgraphics(fig, fullfile(OUT_DIR, [Date '.png']), 'Resolution', 150);
exportgraphics(fig, PDF_FILE, 'Append', isfile(PDF_FILE), 'ContentType', 'image', 'Resolution', 150);
close(fig);
end


%% ============================================================ Edge background (Method A, verbatim)
function [BackgroundLine, Downwind, Intersection, orient] = edge_background(Lat, Lon, CH4, Leg, MeanU, MeanV, CountyCorners)
if (MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0) %#ok<*OR2,*AND2>
    PointB1=CountyCorners(3,:); PointB2=CountyCorners(1,:);
else
    PointB1=CountyCorners(2,:); PointB2=CountyCorners(4,:);
end
PointC=[PointB1(1)+MeanV, PointB1(2)+MeanU*(111/86)]; PointA=[PointB1(1)-MeanV, PointB1(2)-MeanU*(111/86)];
Line1a=[PointA(2) PointB1(2) PointC(2);PointA(1) PointB1(1) PointC(1)];
PointC=[PointB2(1)+MeanV, PointB2(2)+MeanU*(111/86)]; PointA=[PointB2(1)-MeanV, PointB2(2)-MeanU*(111/86)];
Line1b=[PointA(2) PointB2(2) PointC(2);PointA(1) PointB2(1) PointC(1)];
Line2=[Lon(Leg)';Lat(Leg)'];
Intersection1=flipud(InterX(Line1a,Line2))'; Intersection2=flipud(InterX(Line1b,Line2))';
if ((MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0))
    if isempty(Intersection1), Intersection1=[max(Lat(Leg)) min(Lon(Leg))]; end
    if isempty(Intersection2), Intersection2=[min(Lat(Leg)) max(Lon(Leg))]; end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05,     Upwind1=find(Lon(Leg)<Intersection1(2));
    elseif abs(Intersection1(2)-Intersection2(2)) < 0.05, Upwind1=find(Lat(Leg)>Intersection1(1));
    else,                                                 Upwind1=find(Lat(Leg)>Intersection1(1) | Lon(Leg)<Intersection1(2)); end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05,     Upwind2=find(Lon(Leg)>Intersection2(2));
    elseif abs(Intersection2(2)-Intersection2(2)) < 0.05, Upwind2=find(Lat(Leg)<Intersection2(1));
    else,                                                 Upwind2=find(Lat(Leg)<Intersection2(1) | Lon(Leg)>Intersection2(2)); end
else
    if isempty(Intersection1), Intersection1=[min(Lat(Leg)) min(Lon(Leg))]; end
    if isempty(Intersection2), Intersection2=[max(Lat(Leg)) max(Lon(Leg))]; end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05,     Upwind1=find(Lon(Leg)<Intersection1(2));
    elseif abs(Intersection1(2)-Intersection2(2)) < 0.05, Upwind1=find(Lat(Leg)<Intersection1(1));
    else,                                                 Upwind1=find(Lat(Leg)<Intersection1(1) | Lon(Leg)<Intersection1(2)); end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05,     Upwind2=find(Lon(Leg)>Intersection2(2));
    elseif abs(Intersection2(2)-Intersection2(2)) < 0.05, Upwind2=find(Lat(Leg)>Intersection2(1));
    else,                                                 Upwind2=find(Lat(Leg)>Intersection2(1) | Lon(Leg)>Intersection2(2)); end
end
Intersection=[Intersection1;Intersection2];
Upwind1=Leg(Upwind1); Upwind2=Leg(Upwind2);
Downwind=Leg(~ismember(Leg,[Upwind1 Upwind2]));
Upwind1Background=median(CH4(Upwind1),'omitnan'); Upwind2Background=median(CH4(Upwind2),'omitnan');
mTemp=round(length(Leg)/10);
if isempty(Upwind1)
    if Leg(1)==Upwind2(1), Upwind1Background=median(CH4(Leg(end-mTemp:end)),'omitnan'); else, Upwind1Background=median(CH4(Leg(1:mTemp)),'omitnan'); end
elseif isempty(Upwind2)
    if Leg(1)==Upwind1(1), Upwind2Background=median(CH4(Leg(end-mTemp:end)),'omitnan'); else, Upwind2Background=median(CH4(Leg(1:mTemp)),'omitnan'); end
end
LatDistance=distance(min(Lat(Leg)),mean(Lon(Leg),'omitnan'),max(Lat(Leg)),mean(Lon(Leg),'omitnan'));
LonDistance=distance(mean(Lat(Leg),'omitnan'),min(Lon(Leg)),mean(Lat(Leg),'omitnan'),max(Lon(Leg)));
if LonDistance>LatDistance
    if Lon(Leg(1))>Lon(Leg(end)), [Upwind1Background,Upwind2Background]=deal(Upwind2Background,Upwind1Background); end
    LineCoeff=polyfit([Lon(Leg(1)),Lon(Leg(end))],[Upwind1Background,Upwind2Background],1);
    BackgroundLine=LineCoeff(1).*Lon(Leg)+LineCoeff(2); orient='east';
else
    if (MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0), [Upwind1Background,Upwind2Background]=deal(Upwind2Background,Upwind1Background); end
    if Lat(Leg(1))>Lat(Leg(end)), [Upwind1Background,Upwind2Background]=deal(Upwind2Background,Upwind1Background); end
    LineCoeff=polyfit([Lat(Leg(1)),Lat(Leg(end))],[Upwind1Background,Upwind2Background],1);
    BackgroundLine=LineCoeff(1).*Lat(Leg)+LineCoeff(2); orient='north';
end
end

%% ============================================================ Upwind background (Method B, verbatim)
function BackgroundLine = upwind_background(Lat, Lon, CH4, Leg, UpwindLeg, MeanU, MeanV, ~)
UpwindLine=[Lon(UpwindLeg), Lat(UpwindLeg)];
BackgroundLine=NaN(1,length(Leg));
for UpwindLoop=1:length(Leg)
    PointB=[Lat(Leg(UpwindLoop)), Lon(Leg(UpwindLoop))];
    PointC=[PointB(1)+MeanV, PointB(2)+MeanU*(111/86)];
    PointA=[PointB(1)-MeanV, PointB(2)-MeanU*(111/86)];
    Line1A=[PointA(2) PointB(2) PointC(2);PointA(1) PointB(1) PointC(1)];
    UpwindIntersection=flipud(InterX(Line1A,UpwindLine'))';
    if ~isempty(UpwindIntersection)
        TempDistance=distance(UpwindIntersection(1),UpwindIntersection(2),UpwindLine(:,2),UpwindLine(:,1));
        GoodSpot=find(TempDistance==min(TempDistance));
        BackgroundLine(UpwindLoop)=CH4(UpwindLeg(GoodSpot));
    end
end
BackgroundLine=fill_gaps(BackgroundLine);
DownwindLine=CH4(Leg);
R=NaN(1,101); ShiftCounter=0;
for Shift=-50:50
    ShiftCounter=ShiftCounter+1;
    TempBackground=NaN(length(BackgroundLine)+100,1); TempDownwind=NaN(length(BackgroundLine)+100,1);
    TempDownwind(51:end-50)=DownwindLine; TempBackground(51+Shift:end-50+Shift)=BackgroundLine;
    TempCorr=corrcoef(TempBackground,TempDownwind,'Rows','complete'); R(ShiftCounter)=TempCorr(2);
end
Shift=find(R==max(R),1)-51;
BackgroundLineDoppel=NaN(length(BackgroundLine)+100,1);
BackgroundLineDoppel(51+Shift:end-50+Shift)=BackgroundLine;
BackgroundLine=BackgroundLineDoppel(51:end-50)';
BackgroundLine=fill_gaps(BackgroundLine);
BackgroundLine(isnan(BackgroundLine))=median(BackgroundLine,'omitnan');
end

%% ============================================================ helpers
function y = fill_gaps(y)
for w = [120 240 360], Temp = movmean(y, w, 'omitnan'); Bad = isnan(y); y(Bad) = Temp(Bad); end
end
function S = pick(SolX, Date)
S = []; if isempty(SolX), return; end
k = find(strcmp({SolX.Date}, Date), 1); if ~isempty(k), S = SolX(k); end
end
function v = fluxOf(S, k)
v = NaN; if ~isempty(S) && isfield(S,'Flux') && numel(S.Flux) >= k, v = S.Flux(k); end
end
function s = n0(v)
if isnan(v), s = 'n/a'; else, s = sprintf('%.0f', v); end
end
function s = hhmm(t)
s = sprintf('%02d:%02d', floor(t/3600), round(mod(t,3600)/60));
end
function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
function draw_roads()
if ~isempty(which('PlotINFLUX')), PlotINFLUX; return; end
if isempty(which('Indianapolis_Interstate_Mileposts.shp')), return; end
wsv = warning('off','all');
S  = shaperead('Indianapolis_Interstate_Mileposts.shp');
ms = defaultm('tranmerc'); ms.geoid = almanac('earth','grs80','survey foot');
ms.origin = [37.5 -85.66666666666667 0]; ms.falseeasting = 328083.3333333333;
ms.falsenorthing = 820208.3333333331; ms.scalefactor = 0.9999666666666667; ms = defaultm(ms);
[LATr, LONr] = minvtran(ms, [S.X], [S.Y]);
scatterm(LATr, LONr, 4, [0.25 0.25 0.25], 'filled');
warning(wsv);
end

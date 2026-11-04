%% Figure12_v2.m
%  Revised Figure 12 (Figure 13 in the submission): WRF-Chem CH4 enhancement
%  maps for 11/19/2014 at 16, 17 and 18 UTC (top) and model vs. aircraft along
%  the contemporaneous flight segments (bottom).
%
%  Changes from CoalMineCase_6panel_hours (the submitted figure):
%    R1.12  colour scale and lower panels in ppb instead of ppm
%    R2.21  coal mine locations (triangles), latitude/longitude grid with
%           labels, and a 100 km scale bar on the maps
%    R1.6   lower panels keep UTC hour on the bottom axis (the segments are
%           closed loops, so time is the natural coordinate) and gain a top
%           axis with the distance flown along the track in km
%  Everything else (data files, model sampling at the aircraft positions, the
%  constant offset dx applied to the observations) is unchanged.
%
%  Mine coordinates: GHGRP / Global Energy Monitor (Sugar Camp = MC#1, Hamilton
%  County Coal No. 1, Pond Creek = Mach #1, New Era and New Future at Galatia,
%  Pattiki); set SHOW_MINES = false to drop them.

clear; close all

addpath('/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/UserSaveData/2014-11-19')
addpath('/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/Code')
addpath('/storage/work/nvb5011/INFLUX_Paper_2/Thomas_2016/Shapefiles')

dx = -0.010;                      % ppm offset applied to the observations (as before)
PPB = 1000;                       % ppm -> ppb
CMAX = 80;                        % colour scale top (ppb); was 0.08 ppm
YL  = [40 90];                    % lower panels (ppb); was [0.04 0.09] ppm
SHOW_MINES = true;
% southern Illinois underground coal mines (lat, lon), 2014 GHGRP reporters
mines = [38.0264 -88.7679;        % Sugar Camp / MC#1 (Franklin Co.)
         38.1702 -88.6014;        % Hamilton County Coal No. 1 (Hamilton Co.)
         37.8400 -88.8300;        % Pond Creek / Mach #1 (Williamson Co.)
         37.8643 -88.6145;        % New Era (Saline Co.)
         37.8297 -88.6399;        % New Future (Saline Co.)
         38.0648 -88.0896];       % Pattiki (White Co.)
MARION = [39.7806 -86.1401];
dMine = distance(mean(mines(:,1)), mean(mines(:,2)), MARION(1), MARION(2))*111.19;
fprintf('Mine cluster centre to Marion County centre: %.0f km\n', dMine);

mapFiles = {'MapVariables-16ZINDY3_16.mat', ...
            'MapVariables-17ZINDY3_17.mat', ...
            'MapVariables-18ZINDY3_18.mat'};
tsFiles  = {'TimeSeriesVariablesCOAL16.mat', ...
            'TimeSeriesVariablesCOAL17.mat', ...
            'TimeSeriesVariablesCOAL18.mat'};
titles  = {'16 UTC','17 UTC','18 UTC'};
letters = {'(a)','(b)','(c)','(d)','(e)','(f)'};

%% ------------------------------------------------------------ shapefiles
States = shaperead('usastatelo','UseGeoCoords',true);
MCB    = shaperead('Marion_County_Boundary.shp');
S      = shaperead('Indianapolis_Interstate_Mileposts.shp');
mstruct = defaultm('tranmerc');
mstruct.geoid = almanac('earth','grs80','survey foot');
mstruct.origin = [37.5 -85.66666666666667 0];
mstruct.falseeasting = 328083.3333333333;
mstruct.falsenorthing = 820208.3333333331;
mstruct.scalefactor = 0.9999666666666667;
mstruct = defaultm(mstruct);
XX = nan(1,numel(S)); YY = nan(1,numel(S));
for i = 1:numel(S), XX(i) = S(i).X; YY(i) = S(i).Y; end
[RoadLat, RoadLon] = minvtran(mstruct, XX, YY);

%% ------------------------------------------------------------ layout
fig = figure('Color','w','Position',[100 100 1120 720]);
topPos = [0.08  0.58 0.23 0.30;
          0.385 0.58 0.23 0.30;
          0.69  0.58 0.23 0.30];
botPos = [0.08  0.11 0.23 0.26;
          0.385 0.11 0.23 0.26;
          0.69  0.11 0.23 0.26];
mapAx  = gobjects(1,3);
lineAx = gobjects(1,3);

for k = 1:3
    M = load(mapFiles{k});
    T = load(tsFiles{k});

    % ===================================================== map
    subplot(2,3,k)
    worldmap([37 41],[-90.3 -84.5]);
    set(gca,'Position',topPos(k,:))
    surfm(M.ModelLat, M.ModelLon, M.WRFVariable*PPB)
    hold on
    shading flat
    colormap(parula(256))
    clim([0 CMAX])

    geoshow(States,'FaceColor','none','EdgeColor',[0.35 0.35 0.35],'LineWidth',0.7)
    scatterm(RoadLat, RoadLon, 8, '.', 'MarkerEdgeColor', [0.45 0.45 0.45])
    plotm(MCB.Y, MCB.X, 'k', 'LineWidth', 1.0)

    scatterm(M.ObsLat, M.ObsLon, 90, 'k', 'filled')
    scatterm(M.ObsLat, M.ObsLon, 42, (M.ObsMeasurement + dx)*PPB, 'filled')

    WindRes = 15;
    [m,n] = size(M.ModelLat);
    z  = 1:WindRes:m;  zz = 1:WindRes:n;
    q = quiverm(M.ModelLat(z,zz), M.ModelLon(z,zz), M.V(z,zz), M.U(z,zz));
    set(q,'Color',[1 1 1],'LineWidth',0.9)

    if SHOW_MINES
        plotm(mines(:,1), mines(:,2), '^', 'MarkerSize', 7, 'MarkerFaceColor', [0.95 0.15 0.15], ...
              'MarkerEdgeColor', 'k', 'LineWidth', 0.7)
        if k == 1
            textm(37.55, -88.45, 'coal mines', 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.95 0.15 0.15], ...
                  'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'Margin', 1)
        end
    end

    % latitude / longitude grid with labels (R2.21)
    setm(gca, 'Grid','on', 'GLineStyle',':', 'GColor',[0.25 0.25 0.25], 'GLineWidth',0.5, ...
         'MLineLocation',2, 'PLineLocation',1, ...
         'MeridianLabel','on', 'ParallelLabel','on', 'MLabelLocation',2, 'PLabelLocation',1, ...
         'MLabelParallel','south', 'PLabelMeridian','west', 'LabelFormat','compass', 'FontSize',8, ...
         'Frame','on', 'FEdgeColor','k', 'FLineWidth',1.0);

    % 100 km scale bar (R2.21), bottom right of each map
    lat0 = 37.35; lon0 = -86.9;
    dlon = 100/(111.32*cosd(lat0));
    plotm([lat0 lat0], [lon0 lon0+dlon], 'k-', 'LineWidth', 3)
    plotm([lat0 lat0], [lon0 lon0+dlon], 'w-', 'LineWidth', 1)
    textm(lat0+0.12, lon0+dlon/2, '100 km', 'FontSize', 8, 'HorizontalAlignment', 'center', ...
          'BackgroundColor', 'w', 'Margin', 0.5)

    title(titles{k}, 'FontSize', 13, 'FontWeight', 'bold')
    mapAx(k) = gca;

    % ===================================================== lower panels
    lineAx(k) = axes('Position', botPos(k,:));
    hold(lineAx(k),'on')
    obs_plot = (M.ObsMeasurement(:) + dx)*PPB;
    model_plot = griddata(M.ModelLon(:), M.ModelLat(:), M.WRFVariable(:), M.ObsLon(:), M.ObsLat(:), 'linear');
    nanidx = isnan(model_plot);
    if any(nanidx)
        model_plot(nanidx) = griddata(M.ModelLon(:), M.ModelLat(:), M.WRFVariable(:), ...
                                      M.ObsLon(nanidx), M.ObsLat(nanidx), 'nearest');
    end
    model_plot = model_plot*PPB;
    xobs = T.ObsSecs(:) ./ 3600;

    plot(xobs, obs_plot, '-o', 'Color', [0.15 0.55 0.80], 'MarkerFaceColor', [0.15 0.55 0.80], ...
         'MarkerEdgeColor', 'k', 'MarkerSize', 3.5, 'LineWidth', 1.0)
    plot(xobs, model_plot, '-', 'Color', 'r', 'LineWidth', 1.8)
    box on; grid on
    set(gca, 'FontSize', 10, 'LineWidth', 1.0, 'TickDir', 'out', 'FontName', 'Arial')
    xlabel('Hour (UTC)', 'FontSize', 11)
    if k == 1, ylabel('CH_4 enhancement (ppb)', 'FontSize', 11); end
    xlim([min(xobs) max(xobs)]); ylim(YL)
    if k == 1
        legend({'Observations','Model at observation locations'}, 'Location','northwest', 'Box','off', 'FontSize',9)
    end

    % top axis: distance flown along the track (km), ticks placed at the
    % times when the aircraft had flown that far
    la = M.ObsLat(:); lo = M.ObsLon(:);
    seg = distance(la(1:end-1), lo(1:end-1), la(2:end), lo(2:end))*111.19;   % km between samples
    seg(isnan(seg)) = 0;
    dist = [0; cumsum(seg)];
    stepKm = 20; if max(dist) > 130, stepKm = 40; end               % at most ~6 ticks
    kmTicks = 0:stepKm:floor(max(dist)/stepKm)*stepKm;
    [dU, iu] = unique(dist);
    tTicks = interp1(dU, xobs(iu), kmTicks);
    ax2 = axes('Position', botPos(k,:), 'Color','none', 'XAxisLocation','top', 'YAxisLocation','right', ...
               'YTick', [], 'XLim', xlim(lineAx(k)), 'XTick', tTicks, 'XTickLabel', string(kmTicks), ...
               'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'Box', 'off', 'YColor', 'none');
    xlabel(ax2, 'Distance flown along the track (km)', 'FontSize', 9)
    fprintf('%s: segment %.1f to %.1f UTC, %.0f km flown\n', titles{k}, min(xobs), max(xobs), max(dist));
end

%% ------------------------------------------------------------ colour bar
cb = colorbar(mapAx(3), 'southoutside');
cb.Position = [0.30 0.48 0.40 0.022];
cb.Ticks = 0:10:CMAX;
cb.FontSize = 9;
cb.LineWidth = 0.8;
cb.Label.String = 'CH_4 enhancement (ppb)';
cb.Label.FontSize = 12;
cb.Label.FontWeight = 'bold';
for k = 1:3
    colormap(mapAx(k), parula(256))
    clim(mapAx(k), [0 CMAX])
end

%% ------------------------------------------------------------ panel letters
for k = 1:3
    annotation('textbox', [topPos(k,1)-0.060, topPos(k,2)+topPos(k,4)-0.010, 0.04, 0.03], ...
        'String', letters{k}, 'LineStyle','none', 'FontSize',14, 'FontWeight','bold', 'FontName','Arial');
    annotation('textbox', [botPos(k,1)-0.060, botPos(k,2)+botPos(k,4)+0.02, 0.04, 0.03], ...
        'String', letters{k+3}, 'LineStyle','none', 'FontSize',14, 'FontWeight','bold', 'FontName','Arial');
end

%% ------------------------------------------------------------ save
set(fig,'InvertHardcopy','off')
exportgraphics(fig,'Figure12_v2_CoalMineCase.png','Resolution',300)
fprintf('Saved Figure12_v2_CoalMineCase.png\n');

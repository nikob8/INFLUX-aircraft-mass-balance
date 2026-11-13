%% plot_Figure6_emissions.m
%  Figure 6: CH4 emissions used in the WRF-Chem simulations (Section 2.6).
%    (a) 9 km domain and (b) 3 km domain: sum of the five regional tracers (oil and gas,
%        coal, agriculture, landfills, other) from the EPA gridded inventory
%        (Maasakkers et al., 2016); Marion County is excluded from these grids.
%    (c) Marion County at 1 km (Lamb et al., 2016): natural gas, landfill and wastewater.
%
%  Inputs: the text files of data/emissions/ of the repository (columns latitude,
%  longitude, emission of the grid cell in kg CH4 per hour; converted here to mol/s):
%    Emissions9km/INFLUX_<sector>9, Emissions3km/INFLUX_<sector>3,
%    ThomasEmissions/ThomasGas, ThomasTrash (landfill), ThomasWastewater.
%  In panel (c) the landfill file is multiplied by 1.4 and the wastewater file divided
%  by 3.5, as in the figure of the paper.
%  Optional: Indianapolis_Interstate_Mileposts.shp (roads of panel c, not in the repository).
%  Needs the Mapping Toolbox; clim needs MATLAB R2022a or later (caxis in older versions).
%  Output: INFLUX_emissions_log10_turbo.png

clc
clear
close all

%% ===================== PATHS =====================
DATA_DIR  = fullfile('..', 'data', 'emissions');      % folder with Emissions9km, Emissions3km, ThomasEmissions
ROADS_SHP = 'Indianapolis_Interstate_Mileposts.shp';  % optional
addpath(fullfile(DATA_DIR, 'Emissions9km'))
addpath(fullfile(DATA_DIR, 'ThomasEmissions'))
addpath(fullfile(DATA_DIR, 'Emissions3km'))

%% ===================== STYLE =====================
fontName = 'Arial';
titleFS  = 15;
tickFS   = 10;
stateLW  = 1.15;
lineW    = 1.0;

roadMS_outer = 4.0;   % subtle grey underlay
roadMS_inner = 2.4;   % white roads on top
marker9      = 18;
marker3      = 16;

cmap = turbo(256);
bgGrey = [0.82 0.82 0.82];

frameLW = 1.6;
labelX  = 0.015;
labelY  = 0.985;

% log-scale plotting limits in physical units
zmin9 = 0.1;   zmax9 = 30;
zmin3 = 0.1;   zmax3 = 15;
zmin1 = 0.05;  zmax1 = 3;

%% ===================== INDIANAPOLIS ROADS =====================
haveRoads = isfile(ROADS_SHP);
if haveRoads
S = shaperead(ROADS_SHP);

mstruct = defaultm('tranmerc');
mstruct.geoid = almanac('earth','grs80','survey foot');
mstruct.origin = [37.5 -85.66666666666667 0];
mstruct.falseeasting = 328083.3333333333;
mstruct.falsenorthing = 820208.3333333331;
mstruct.scalefactor = 0.9999666666666667;
mstruct = defaultm(mstruct);

XX = zeros(1,1469);
YY = zeros(1,1469);
for i = 1:1469
    XX(i) = S(i).X;
    YY(i) = S(i).Y;
end

[LAT_road, LON_road] = minvtran(mstruct, XX, YY);
end

%% ===================== THOMAS / MARION COUNTY INVENTORIES =====================
A = importdata('ThomasGas');
lat_mc = A(:,1);
lon_mc = A(:,2);
TG = A(:,3);

lat_mesh = meshgrid(lat_mc(1:34:1156),1:34);
lon_mesh = meshgrid(lon_mc(1:34),1:34);
TG_mesh = reshape(TG,34,34);
TG_mesh2 = TG_mesh * 1000 * 1/3600 * 1/16;

B = importdata('ThomasTrash');
TT = B(:,3) * 1.4;
TT_mesh = reshape(TT,34,34);
TT_mesh2 = TT_mesh * 1000 * 1/3600 * 1/16;

C = importdata('ThomasWastewater');
lat1 = C(:,1);
lon1 = C(:,2);
TW = C(:,3);
TW_mesh = reshape(TW,34,34);
TW_mesh2 = TW_mesh/3 * 1000 * 1/3600 * 1/16;

T_total  = TW_mesh/3.5 + TT_mesh + TG_mesh;
T_total2 = T_total * 1000 * 1/3600 * 1/16;

%% ===================== 3-km INVENTORY =====================
A = importdata('INFLUX_OilAndGas3');
lat3 = A(:,1);
lon3 = A(:,2);
EG3  = A(:,3);

B = importdata('INFLUX_Landfill3');
ET3  = B(:,3);

C = importdata('INFLUX_Coal3');
EC3  = C(:,3);

D = importdata('INFLUX_Ag3');
EA3  = D(:,3);

E = importdata('INFLUX_Other3');
EO3  = E(:,3);

total_E3 = EO3 + EA3 + EC3 + ET3 + EG3;
total_E3 = total_E3 * 1000 * 1/3600 * 1/16;   % kg/hr -> mol/s

%% ===================== 9-km INVENTORY =====================
A = importdata('INFLUX_OilAndGas9');
lat9 = A(:,1);
lon9 = A(:,2);
EG9  = A(:,3);

B = importdata('INFLUX_Landfill9');
ET9  = B(:,3);

C = importdata('INFLUX_Coal9');
EC9  = C(:,3);

D = importdata('INFLUX_Ag9');
EA9  = D(:,3);

E = importdata('INFLUX_Other9');
EO9  = E(:,3);

total_E9 = EO9 + EA9 + EC9 + ET9 + EG9;
total_E9 = total_E9 * 1000 * 1/3600 * 1/16;   % kg/hr -> mol/s

%% ===================== LOG SCALE DATA =====================
% Panels (a) and (b): keep zero / excluded cells unplotted
valid9 = total_E9 > 0;
valid3 = total_E3 > 0;

total_E9_log = log10(max(total_E9(valid9), zmin9));
total_E3_log = log10(max(total_E3(valid3), zmin3));

% Panel (c): keep smooth filled appearance
T_total2_log = log10(max(T_total2, zmin1));

%% ===================== STATES SHAPEFILE =====================
States = shaperead('usastatelo','UseGeoCoords',true);

%% ===================== FIGURE =====================
figure('Color','w','Position',[70 110 1550 560])

%% ===================== PANEL (a): 9-km =====================
subplot(1,3,1)
worldmap([min(lat9) max(lat9)], [min(lon9) max(lon9)]);

setm(gca,'FFaceColor',bgGrey)
hold on

scatterm(lat9(valid9), lon9(valid9), marker9, total_E9_log, ...
    'filled', 's', 'MarkerEdgeColor', 'none')

geoshow(States, ...
    'FaceColor','none', ...
    'EdgeColor',[1 1 1], ...
    'LineWidth',stateLW)

gridm('off')
mlabel off
plabel off
tightmap

colormap(gca, cmap)
clim(log10([zmin9 zmax9]))

title({'9-km WRF-Chem'; ...
       'Domain with CH_4 Emissions (mol s^{-1})'}, ...
       'FontSize', titleFS, ...
       'FontWeight', 'bold', ...
       'FontName', fontName)

text(0.14, 0.985, '(a)', ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'Color','w', ...
    'FontSize',23, ...
    'FontWeight','bold')

%framem on
%setm(gca, 'Frame', 'on', 'FEdgeColor', 'k', 'FLineWidth', 2.2)
%hFrame = handlem('Frame');
%uistack(hFrame, 'top')

cb1 = colorbar;
cb1.FontSize = tickFS;
cb1.TickDirection = 'out';
cb1.LineWidth = 0.8;
cb1.Ticks = log10([0.1 0.3 1 3 10 30]);
cb1.TickLabels = {'0.1','0.3','1','3','10','30'};

framem on
setm(gca, 'Frame', 'on', 'FEdgeColor', 'k', 'FLineWidth', 2.0)

drawnow
hF = handlem('Frame');

for ii = 1:numel(hF)
    if isprop(hF(ii), 'FaceColor')
        set(hF(ii), 'FaceColor', 'none')
    end
    if isprop(hF(ii), 'EdgeColor')
        set(hF(ii), 'EdgeColor', 'k')
    end
    if isprop(hF(ii), 'Color')
        set(hF(ii), 'Color', 'k')
    end
    if isprop(hF(ii), 'LineWidth')
        set(hF(ii), 'LineWidth', 2.0)
    end
end

uistack(hF, 'top')

%% ===================== PANEL (b): 3-km =====================
subplot(1,3,2)
worldmap([min(lat3) max(lat3)], [min(lon3) max(lon3)]);


setm(gca,'FFaceColor',bgGrey)
hold on

scatterm(lat3(valid3), lon3(valid3), marker3, total_E3_log, ...
    'filled', 's', 'MarkerEdgeColor', 'none')

geoshow(States, ...
    'FaceColor','none', ...
    'EdgeColor',[1 1 1], ...
    'LineWidth',stateLW)

gridm('off')
mlabel off
plabel off
tightmap

colormap(gca, cmap)
clim(log10([zmin3 zmax3]))

title({'3-km WRF-Chem'; ...
       'Domain with CH_4 Emissions (mol s^{-1})'}, ...
       'FontSize', titleFS, ...
       'FontWeight', 'bold', ...
       'FontName', fontName)

text(0.07, 0.985, '(b)', ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'Color','w', ...
    'FontSize',23, ...
    'FontWeight','bold')

%framem on
%setm(gca, 'Frame', 'on', 'FEdgeColor', 'k', 'FLineWidth', 2.2)
%hFrame = handlem('Frame');
%uistack(hFrame, 'top')

cb2 = colorbar;
cb2.FontSize = tickFS;
cb2.TickDirection = 'out';
cb2.LineWidth = 0.8;
cb2.Ticks = log10([0.1 0.3 1 3 10 15]);
cb2.TickLabels = {'0.1','0.3','1','3','10','15'};

framem on
setm(gca, 'Frame', 'on', 'FEdgeColor', 'k', 'FLineWidth', 2.0)

drawnow
hF = handlem('Frame');

for ii = 1:numel(hF)
    if isprop(hF(ii), 'FaceColor')
        set(hF(ii), 'FaceColor', 'none')
    end
    if isprop(hF(ii), 'EdgeColor')
        set(hF(ii), 'EdgeColor', 'k')
    end
    if isprop(hF(ii), 'Color')
        set(hF(ii), 'Color', 'k')
    end
    if isprop(hF(ii), 'LineWidth')
        set(hF(ii), 'LineWidth', 2.0)
    end
end

uistack(hF, 'top')

%% ===================== PANEL (c): MARION COUNTY =====================
subplot(1,3,3)
worldmap([min(lat1) max(lat1)], [min(lon1) max(lon1)]);

setm(gca,'FFaceColor',bgGrey)
hold on

surfm(lat_mesh, lon_mesh', T_total2_log)

% roads: subtle grey underlay + white top line
if haveRoads
plotm(LAT_road, LON_road, '.', ...
      'Color',[0.35 0.35 0.35], ...
      'MarkerSize', roadMS_outer)

plotm(LAT_road, LON_road, '.', ...
      'Color',[1 1 1], ...
      'MarkerSize', roadMS_inner)
end

gridm('off')
mlabel off
plabel off
tightmap

colormap(gca, cmap)
clim(log10([zmin1 zmax1]))

title({'Marion County, Indiana,'; ...
       'CH_4 1-km Emissions (mol s^{-1})'}, ...
       'FontSize', titleFS, ...
       'FontWeight', 'bold', ...
       'FontName', fontName)

text(labelX, labelY, '(c)', ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'Color','w', ...
    'FontSize',23, ...
    'FontWeight','bold')

framem on
setm(gca, 'Frame', 'on', 'FEdgeColor', 'k', 'FLineWidth', frameLW)

cb3 = colorbar;
cb3.FontSize = tickFS;
cb3.TickDirection = 'out';
cb3.LineWidth = 0.8;
cb3.Ticks = log10([0.05 0.1 0.3 1 2 3]);
cb3.TickLabels = {'0.05','0.1','0.3','1','2','3'};

%framem on
%setm(gca, 'Frame', 'on', 'FEdgeColor', 'k', 'FLineWidth', 2.0)

drawnow
hF = handlem('Frame');

for ii = 1:numel(hF)
    if isprop(hF(ii), 'FaceColor')
        set(hF(ii), 'FaceColor', 'none')
    end
    if isprop(hF(ii), 'EdgeColor')
        set(hF(ii), 'EdgeColor', 'k')
    end
    if isprop(hF(ii), 'Color')
        set(hF(ii), 'Color', 'k')
    end
    if isprop(hF(ii), 'LineWidth')
        set(hF(ii), 'LineWidth', 2.0)
    end
end

uistack(hF, 'top')

%% ===================== FINAL POLISH =====================
axs = findall(gcf,'type','axes');
set(axs, ...
    'FontName', fontName, ...
    'FontSize', tickFS, ...
    'LineWidth', lineW, ...
    'Color', 'w')

set(gcf,'Color','w')

%% ===================== EXPORT =====================
exportgraphics(gcf, 'INFLUX_emissions_log10_turbo.png', 'Resolution', 800)
%exportgraphics(gcf, 'INFLUX_emissions_log10_turbo.png', 'Resolution', 600)

% If exportgraphics does not work on your MATLAB version:
% print(gcf, '-dpng', '-r600', 'INFLUX_emissions_log10_turbo.png')
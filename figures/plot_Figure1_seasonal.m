%% plot_Figure1_seasonal.m
%  Revised Figure 1 (Reviewer 1, Figure 1 comment; Reviewer 2, Line 173)
%    (a) number of flights in each meteorological season, stacked by year
%    (b) wind rose of the mean (mass balance) wind direction on each flight
%        day, stacked by mean wind speed
%
%  Reads Indy_Mass_Balance_Table_withWD.xlsx (sheet FlightMasterSheet):
%    column A  date rows 'yyyy/mm/dd'
%    column E  Mean Flight Wind Speed (m/s)
%    column P  Mean Flight Wind Direction (deg, direction the wind blows FROM)
%  The wind-direction column is found by its header, so it can move.
%
%  The wind rose is drawn with patches on ordinary axes (no polaraxes), so the
%  speed classes can be stacked and the look is the same in any MATLAB version.
%  Base MATLAB only (R2019a+ for readcell).
%
%  Expected: seasons DJF/MAM/JJA/SON = 5/8/13/11 flights;
%            wind speed <5 / 5-10 / >10 m/s = 13/16/8 days.

clear; clc;
xlsFile = 'Indy_Mass_Balance_Table_withWD.xlsx';
outName = 'Fig1_flights_season_windrose';

%% ------------------------------------------------------------ read
C   = readcell(xlsFile, 'Sheet', 'FlightMasterSheet');
hdr = C(1,:);
isTxt = @(h) ischar(h) || isstring(h);
wdCol = find(cellfun(@(h) isTxt(h) && ~isempty(regexpi(char(h), 'wind\s*direction', 'once')), hdr), 1);
wsCol = find(cellfun(@(h) isTxt(h) && ~isempty(regexpi(char(h), 'wind\s*speed',     'once')), hdr), 1);
assert(~isempty(wdCol), 'No "Wind Direction" column found in %s', xlsFile);
assert(~isempty(wsCol), 'No "Wind Speed" column found in %s', xlsFile);

dates = datetime.empty(0,1); wdir = []; wspd = [];
for r = 2:size(C,1)
    lab = C{r,1};
    if isa(lab,'datetime')
        dt = lab;
    elseif isTxt(lab) && ~isempty(regexp(strtrim(char(lab)), '^\d{4}/\d{2}/\d{2}$', 'once'))
        dt = datetime(strtrim(char(lab)), 'InputFormat', 'yyyy/MM/dd');
    else
        continue
    end
    dates(end+1,1) = dt;                     %#ok<SAGROW>
    wdir(end+1,1)  = toNum(C{r, wdCol});     %#ok<SAGROW>
    wspd(end+1,1)  = toNum(C{r, wsCol});     %#ok<SAGROW>
end
nF = numel(dates);
fprintf('Found %d flight days (%s to %s); %d with wind direction\n', nF, ...
    datestr(min(dates),'yyyy-mm-dd'), datestr(max(dates),'yyyy-mm-dd'), nnz(~isnan(wdir)));

%% ------------------------------------------------------------ season x year
mo = month(dates); yr = year(dates);
sIdx = 4*ones(nF,1);                                   % 4 = SON
sIdx(ismember(mo,[12 1 2])) = 1; sIdx(ismember(mo,3:5)) = 2; sIdx(ismember(mo,6:8)) = 3;
years = min(yr):max(yr);
counts = zeros(4, numel(years));
for s = 1:4
    for y = 1:numel(years), counts(s,y) = sum(sIdx == s & yr == years(y)); end
end
seasonNames = {'Winter (DJF)','Spring (MAM)','Summer (JJA)','Fall (SON)'};
fprintf('Flights per season (DJF MAM JJA SON): %s\n', mat2str(sum(counts,2)'));

%% ------------------------------------------------------------ figure
fig = figure('Color','w','Units','pixels','Position',[100 100 1250 470]);

% (a) stacked bars, one colour per year (categorical, colour-blind safe)
yearCol = [0.00 0.45 0.70; 0.90 0.62 0.00; 0.00 0.62 0.45; 0.80 0.47 0.65; ...
           0.34 0.71 0.91; 0.84 0.37 0.00; 0.60 0.60 0.60; 0.30 0.30 0.30];
ax1 = subplot(1,2,1); hold(ax1,'on'); box(ax1,'on');
hb = bar(ax1, counts, 'stacked', 'BarWidth', 0.65);
for y = 1:numel(years)
    hb(y).FaceColor = yearCol(y,:); hb(y).EdgeColor = 'k';
end
tot = sum(counts,2);
for s = 1:4
    text(ax1, s, tot(s)+0.35, sprintf('n = %d', tot(s)), 'HorizontalAlignment','center', 'FontSize',11);
end
set(ax1, 'XTick',1:4, 'XTickLabel',seasonNames, 'FontSize',11, 'TickDir','out');
xlim(ax1, [0.4 4.6]); ylim(ax1, [0 max(tot)+2]);
ylabel(ax1, 'Number of flights');
legend(ax1, hb, arrayfun(@num2str, years, 'UniformOutput',false), ...
       'Location','northwest', 'NumColumns',2, 'Box','on');
title(ax1, sprintf('(a) Flights per season (%d flight days)', nF), 'FontWeight','normal');

% (b) wind rose, 16 sectors of 22.5 deg centred on N, NNE, ..., stacked by speed
ax2 = subplot(1,2,2); hold(ax2,'on'); axis(ax2,'equal'); axis(ax2,'off');
binW    = 22.5;
centers = 0:binW:360-binW;
sector  = mod(round(wdir/binW), numel(centers)) + 1;          % nearest centre
spdEdge = [0 5 10 Inf];
spdLab  = {'< 5 m s^{-1}', '5-10 m s^{-1}', '> 10 m s^{-1}'};
spdCol  = [0.78 0.86 0.94; 0.42 0.68 0.84; 0.03 0.32 0.61];
N = zeros(numel(centers), numel(spdEdge)-1);
for c = 1:numel(spdEdge)-1
    m = wspd >= spdEdge(c) & wspd < spdEdge(c+1) & ~isnan(wdir);
    N(:,c) = accumarray(sector(m), 1, [numel(centers) 1]);
end
rMax = max(sum(N,2));  rMax = 2*ceil((rMax+0.5)/2);           % even ring spacing
% rings and spokes
th = linspace(0, 2*pi, 361);
for rr = 2:2:rMax
    plot(ax2, rr*sin(th), rr*cos(th), '-', 'Color',[0.82 0.82 0.82], 'LineWidth',0.6);
    text(ax2, rr*sind(112.5), rr*cosd(112.5), sprintf('%d', rr), 'FontSize',9, ...
         'Color',[0.35 0.35 0.35], 'HorizontalAlignment','left', 'VerticalAlignment','top');
end
plot(ax2, rMax*sin(th), rMax*cos(th), 'k-', 'LineWidth',0.9);
for a = 0:45:315
    plot(ax2, [0 rMax*sind(a)], [0 rMax*cosd(a)], '-', 'Color',[0.82 0.82 0.82], 'LineWidth',0.6);
end
% sectors (direction FROM which the wind blows; N up, clockwise)
hS = gobjects(1, size(N,2));
for k = 1:numel(centers)
    r0 = 0;
    for c = 1:size(N,2)
        if N(k,c) == 0, continue; end
        r1 = r0 + N(k,c);
        a  = linspace(centers(k)-binW/2, centers(k)+binW/2, 20);
        xs = [r0*sind(a), fliplr(r1*sind(a))];
        ys = [r0*cosd(a), fliplr(r1*cosd(a))];
        hS(c) = patch(ax2, xs, ys, spdCol(c,:), 'EdgeColor','k', 'LineWidth',0.6);
        r0 = r1;
    end
end
compass = {'N','NE','E','SE','S','SW','W','NW'};
for i = 1:8
    a = (i-1)*45;
    text(ax2, 1.12*rMax*sind(a), 1.12*rMax*cosd(a), compass{i}, 'FontSize',11, ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
xlim(ax2, 1.22*rMax*[-1 1]); ylim(ax2, 1.22*rMax*[-1 1]);
ok = isgraphics(hS);
legend(ax2, hS(ok), spdLab(ok), 'Location','eastoutside', 'Box','on', 'FontSize',10);
title(ax2, '(b) Mean wind direction on flight days (number of days)', 'FontWeight','normal', 'FontSize',11);
fprintf('Wind speed classes (<5, 5-10, >10 m/s): %s days\n', mat2str(sum(N,1)));

%% ------------------------------------------------------------ save
set(fig, 'PaperPositionMode','auto');
print(fig, [outName '.png'], '-dpng', '-r300');
print(fig, [outName '.eps'], '-depsc');
fprintf('Saved %s.png / .eps\n', outName);

%% ------------------------------------------------------------ helper
function v = toNum(x)
    if isnumeric(x) && isscalar(x), v = double(x);
    elseif ischar(x) || isstring(x), v = str2double(strtrim(char(x)));
    else, v = NaN;
    end
end

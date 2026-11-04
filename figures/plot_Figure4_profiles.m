%% plot_Figure4_profiles.m
%  Revised Figure 4: potential temperature and CH4 profiles on 10/04/2012,
%  (a) the 17:00 UTC profile near the upwind transect, (b) the 19:05 UTC
%  profile near the last downwind transect.
%  Reviewer 1 (R1.8): both panels now share IDENTICAL axis limits
%  (CH4, potential temperature and height), so the ABL growth is read directly.
%
%  Profiles are binned exactly as in the Method scripts (50 m bins of height
%  above ground, same terrain lookup); the ABL heights are the ones the Method
%  scripts use (ManualABLHeights if set for this day, otherwise the same
%  theta/CH4 gradient test).

clear; close all; clc;

%% ------------------------------------------------------------ settings
Date    = '20121004';
MB_ROOT = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';
addpath(MB_ROOT); addpath(fullfile(MB_ROOT,'LegAndProfileList'));
addpath('/storage/work/zrb5027/WRFGUI_V3/Observations/INFLUX');
addpath('/storage/work/zrb5027/WRFGUI_V3/Code');            % ElevationData.nc lives with the Method-script helpers
ELEV = which('ElevationData.nc');
assert(~isempty(ELEV), 'ElevationData.nc not found: add its folder with addpath');
outName = sprintf('Fig4_profiles_%s', Date);
Z_TOP   = 2500;                 % top of the height axis (m), as in the submitted figure

%% ------------------------------------------------------------ data (as in the Method scripts)
L = load(fullfile(MB_ROOT,'LegAndProfileList',Date));
ProfileList = L.ProfileList;
x = load(['INFLUX_',Date]); x(x<-998) = NaN;
Time=x(:,1); Lat=x(:,2); Lon=x(:,3); Alt=x(:,4); P=x(:,9); T=x(:,10);
CH4 = movmean(x(:,6),30,'omitnan');
if isfield(L,'CH4Flagged'), CH4 = L.CH4Flagged; end
if median(CH4,'omitnan') < 100, CH4 = CH4*1000; end      % ppm -> ppb

Elev = ncread(ELEV,'data'); Elev(Elev<0) = 0;
WLat = ncread(ELEV,'lat');
WLon = ncread(ELEV,'lon'); WLon(WLon>180) = WLon(WLon>180) - 360;
Elevation = NaN(size(Alt));
for i = 1:length(Alt)
    Row = find(abs(WLon-Lon(i))<0.2, 1); Col = find(abs(WLat-Lat(i))<0.2, 1);
    if ~isempty(Row) && ~isempty(Col), Elevation(i) = Elev(Row,Col); end
end
AltAGL = Alt - Elevation;

%% ------------------------------------------------------------ profiles
nP = size(ProfileList,1);
Hedge = -25:50:5025;  Heights = 0:50:5000;
pr = struct('t',{},'th',{},'ch4',{},'abl',{});
for p = 1:nP
    idx = ProfileList(p,1):ProfileList(p,2);
    z = AltAGL(idx); th = T(idx).*(1000./P(idx)).^(2/7); ch = CH4(idx);
    thA = NaN(1,numel(Heights)); chA = thA;
    for h = 1:numel(Heights)
        g = z>=Hedge(h) & z<Hedge(h+1);
        thA(h) = mean(th(g),'omitnan'); chA(h) = mean(ch(g),'omitnan');
    end
    if isfield(L,'ManualABLHeights')
        abl = L.ManualABLHeights(p);
    else
        dTh = (thA(1:end-2)-thA(3:end))/100; dCH = (chA(1:end-2)-chA(3:end))/100/1000;   % ppm per m, as in the scripts
        k = find(dCH>2e-5 & dTh<-5e-3, 1); abl = Heights(k+1);
    end
    pr(p) = struct('t',mean(Time(idx),'omitnan'),'th',thA,'ch4',chA,'abl',abl);
    fprintf('Profile %d: %s UTC, ABL height %g m\n', p, hhmm(pr(p).t), abl);
end
assert(nP >= 2, 'Expected two profiles for %s', Date);
pU = 1; pD = nP;                          % first = upwind (17 UTC), last = downwind (19 UTC)

%% ------------------------------------------------------------ common limits
inZ = Heights <= Z_TOP;
allCH = [pr(pU).ch4(inZ) pr(pD).ch4(inZ)]; allTH = [pr(pU).th(inZ) pr(pD).th(inZ)];
chLim = [floor(min(allCH)/10)*10, ceil(max(allCH)/10)*10];
thLim = [floor(min(allTH)), ceil(max(allTH))];
fprintf('Common limits: CH4 %g-%g ppb, theta %g-%g K, height 0-%g m\n', chLim, thLim, Z_TOP);

%% ------------------------------------------------------------ figure
fig = figure('Color','w','Units','pixels','Position',[80 80 1150 470]);
tags = {'a.)','b.)'}; who = {'Upwind','Downwind'}; ps = [pU pD];
for q = 1:2
    p = ps(q);
    pos = [0.07+0.49*(q-1) 0.13 0.40 0.66];
    ax1 = axes('Parent',fig,'Position',pos); hold(ax1,'on');
    plot(ax1, pr(p).ch4, Heights, '-k', 'LineWidth', 1.2);
    set(ax1,'XLim',chLim,'YLim',[0 Z_TOP],'FontSize',12,'TickDir','out','Box','off');
    xlabel(ax1,'Methane Mole Fraction (ppb)'); ylabel(ax1,'Height AGL (m)');
    ax2 = axes('Parent',fig,'Position',pos,'Color','none','XAxisLocation','top', ...
        'YAxisLocation','right','YTick',[],'XColor','r','FontSize',12,'TickDir','out','Box','off');
    hold(ax2,'on');
    plot(ax2, pr(p).th, Heights, '-r', 'LineWidth', 1.2);
    set(ax2,'XLim',thLim,'YLim',[0 Z_TOP]);
    xlabel(ax2,'Potential Temperature (K)');
    % ABL height of this profile
    plot(ax1, chLim, pr(p).abl*[1 1], '-k', 'LineWidth', 3);
    text(ax1, mean(chLim), pr(p).abl, sprintf('ABL Height = %g m', pr(p).abl), ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',11,'BackgroundColor','w');
    if q == 1     % downwind ABL height shown on the upwind panel (Adjusted ABL Upwind method)
        plot(ax1, chLim, pr(pD).abl*[1 1], '--k', 'LineWidth', 1.5);
        text(ax1, mean(chLim), pr(pD).abl, sprintf('Adjusted ABL Height = %g m', pr(pD).abl), ...
            'HorizontalAlignment','center','VerticalAlignment','top','FontSize',11,'BackgroundColor','w');
    end
    text(ax1, 0.06, 0.95, tags{q}, 'Units','normalized','FontSize',14,'FontWeight','bold');
    title(ax2, sprintf('%s/%s/%s: %s UTC (%s)', Date(5:6), Date(7:8), Date(1:4), hhmm(pr(p).t), who{q}), ...
        'FontSize',14,'FontWeight','bold','Color','k');
end

set(fig,'PaperPositionMode','auto');
print(fig, [outName '.png'], '-dpng', '-r300');
print(fig, [outName '.eps'], '-depsc');
fprintf('Saved %s.png / .eps\n', outName);

function s = hhmm(t)
    s = sprintf('%02d:%02d', floor(t/3600), round(mod(t,3600)/60));
end

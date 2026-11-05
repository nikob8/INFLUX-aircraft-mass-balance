%% plot_FigS1_Fig10_trends.m
%  Makes two figures from Zach's TempAndTime*Plots.m scripts, with the fluxes
%  recomputed with cos(theta) and the p-value of each slope printed (R1.10):
%     Figure S1  : citywide Edge emission rate vs. flight date and vs. takeoff
%                  temperature (was Figure 9 in the submission)
%     Figure 10  : SSLF emission rate vs. flight date and vs. temperature
%                  (was Figure 11 in the submission)
%
%  Inputs
%     TempsAndTimes.mat in INFLUXMassBalance (one row per flight day)
%        TimeMatrix(:,1) old flux, TimeMatrix(:,2) day number
%        TempMatrix(:,1) old flux, TempMatrix(:,2) takeoff temperature (deg F)
%        Only the day numbers and temperatures are used; the fluxes are replaced.
%        The takeoff temperature of a flight day is the same for the city and
%        the landfill, so this one file serves both figures (the older copy in
%        INFLUXMassBalance/Landfill reflects an earlier plume selection).
%     All_Solutions.mat (cos theta Edge fluxes), Landfill/<Date>_Loop<k>.mat
%
%  The rows of TimeMatrix are matched to flight dates by their day spacing, and
%  the script stops if the match is not exact.

clear; close all; clc;

%% ------------------------------------------------------------ settings
MB_ROOT  = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';
LF_DIR   = fullfile(pwd, 'Landfill');                    % Loop files written by MassBalanceLandfill.m
TT_FILE  = fullfile(MB_ROOT, 'TempsAndTimes.mat');       % takeoff temperatures and day numbers of all flight days
SOL_FILE = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB/All_Solutions.mat';
FN = 'Arial';

%% ------------------------------------------------------------ cos(theta) fluxes
S = load(SOL_FILE, 'SolA'); SolA = S.SolA;
cityDates = {SolA.Date};
cityFlux  = arrayfun(@(s) mean(s.Flux), SolA);                       % daily Edge means, cos theta

lfDates = {}; lfFlux = [];
first = dir(fullfile(LF_DIR, '*Loop1*'));
dates = sort(unique(cellfun(@(n) n(1:8), {first.name}, 'UniformOutput', false)));
for i = 1:numel(dates)
    D = dates{i}; files = dir(fullfile(LF_DIR, [D '_Loop*.mat'])); fl = [];
    for j = 1:numel(files)
        L = load(fullfile(files(j).folder, files(j).name));
        if L.Flag > 0
            k  = sscanf(regexprep(files(j).name, '.*_Loop', ''), '%d');
            kk = find(strcmp(cityDates, D), 1);
            assert(~isempty(kk) && numel(SolA(kk).CosTheta) >= k, 'No CosTheta for %s transect %d', D, k);
            fl(end+1) = L.Flux * SolA(kk).CosTheta(k); %#ok<AGROW>
        end
    end
    if ~isempty(fl), lfDates{end+1} = D; lfFlux(end+1) = mean(fl); end %#ok<AGROW>
end
fprintf('%d city days, %d landfill days with cos(theta) fluxes\n', numel(cityFlux), numel(lfFlux));

%% ------------------------------------------------------------ the two figures
origin = make_trend_figure(TT_FILE, cityDates, cityFlux, ...
    'Citywide CH_4 emission rate (mol s^{-1})', 'FigS1_city_trends', FN, []);
make_trend_figure(TT_FILE, lfDates, lfFlux, ...
    'SSLF CH_4 emission rate (mol s^{-1})', 'Fig10_SSLF_trends', FN, origin);


%% ============================================================ one figure
function origin = make_trend_figure(matFile, dateList, newFlux, ylab, outName, FN, origin)
% origin = datenum of day number 0 in TimeMatrix(:,2); [] = derive it here
sameKind = isempty(origin);          % the old fluxes in the file are the city fluxes: compare only in the city call
M = load(matFile);
TimeMatrix = M.TimeMatrix; TempMatrix = M.TempMatrix;
n = size(TimeMatrix, 1);
% rows of TempMatrix follow TimeMatrix (same old fluxes in column 1)
[~, iT] = ismember(round(TempMatrix(:,1), 6), round(TimeMatrix(:,1), 6));
assert(all(iT > 0), 'TempMatrix rows could not be matched to TimeMatrix rows');
tempF = NaN(n,1); tempF(iT) = TempMatrix(:,2);
tempC = (tempF - 32)*5/9;
dayNum = TimeMatrix(:,2); oldFlux = TimeMatrix(:,1);
dn = datenum(dateList, 'yyyymmdd'); dn = dn(:); newFlux = newFlux(:);

if isempty(origin)                     % derive from this file: rows and dates both chronological
    assert(n == numel(newFlux), '%s has %d rows but %d flight days have fluxes', matFile, n, numel(newFlux));
    [dayNum, ro] = sort(dayNum); tempC = tempC(ro); oldFlux = oldFlux(ro);
    [dn, od] = sort(dn); newFlux = newFlux(od);
    dd = max(abs(diff(dayNum) - diff(dn)));
    assert(dd <= 1, 'Day spacing in %s does not match the flight dates (max difference %g days)', matFile, dd);
    origin = round(mean(dn - dayNum));
    fprintf('\n%s: day number 0 = %s\n', outName, datestr(origin));
    rowDate = dayNum + origin;
else                                   % convert day numbers with the known origin, match by date
    rowDate = dayNum + origin;
    keep = false(n,1); matched = NaN(n,1);
    for i = 1:n                        % exact date match: the rows are the flight days themselves
        [dmin, j] = min(abs(dn - rowDate(i)));
        if dmin <= 0.5, keep(i) = true; matched(i) = j; end
    end
    assert(numel(unique(matched(keep))) == nnz(keep), 'a flight day was matched to more than one row of %s', matFile);
    fprintf('\n%s: %d of %d rows matched to the %d flight days with fluxes\n', outName, nnz(keep), n, numel(dn));
    if any(~keep)
        fprintf('  rows with no matching flight day (dropped): %s\n', strjoin(cellstr(datestr(rowDate(~keep), 'yyyy-mm-dd'))', ' '));
    end
    missing = setdiff(1:numel(dn), matched(keep));
    if ~isempty(missing)
        fprintf('  flight days with a flux but no row in the file (dropped): %s\n', strjoin(cellstr(datestr(dn(missing), 'yyyy-mm-dd'))', ' '));
    end
    rowDate = rowDate(keep); tempC = tempC(keep); oldFlux = oldFlux(keep);
    dn = dn(matched(keep)); newFlux = newFlux(matched(keep)); n = nnz(keep);
end
if sameKind, fprintf('  mean |old - new| flux = %.1f mol/s (n = %d)\n', mean(abs(oldFlux - newFlux)), n); end

x1 = datetime(dn, 'ConvertFrom', 'datenum'); y = newFlux(:);
yrs = years(x1 - x1(1));
a = fitlm(yrs, y);  sA = a.Coefficients.Estimate(2); eA = a.Coefficients.SE(2); pA = a.Coefficients.pValue(2);
ok = ~isnan(tempC);
b = fitlm(tempC(ok), y(ok)); sB = b.Coefficients.Estimate(2); eB = b.Coefficients.SE(2); pB = b.Coefficients.pValue(2);
fprintf('  vs time:        slope %.2f ± %.2f mol/s per year, p = %.3f (n = %d)\n', sA, eA, pA, n);
fprintf('  vs temperature: slope %.2f ± %.2f mol/s per degC, p = %.3f (n = %d)\n', sB, eB, pB, nnz(ok));

fig = figure('Color','w','Units','pixels','Position',[100 100 1150 450]);
% (a) time
ax1 = subplot(1,2,1); hold(ax1,'on'); box(ax1,'on');
scatter(ax1, x1, y, 55, [0.85 0.15 0.15], 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.4);
xf = linspace(min(x1), max(x1), 50);
plot(ax1, xf, a.Coefficients.Estimate(1) + sA*years(xf - x1(1)), '-k', 'LineWidth', 2);
xlabel(ax1, 'Flight date'); ylabel(ax1, ylab);
text(ax1, 0.03, 0.96, sprintf('slope = %.1f \\pm %.1f mol s^{-1} yr^{-1}\np = %.2f', sA, eA, pA), ...
     'Units','normalized', 'VerticalAlignment','top', 'FontSize', 11, 'FontName', FN, 'BackgroundColor','w', 'EdgeColor','k', 'Margin', 3);
text(ax1, 0.03, 0.06, '(a)', 'Units','normalized', 'FontSize', 13, 'FontWeight','bold');
% (b) temperature
ax2 = subplot(1,2,2); hold(ax2,'on'); box(ax2,'on');
scatter(ax2, tempC(ok), y(ok), 55, [0.1 0.35 0.8], 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.4);
tf = linspace(min(tempC(ok)), max(tempC(ok)), 50);
plot(ax2, tf, b.Coefficients.Estimate(1) + sB*tf, '-k', 'LineWidth', 2);
xlabel(ax2, 'Temperature at takeoff (\circC)'); ylabel(ax2, ylab);
text(ax2, 0.03, 0.96, sprintf('slope = %.2f \\pm %.2f mol s^{-1} \\circC^{-1}\np = %.2f', sB, eB, pB), ...
     'Units','normalized', 'VerticalAlignment','top', 'FontSize', 11, 'FontName', FN, 'BackgroundColor','w', 'EdgeColor','k', 'Margin', 3);
text(ax2, 0.03, 0.06, '(b)', 'Units','normalized', 'FontSize', 13, 'FontWeight','bold');
yl = [min(0, floor(min(y)/10)*10 - 10), ceil(max(y)/10)*10 + 10];
for ax = [ax1 ax2], ylim(ax, yl); set(ax, 'FontSize', 12, 'FontName', FN, 'TickDir','out'); grid(ax,'on'); ax.GridAlpha = 0.15; end
set(fig, 'PaperPositionMode', 'auto');
print(fig, [outName '.png'], '-dpng', '-r300');
print(fig, [outName '.eps'], '-depsc');
fprintf('  saved %s.png / .eps\n', outName);
end

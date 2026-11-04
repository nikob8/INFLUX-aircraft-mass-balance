%% plot_Figure13_diagnostic.m
%  Revised Figure 13 (Figure 14 in the submission): diagnostic ratio vs. the
%  Edge emission rate for the 37 flight days, coloured by mean wind speed.
%  Reviewer 1 (R1.13): the region below the threshold ratio is shaded and the
%  mean and standard deviation of the two subsets are drawn (bands on either
%  side of the threshold) and printed; the overall mean is a dashed line.
%
%  Ratios and wind speeds: Diag4.mat (Diag_All2: wind speed, wind direction,
%  diagnostic ratio, Edge flux of the submission, yyyymmdd).
%  Emission rates: daily means of the Edge transects WITH cos(theta), read from
%  All_Solutions.mat, so the figure matches Table S1 v3.  The script prints the
%  numbers for Table 3 (all / good / good with wind filter / bad).

clear; close all; clc;

%% ------------------------------------------------------------ settings
DIAG_FILE = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/WRF/Diag_ratio.mat';   % written by analyses/Diagnostic_ratio.m (same layout as the old Diag4.mat)
SOL_FILE  = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB/All_Solutions.mat';
THRESH    = 2;          % diagnostic ratio threshold
WIND_MIN  = 3;          % m/s, wind filter used in Table 3
outName   = 'Fig13_diagnostic_ratio';
FN = 'Arial';

%% ------------------------------------------------------------ data
assert(isfile(DIAG_FILE), 'Diag4.mat not found at %s (fix DIAG_FILE)', DIAG_FILE);
assert(isfile(SOL_FILE),  'All_Solutions.mat not found at %s (fix SOL_FILE)', SOL_FILE);
D = load(DIAG_FILE); D = D.Diag_All2;
ws = D(:,1); ratio = D(:,3); fluxOld = D(:,4);
dates = arrayfun(@(v) sprintf('%d', v), D(:,5), 'UniformOutput', false);
S = load(SOL_FILE, 'SolA'); SolA = S.SolA;
flux = NaN(size(ratio));
for i = 1:numel(dates)
    k = find(strcmp({SolA.Date}, dates{i}), 1);
    assert(~isempty(k), 'No Solution-MethodA for %s in All_Solutions.mat', dates{i});
    flux(i) = mean(SolA(k).Flux);                       % daily mean of the Edge transects, with cos(theta)
end
fprintf('Daily Edge means with cos(theta); largest change from the submitted values: %.1f mol/s\n', max(abs(flux - fluxOld)));

good = ratio > THRESH;  bad = ~good;  gw = good & ws >= WIND_MIN;
st = @(v) [numel(v) mean(v) std(v) std(v)/sqrt(numel(v))];
T  = [st(flux); st(flux(good)); st(flux(gw)); st(flux(bad))];
rowName = {'All flights', 'Ratio > 2 (good)', 'Ratio > 2 and wind >= 3 m/s', 'Ratio <= 2 (bad)'};
fprintf('\nTable 3 with cos(theta):\n');
for r = 1:4, fprintf('  %-30s n = %2d   mean %5.1f   SD %5.1f   SE %4.1f\n', rowName{r}, T(r,:)); end
fprintf('Good days: %s\n', strjoin(dates(good)', ' '));
fprintf('Good days removed by the wind filter: %s\n', strjoin(dates(good & ws < WIND_MIN)', ' '));

%% ------------------------------------------------------------ figure
fig = figure('Color','w','Position',[100 100 860 600]);
ax = axes(fig); hold(ax,'on'); box(ax,'on');
xl = [0 ceil(max(ratio)) + 0.5];
yl = [floor((min(flux) - 10)/10)*10, ceil((max(flux) + 10)/10)*10];

% region below the threshold
patch(ax, [xl(1) THRESH THRESH xl(1)], [yl(1) yl(1) yl(2) yl(2)], [0.91 0.91 0.91], 'EdgeColor','none');

% mean and +-1 SD of each subset, drawn across its side of the threshold
cBand = [0.10 0.35 0.75];
for r = [4 2]                                   % 4 = ratio <= 2, 2 = ratio > 2
    if r == 4, xr = [xl(1) THRESH]; else, xr = [THRESH xl(2)]; end
    m = T(r,2); s = T(r,3);
    hBand = patch(ax, [xr(1) xr(2) xr(2) xr(1)], [m-s m-s m+s m+s], cBand, 'FaceAlpha', 0.12, 'EdgeColor','none');
    plot(ax, xr, [m m], '-', 'Color', cBand, 'LineWidth', 2);
    plot(ax, xr, [m-s m-s], ':', 'Color', cBand, 'LineWidth', 1.2);
    plot(ax, xr, [m+s m+s], ':', 'Color', cBand, 'LineWidth', 1.2);
end
hAll = plot(ax, xl, T(1,2)*[1 1], '--k', 'LineWidth', 1.3);          % mean of all flight days
xline(ax, THRESH, '-k', 'LineWidth', 1.5);

% the flight days
hPts = scatter(ax, ratio, flux, 80, ws, 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.3);
colormap(ax, jet); caxis(ax, [min(ws) 12]);
cb = colorbar(ax); cb.Label.String = 'Mean wind speed (m s^{-1})'; cb.Label.FontSize = 13; cb.FontSize = 11;

xlim(ax, xl); ylim(ax, yl); grid(ax, 'on'); ax.GridAlpha = 0.15;
set(ax, 'FontSize', 13, 'LineWidth', 1.2, 'TickDir', 'out', 'FontName', FN, 'Layer', 'top');
xlabel(ax, 'Diagnostic ratio', 'FontSize', 15, 'FontName', FN);
ylabel(ax, 'Downwind Edge method emission rate (mol s^{-1})', 'FontSize', 15, 'FontName', FN);

% subset statistics
txtL = sprintf('Ratio \\leq %d  (n = %d)\nmean %.0f, SD %.0f mol s^{-1}', THRESH, T(4,1), T(4,2), T(4,3));
txtR = sprintf('Ratio > %d  (n = %d)\nmean %.0f, SD %.0f mol s^{-1}\nwind \\geq %d m s^{-1}: n = %d, mean %.0f, SD %.0f', ...
               THRESH, T(2,1), T(2,2), T(2,3), WIND_MIN, T(3,1), T(3,2), T(3,3));
text(ax, 0.02, 0.97, txtL, 'Units','normalized', 'FontSize', 11, 'FontName', FN, 'VerticalAlignment','top', ...
     'BackgroundColor','w', 'EdgeColor','k', 'Margin', 4);
text(ax, 0.98, 0.97, txtR, 'Units','normalized', 'FontSize', 11, 'FontName', FN, 'VerticalAlignment','top', ...
     'HorizontalAlignment','right', 'BackgroundColor','w', 'EdgeColor','k', 'Margin', 4);
legend(ax, [hPts hBand hAll], {'flight day', 'subset mean \pm 1 SD', 'mean of all flight days'}, ...
       'Location', 'southeast', 'FontSize', 11, 'FontName', FN);

set(fig, 'PaperPositionMode', 'auto');
print(fig, [outName '.png'], '-dpng', '-r300');
print(fig, [outName '.eps'], '-depsc');
fprintf('Saved %s.png / .eps\n', outName);

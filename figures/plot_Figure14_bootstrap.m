%% plot_Figure14_bootstrap.m
%  Figure 14: expected standard error of the mean citywide emission rate as a
%  function of the number of flights, for two sampling pools: all 37 flight
%  days and the 23 days with a diagnostic ratio above 2 (Section 3.2.2).
%
%  Method: for n = 2 ... NMAX flights, draw n daily Downwind Edge estimates
%  with replacement from the pool, compute their mean, repeat NBOOT times and
%  take the standard deviation of the NBOOT means as the standard error of an
%  n-flight campaign.  The dashed lines are the analytic sigma/sqrt(n) of each
%  pool.  The script prints the numbers quoted in the text (standard error for
%  5 flights from each pool, and the number of unfiltered flights that gives
%  the standard error of 5 filtered flights).
%
%  Inputs: All_Solutions.mat (daily means of the Downwind Edge transects, with
%  cos theta) and Diag4.mat (diagnostic ratio of each flight day).

clear; close all; clc;

%% ------------------------------------------------------------ settings
SOL_FILE  = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB/All_Solutions.mat';
DIAG_FILE = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/WRF/Diag_ratio.mat';   % written by analyses/Diagnostic_ratio.m (same layout as the old Diag4.mat)
THRESH = 2; NMAX = 20; NBOOT = 5000; SEED = 1;
outName = 'Fig14_bootstrap_SE'; FN = 'Arial';

%% ------------------------------------------------------------ daily estimates and the two pools
S = load(SOL_FILE, 'SolA'); SolA = S.SolA;
D = load(DIAG_FILE); D = D.Diag_All2;                 % columns: wind speed, wind direction, ratio, old flux, yyyymmdd
dates = arrayfun(@(v) sprintf('%d', v), D(:,5), 'UniformOutput', false);
flux = NaN(numel(dates),1);
for i = 1:numel(dates)
    k = find(strcmp({SolA.Date}, dates{i}), 1); assert(~isempty(k), 'no Solution for %s', dates{i});
    flux(i) = mean(SolA(k).Flux);
end
poolAll = flux; poolGood = flux(D(:,3) > THRESH);
fprintf('pools: all %d days (mean %.1f, SD %.1f), ratio > %g: %d days (mean %.1f, SD %.1f)\n', ...
    numel(poolAll), mean(poolAll), std(poolAll), THRESH, numel(poolGood), mean(poolGood), std(poolGood));

%% ------------------------------------------------------------ bootstrap
rng(SEED);
n = 2:NMAX; seAll = NaN(size(n)); seGood = NaN(size(n));
for j = 1:numel(n)
    mAll  = mean(poolAll(randi(numel(poolAll),  n(j), NBOOT)), 1);
    mGood = mean(poolGood(randi(numel(poolGood), n(j), NBOOT)), 1);
    seAll(j) = std(mAll); seGood(j) = std(mGood);
end
sigAll = std(poolAll, 1); sigGood = std(poolGood, 1);     % population SD of each pool (what the bootstrap converges to)
anAll = sigAll./sqrt(n); anGood = sigGood./sqrt(n);

% numbers for the text
i5 = find(n == 5);
[sa, ia] = unique(seAll); nEq = interp1(sa, n(ia), seGood(i5));   % unfiltered flights with the SE of 5 filtered flights
fprintf('5 flights: SE = %.1f mol/s (all) vs %.1f mol/s (ratio > %g)\n', seAll(i5), seGood(i5), THRESH);
fprintf('unfiltered flights needed for the SE of 5 filtered flights: %.1f\n', nEq);
fprintf('analytic sigma/sqrt(n) at n = 5: %.1f (all), %.1f (good); n for all to reach %.1f: %.1f\n', anAll(i5), anGood(i5), anGood(i5), (sigAll/anGood(i5))^2);

%% ------------------------------------------------------------ figure
cA = [0.20 0.45 0.75]; cG = [0.85 0.45 0.15];
fig = figure('Color','w','Position',[100 100 760 560]); ax = axes(fig); hold(ax,'on'); box(ax,'on');
p1 = plot(ax, n, seAll,  '-', 'LineWidth', 2.4, 'Color', cA);
p2 = plot(ax, n, seGood, '-', 'LineWidth', 2.4, 'Color', cG);
p3 = plot(ax, n, anAll,  '--', 'LineWidth', 1.0, 'Color', cA);
plot(ax, n, anGood, '--', 'LineWidth', 1.0, 'Color', cG);
% guides: 5 flights from each pool, and the unfiltered sample size with the same SE as 5 filtered flights
plot(ax, [5 5], [0 seAll(i5)], ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
plot(ax, [2 nEq], [seGood(i5) seGood(i5)], ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
plot(ax, [nEq nEq], [0 seGood(i5)], ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
text(ax, 5.3, seAll(i5) + 3.5, sprintf('%.0f vs %.0f mol s^{-1}\nfor 5 flights', seAll(i5), seGood(i5)), 'FontSize', 11, 'FontName', FN);
text(ax, nEq + 0.3, seGood(i5) + 2.6, sprintf('%.0f unfiltered flights give\nthe SE of 5 filtered flights', nEq), 'FontSize', 11, 'FontName', FN);
xlabel(ax, 'Number of flights', 'FontSize', 15, 'FontName', FN);
ylabel(ax, 'Standard error of the mean (mol s^{-1})', 'FontSize', 15, 'FontName', FN);
legend(ax, [p1 p2 p3], {sprintf('all %d flight days', numel(poolAll)), sprintf('%d days with ratio > %g', numel(poolGood), THRESH), '\sigma / n^{1/2} of each pool'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 12, 'FontName', FN);
xlim(ax, [2 NMAX]); xticks(ax, 2:2:NMAX); ylim(ax, [0 ceil(max(seAll)/5)*5]);
set(ax, 'FontSize', 13, 'FontName', FN, 'LineWidth', 1.1, 'TickDir', 'out', 'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.15);
set(fig, 'PaperPositionMode', 'auto');
print(fig, [outName '.png'], '-dpng', '-r300'); print(fig, [outName '.eps'], '-depsc');
fprintf('Saved %s.png / .eps\n', outName);

%% plot_Figure11_SSLF_comparison.m
%  Figure 11: estimates of the South Side Landfill (SSLF) emission rate (Section 3.1.3).
%
%  Inputs (current folder):
%    Landfill_daily.csv   (data/ folder of the repository; or Table S1,
%                         Indy_Mass_Balance_Table_v5.csv) daily SSLF emission rates of
%                         the 31 flight days with a quantified plume (Collect_Landfill.m);
%                         their mean is the 26 +/- 5 mol/s of the paper.
%    EPA_SSLF_GHGRP.csv   EPA Greenhouse Gas Reporting Program estimates for the SSLF,
%                         2011 to 2023, from Equations HH-6 and HH-8 (mol/s), with the
%                         method used for the reported value. If it is missing, EPA.mat
%                         (variables a = HH-6 and b = HH-8, 2011 to 2023, mol/s; the file of
%                         the original figure) is read instead and the csv is written from it.
%  Output: Fig11_SSLF_comparison.png (300 dpi) and .pdf; the annual means are printed.
%
%  Published values drawn in the figure:
%    Barkley et al. (2026), tower network observations 2013 to 2024: 1,476 kg/h from all
%        towers, 1,330 to 1,760 kg/h across towers and sensitivity tests.
%    Nesser et al. (2024), TROPOMI, 2019: 77 mol/s, ensemble range 64 to 103 mol/s.
%  Annual means are shown with 2 sigma = 2 standard errors of the daily values of that year.
%  Uses base MATLAB only (no toolboxes).

clear; close all; clc;

K = 16.04*3.6;                              % 1 mol/s of CH4 = 57.74 kg/h
YEARS_EPA = (2011:2023)';
REP_HH6 = [2011 2016];  REP_HH8 = [2016 2023];   % method of the reported GHGRP value (as in the original figure)
BARK = struct('mean', 1476/K, 'lo', 1330/K, 'hi', 1760/K, 't0', 2013, 't1', 2024);
NESS = struct('year', 2019, 'mean', 77, 'lo', 64, 'hi', 103);

C_STUDY = [0.122 0.306 0.612];  C_BARK = [0.482 0.196 0.580];  C_NESS = [0.878 0.510 0.078];
C_HH6 = [0.173 0.627 0.173];    C_HH8 = [0.839 0.153 0.157];
FS = 11;
[FONT, UNI] = pick_font();
set(0, 'DefaultAxesFontName', FONT, 'DefaultTextFontName', FONT);
if UNI
    INT = 'none';  sCH4 = ['CH' char(8324)];  sPerS = [' s' char(8315) char(185)];  sPerH = [' h' char(8315) char(185)];
else
    INT = 'tex';   sCH4 = 'CH_4';               sPerS = ' s^{-1}';                     sPerH = ' h^{-1}';
end
fprintf('Font used for all text: %s (%s)\n', FONT, ternary(UNI, 'Unicode sub/superscripts', 'TeX sub/superscripts'));

%% ------------------------------------------------------------ EPA values
if isfile('EPA_SSLF_GHGRP.csv')
    E = readtable('EPA_SSLF_GHGRP.csv');
    yrs = E.Year(:);  hh6 = E.HH6_mols(:);  hh8 = E.HH8_mols(:);
else
    assert(isfile('EPA.mat'), 'EPA_SSLF_GHGRP.csv (or EPA.mat) not found in this folder');
    S = load('EPA.mat');                     % a = HH-6, b = HH-8, 2011 to 2023 (mol/s)
    yrs = YEARS_EPA;  hh6 = double(S.a(:));  hh8 = double(S.b(:));
    assert(numel(hh6) == numel(yrs) && numel(hh8) == numel(yrs), 'EPA.mat: expected 13 values (2011 to 2023) in a and b');
    rep = repmat({'HH-8'}, numel(yrs), 1);  rep(yrs < REP_HH8(1)) = {'HH-6'};
    writetable(table(yrs, hh6, hh8, rep, 'VariableNames', {'Year', 'HH6_mols', 'HH8_mols', 'Reported_method'}), 'EPA_SSLF_GHGRP.csv');
    fprintf('wrote EPA_SSLF_GHGRP.csv from EPA.mat (add it to the data/ folder of the repository)\n');
end

%% ------------------------------------------------------------ this study
if isfile('Landfill_daily.csv')                          % one row per flight day
    opts = detectImportOptions('Landfill_daily.csv');
    opts = setvartype(opts, 'Date', 'char');
    T = readtable('Landfill_daily.csv', opts);
    dates = T.Date;  f = T.Flux_Landfill_mean_mols;
elseif isfile('Indy_Mass_Balance_Table_v5.csv')         % Table S1, one row per transect: daily means of the landfill column
    opts = detectImportOptions('Indy_Mass_Balance_Table_v5.csv');
    opts = setvartype(opts, 'Date', 'char');
    T = readtable('Indy_Mass_Balance_Table_v5.csv', opts);
    [dates, ~, g] = unique(T.Date);
    f = accumarray(g, T.Flux_Landfill_mols, [], @(v) mean(v, 'omitnan'));
else
    error('Landfill_daily.csv (or Indy_Mass_Balance_Table_v5.csv) not found: copy it from the data/ folder of the repository.');
end
yr = cellfun(@(d) str2double(d(1:4)), dates);             % dates are yyyy-mm-dd
ok = ~isnan(f);
fprintf('all days: n = %d, mean %.1f +/- %.1f mol/s\n', nnz(ok), mean(f(ok)), 2*std(f(ok))/sqrt(nnz(ok)));
Y = (min(yr):max(yr))';  m = NaN(size(Y));  se2 = NaN(size(Y));  n = zeros(size(Y));
for k = 1:numel(Y)
    v = f(yr == Y(k) & ok);  n(k) = numel(v);
    if n(k) > 0, m(k) = mean(v); end
    if n(k) > 1, se2(k) = 2*std(v)/sqrt(n(k)); end
    fprintf('%d: n = %d, mean %.1f, 2SE %.1f mol/s\n', Y(k), n(k), m(k), se2(k));
end
keepY = n > 0;  Y = Y(keepY);  m = m(keepY);  se2 = se2(keepY);  n = n(keepY);

%% ------------------------------------------------------------ figure
fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [80 80 900 560]);
ax = axes('Parent', fig, 'Position', [0.09 0.11 0.80 0.84]);
hold(ax, 'on');
% Barkley et al. (2026): band = range across towers and sensitivity tests, line = all-tower estimate
hBand = patch(ax, [BARK.t0 BARK.t1 BARK.t1 BARK.t0], [BARK.lo BARK.lo BARK.hi BARK.hi], C_BARK, 'EdgeColor', 'none', 'FaceAlpha', 0.16);
plot(ax, [BARK.t0 BARK.t1], [BARK.mean BARK.mean], '-', 'Color', C_BARK, 'LineWidth', 2.2);
hB = patch(ax, NaN(1,4), NaN(1,4), C_BARK, 'FaceAlpha', 0.16, 'EdgeColor', C_BARK, 'LineWidth', 2);   % legend entry (band + line)
% EPA GHGRP: black underlay = method used for the reported value
k6 = yrs >= REP_HH6(1) & yrs <= REP_HH6(2);  k8 = yrs >= REP_HH8(1) & yrs <= REP_HH8(2);
hR = plot(ax, yrs(k6), hh6(k6), '-', 'Color', [0 0 0], 'LineWidth', 7.5);
plot(ax, yrs(k8), hh8(k8), '-', 'Color', [0 0 0], 'LineWidth', 7.5);
h6 = plot(ax, yrs, hh6, '-', 'Color', C_HH6, 'LineWidth', 2.6);
h8 = plot(ax, yrs, hh8, '-', 'Color', C_HH8, 'LineWidth', 2.6);
% Nesser et al. (2024) and this study
hN = errorbar(ax, NESS.year, NESS.mean, NESS.mean - NESS.lo, NESS.hi - NESS.mean, 'o', 'Color', C_NESS, ...
    'MarkerFaceColor', C_NESS, 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 8, 'LineWidth', 2, 'CapSize', 8);
e2 = se2;  e2(isnan(e2)) = 0;
hS = errorbar(ax, Y, m, e2, 'o', 'Color', C_STUDY, 'MarkerFaceColor', C_STUDY, 'MarkerEdgeColor', [0 0 0], ...
    'MarkerSize', 8, 'LineWidth', 2, 'CapSize', 8);
for k = 1:numel(Y)
    text(ax, Y(k), 2.5, sprintf('n = %d', n(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', FS-2, 'Color', [0.35 0.35 0.35]);
end
uistack(hBand, 'bottom');
set(ax, 'XLim', [2010.4 2024.6], 'YLim', [0 110], 'XTick', 2011:2024, 'FontSize', FS, 'TickDir', 'out', 'Box', 'off', 'Layer', 'top');
xl = arrayfun(@(y) sprintf('%d', y), 2011:2024, 'UniformOutput', false);  xl(2:2:end) = {''};
set(ax, 'XTickLabel', xl, 'YGrid', 'on', 'GridAlpha', 0.15);
ylabel(ax, ['SSLF emission rate (mol' sPerS ')'], 'Interpreter', INT, 'FontSize', FS+1);
% right axis in kg CH4 per hour
ax2 = axes('Parent', fig, 'Position', get(ax, 'Position'), 'Color', 'none', 'YAxisLocation', 'right', ...
    'XAxisLocation', 'top', 'XTick', [], 'XColor', 'none', 'YLim', get(ax, 'YLim')*K, 'FontSize', FS, 'TickDir', 'out', 'Box', 'off');
ylabel(ax2, ['SSLF emission rate (kg ' sCH4 sPerH ')'], 'Interpreter', INT, 'FontSize', FS+1);
line(ax, [2010.4 2024.6], [110 110], 'Color', [0 0 0], 'LineWidth', 0.5);      % close the frame at the top
lg = legend(ax, [hS hB hN h6 h8 hR], {'This study (aircraft mass balance)', 'Barkley et al. (2026), tower network', ...
    'Nesser et al. (2024), TROPOMI', 'EPA GHGRP, Equation HH-6', 'EPA GHGRP, Equation HH-8', 'EPA reported value'}, ...
    'Location', 'northwest', 'FontSize', FS-1, 'Interpreter', 'none');
set(lg, 'FontName', FONT, 'Box', 'on');

if ~isempty(which('exportgraphics'))
    exportgraphics(fig, 'Fig11_SSLF_comparison.png', 'Resolution', 300);
    exportgraphics(fig, 'Fig11_SSLF_comparison.pdf', 'ContentType', 'vector');
    fprintf('Saved Fig11_SSLF_comparison.png / .pdf\n');
else
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, 'Fig11_SSLF_comparison.png', '-dpng', '-r300');
    fprintf('Saved Fig11_SSLF_comparison.png\n');
end


%% ============================================================ local functions
function [font, uni] = pick_font()
% One sans-serif font for all text (see plot_Figure2_Edge.m): an installed font with
% the superscript and subscript characters lets the units be written without TeX.
    font = 'Helvetica';  uni = false;
    if ispc || ismac
        if any(strcmpi(listfonts, 'Arial')), font = 'Arial'; end
        return
    end
    [st, out] = system('fc-list : family');
    if st ~= 0, return; end
    fams = strtrim(strsplit(out, {newline, ','}));
    [st2, out2] = system('fc-list '':charset=b0 b9 207b 2084 2212'' family');
    glyph = {};
    if st2 == 0, glyph = strtrim(strsplit(out2, {newline, ','})); end
    for c = {'FreeSans', 'DejaVu Sans', 'Noto Sans'}
        if any(strcmp(glyph, c{1})), font = c{1}; uni = true; return; end
    end
    for c = {'Liberation Sans', 'Nimbus Sans', 'Nimbus Sans L', 'Arial', 'Helvetica'}
        if any(strcmp(fams, c{1})), font = c{1}; return; end
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

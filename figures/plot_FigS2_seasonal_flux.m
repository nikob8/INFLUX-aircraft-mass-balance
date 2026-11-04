%% plot_FigS2_seasonal_flux.m
%  Figure S2: daily citywide CH4 emission rates (Edge approach) by meteorological season
%  (Reviewer 2, Line 493/502; referenced from Section 3.1.2).
%
%  Reads Indy_Mass_Balance_Table.xlsx (sheet FlightMasterSheet). The daily value is the
%  mean of the per-transect Edge fluxes on that day. Box = interquartile range, line =
%  median, whiskers = full range, points = individual flight days. A Kruskal-Wallis
%  test across seasons is computed inline (no toolbox needed).
%
%  Expected (spreadsheet of 2026-10-26): n = 5, 8, 13, 11 for DJF, MAM, JJA, SON;
%  means = 48, 55, 58, 52 mol/s; overall mean = 54 mol/s; Kruskal-Wallis p = 0.97.

clear; clc;
xlsFile = 'Indy_Mass_Balance_Table.xlsx';
outName = 'FigS2_seasonal_flux';

%% ------------------------------------------------------------ read
C = readcell(xlsFile, 'Sheet', 'FlightMasterSheet');
cLabel = 1; cEdge = 6;
dates = datetime.empty(0,1); sumEdge = []; nT = [];
k = 0;
for r = 2:size(C,1)
    lab = C{r, cLabel};
    if isa(lab,'datetime') || ((ischar(lab) || isstring(lab)) && ~isempty(regexp(strtrim(char(lab)), '^\d{4}/\d{2}/\d{2}$', 'once')))
        k = k + 1;
        if isa(lab,'datetime'), dates(k,1) = lab; else, dates(k,1) = datetime(strtrim(char(lab)), 'InputFormat','yyyy/MM/dd'); end
        sumEdge(k,1) = 0; nT(k,1) = 0;                                   %#ok<SAGROW>
    elseif (ischar(lab) || isstring(lab)) && ~isempty(regexpi(strtrim(char(lab)), '^Transect', 'once')) && k > 0
        v = C{r, cEdge};
        if isnumeric(v) && isscalar(v) && ~isnan(v)
            sumEdge(k) = sumEdge(k) + double(v);  nT(k) = nT(k) + 1;
        end
    end
end
dailyEdge = sumEdge ./ nT;      % daily mean of the per-transect Edge fluxes
assert(all(~isnan(dailyEdge)), 'some flight days have no Edge flux');

mo = month(dates);
sIdx = zeros(size(mo));
sIdx(ismember(mo,[12 1 2])) = 1; sIdx(ismember(mo,3:5)) = 2; sIdx(ismember(mo,6:8)) = 3; sIdx(ismember(mo,9:11)) = 4;
names = {'Winter (DJF)','Spring (MAM)','Summer (JJA)','Fall (SON)'};
groups = arrayfun(@(s) dailyEdge(sIdx == s), 1:4, 'UniformOutput', false);
pKW = kruskalWallisP(groups);
fprintf('n per season: %s | means: %s | overall mean %.0f | Kruskal-Wallis p = %.2f\n', ...
    mat2str(cellfun(@numel, groups)), mat2str(round(cellfun(@mean, groups))), mean(dailyEdge), pKW);

%% ------------------------------------------------------------ figure
fig = figure('Color','w','Units','pixels','Position',[100 100 640 460]); hold on; box on;
boxCol = [159 184 217]/255; ptCol = [31 78 156]/255; w = 0.55;
rng(1);
for s = 1:4
    v = sort(groups{s});
    q = prctileLocal(v, [25 50 75]);
    patch(s + [-1 1 1 -1]*w/2, [q(1) q(1) q(3) q(3)], boxCol, 'EdgeColor','k', 'LineWidth', 1);
    plot(s + [-1 1]*w/2, [q(2) q(2)], 'k-', 'LineWidth', 1.8);
    plot([s s], [q(3) max(v)], 'k-'); plot([s s], [min(v) q(1)], 'k-');
    plot(s + [-1 1]*w/4, [max(v) max(v)], 'k-'); plot(s + [-1 1]*w/4, [min(v) min(v)], 'k-');
    scatter(s + (rand(size(v))-0.5)*0.24, v, 22, ptCol, 'filled', 'MarkerEdgeColor','k', 'LineWidth', 0.3);
    text(s, -52, sprintf('n = %d, mean = %.0f', numel(v), mean(v)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize', 9.5);
end
hm = plot([0.4 4.6], mean(dailyEdge)*[1 1], 'r--', 'LineWidth', 1.3);
set(gca, 'XTick', 1:4, 'XTickLabel', names, 'FontSize', 11);
xlim([0.4 4.6]); ylim([-70 180]);
ylabel('Citywide CH_4 emission rate, Downwind Edge approach (mol s^{-1})');
legend(hm, sprintf('All flights: mean = %.0f mol s^{-1}', mean(dailyEdge)), 'Location','northwest', 'Box','on');
title(sprintf('Daily emission rates by season (%d flight days)', numel(dailyEdge)), 'FontWeight','normal');

set(fig, 'PaperPositionMode', 'auto');
print(fig, [outName '.png'], '-dpng', '-r300');
print(fig, [outName '.eps'], '-depsc');
fprintf('Saved %s.png / .eps\n', outName);

%% ============================================================ local functions
function q = prctileLocal(v, p)
    % linear-interpolation percentiles (same convention as MATLAB prctile)
    v = sort(v(:)); n = numel(v); q = zeros(size(p));
    for i = 1:numel(p)
        pos = p(i)/100 * n + 0.5;
        if pos <= 1, q(i) = v(1); elseif pos >= n, q(i) = v(n);
        else, lo = floor(pos); q(i) = v(lo) + (pos-lo) * (v(lo+1) - v(lo)); end
    end
end

function p = kruskalWallisP(groups)
    % Kruskal-Wallis H test with tie correction; chi-square approximation; base MATLAB
    all = []; g = [];
    for i = 1:numel(groups), all = [all; groups{i}(:)]; g = [g; i*ones(numel(groups{i}),1)]; end %#ok<AGROW>
    n = numel(all);
    [vs, idx] = sort(all); rk = zeros(n,1); i = 1;
    while i <= n
        j = i; while j < n && vs(j+1) == vs(i), j = j + 1; end
        rk(idx(i:j)) = (i+j)/2; i = j + 1;
    end
    H = 0;
    for i = 1:numel(groups)
        ni = sum(g == i); H = H + sum(rk(g == i))^2 / ni;
    end
    H = 12/(n*(n+1)) * H - 3*(n+1);
    [~, ~, ic] = unique(all); t = accumarray(ic, 1);
    H = H / (1 - sum(t.^3 - t)/(n^3 - n));
    df = numel(groups) - 1;
    p = 1 - gammainc(H/2, df/2);
end

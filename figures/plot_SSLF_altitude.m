%% plot_SSLF_altitude.m
%  SSLF (South Side Landfill) flux vs. transect altitude  -- new SI figure
%  (addresses Reviewer 1 comment on vertical mixing and Reviewer 2, Line 512)
%
%  Reads Indy_Mass_Balance_Table.xlsx (sheet FlightMasterSheet) from the
%  current folder and produces a three-panel figure:
%    (a) SSLF flux deviation vs. transect-height deviation, same convention
%        as Figure 8 (per-day means removed), for days with >= 2 quantified
%        landfill transects. Linear fit with slope +/- 1 SE and p-value.
%    (b) For days on which the plume was detected on some transects but not
%        on others: within-day height deviation of detected vs. undetected
%        transects (jittered dots, mean +/- SE), with a rank-sum p-value.
%    (c) Fraction of all 147 transects with a detectable plume, binned by
%        transect height / ABL height, with 95% Wilson intervals and a
%        logistic fit.
%
%  Requires only base MATLAB (R2019a or later for readcell). No toolboxes.
%
%  Expected output (from the spreadsheet of 2026-10-26), for checking:
%    (a) n = 80 transects on 25 days; slope = -0.72 +/- 0.57 mol/s per 100 m; p = 0.21
%    (b) 22 days; detected mean = -78 m (n = 58), undetected mean = +103 m (n = 44); p ~ 4.5e-5
%    (c) detection fraction by z/zi bin [0-.25 .25-.5 .5-.75 .75-1 >1]:
%        0.62 0.71 0.54 0.48 0.33  (n = 13 49 52 27 6); logistic slope = -1.63 +/- 0.72, p = 0.02

clear; clc;

xlsFile  = 'Indy_Mass_Balance_Table.xlsx';
outName  = 'FigS_SSLF_altitude';          % output file stem

%% ------------------------------------------------------------ read sheet
C = readcell(xlsFile, 'Sheet', 'FlightMasterSheet');

% column indices in the sheet (1-based)
cLabel = 1; cHeight = 3; cABL = 4; cLF = 10;

D = [];          % struct array: one element per flight day
k = 0;
for r = 2:size(C,1)
    lab = C{r, cLabel};
    if isDateCell(lab)
        k = k + 1;
        D(k).date = dateString(lab);
        D(k).z  = [];  D(k).zi = [];  D(k).lf = [];
    elseif isTransectCell(lab) && k > 0
        D(k).z (end+1) = toNum(C{r, cHeight});   % flight height AGL (m)
        D(k).zi(end+1) = toNum(C{r, cABL});      % ABL height (m)
        D(k).lf(end+1) = toNum(C{r, cLF});       % landfill flux (mol/s), NaN if 'N/A'
    end
end
nDays = numel(D);  nTr = sum(arrayfun(@(d) numel(d.z), D));
fprintf('Parsed %d flight days, %d transects, %d with a quantified SSLF plume on %d days\n', ...
    nDays, nTr, sum(arrayfun(@(d) sum(~isnan(d.lf)), D)), sum(arrayfun(@(d) any(~isnan(d.lf)), D)));

%% ------------------------------------------------------------ panel (a)
xa = []; ya = []; nDaysA = 0;
for d = 1:nDays
    q = ~isnan(D(d).lf);
    if sum(q) < 2, continue; end
    nDaysA = nDaysA + 1;
    xa = [xa, D(d).z(q)  - mean(D(d).z(q))];   %#ok<AGROW>
    ya = [ya, D(d).lf(q) - mean(D(d).lf(q))];  %#ok<AGROW>
end
[slope, seSlope, pA, intercept] = olsFit(xa, ya);
fprintf('(a) n = %d transects on %d days; slope = %.2f +/- %.2f mol/s per 100 m; p = %.3f\n', ...
    numel(xa), nDaysA, 100*slope, 100*seSlope, pA);

%% ------------------------------------------------------------ panel (b)
hDet = []; hNo = []; nDaysB = 0;
for d = 1:nDays
    det = ~isnan(D(d).lf);
    if ~(any(det) && any(~det)), continue; end
    nDaysB = nDaysB + 1;
    dev = D(d).z - mean(D(d).z);
    hDet = [hDet, dev(det)];   %#ok<AGROW>
    hNo  = [hNo,  dev(~det)];  %#ok<AGROW>
end
pB = ranksumP(hDet, hNo);
fprintf('(b) %d days; detected mean = %+.0f m (n = %d); undetected mean = %+.0f m (n = %d); rank-sum p = %.1e\n', ...
    nDaysB, mean(hDet), numel(hDet), mean(hNo), numel(hNo), pB);

%% ------------------------------------------------------------ panel (c)
zAll  = [D.z];  ziAll = [D.zi];  detAll = ~isnan([D.lf]);
zn = zAll ./ ziAll;                                   % height normalized by ABL depth
edges = [0 0.25 0.5 0.75 1.0 1.5];
nb = numel(edges) - 1;
frac = nan(1,nb); lo = frac; hi = frac; nBin = frac; cen = frac;
for b = 1:nb
    m = zn >= edges(b) & zn < edges(b+1);
    nBin(b) = sum(m);  kk = sum(detAll(m));
    [frac(b), lo(b), hi(b)] = wilson(kk, nBin(b));
    cen(b) = mean(edges(b:b+1));
end
[bLog, seLog, pC] = logisticFit(zn(:), double(detAll(:)));
fprintf('(c) detection fraction by bin: %s  (n = %s); logistic slope = %.2f +/- %.2f, p = %.3f\n', ...
    mat2str(round(frac,2)), mat2str(nBin), bLog(2), seLog(2), pC);

%% ------------------------------------------------------------ figure
fig = figure('Color','w','Units','pixels','Position',[100 100 1500 440]);
blue = [31 78 156]/255;  red = [192 57 43]/255;  lightblue = [159 184 217]/255;

% (a)
subplot(1,3,1); hold on; box on;
scatter(xa, ya, 24, blue, 'filled');
xx = [-650 650];
plot(xx, intercept + slope*xx, '-', 'Color', 'r', 'LineWidth', 1.5);
plot(xx, [0 0], 'k-', 'LineWidth', 0.5);
xlim(xx);
xlabel('Meters from mean transect height');
ylabel('SSLF flux deviation (mol s^{-1})');
title(sprintf('(a) SSLF flux vs. altitude (n = %d transects, %d days)', numel(xa), nDaysA), 'FontWeight','normal');
legend(sprintf('Slope = (%.2f \\pm %.2f) per 100 m, p = %.2f', 100*slope, 100*seSlope, pA), ...
       'Location','northeast', 'Box','on');
set(gca,'FontSize',11);

% (b)
subplot(1,3,2); hold on; box on;
rng(1);                                               % reproducible jitter
jit = @(n) (rand(1,n) - 0.5) * 0.36;
scatter(0 + jit(numel(hDet)), hDet, 20, blue, 'filled', 'MarkerFaceAlpha', 0.7);
scatter(1 + jit(numel(hNo)),  hNo,  20, red,  'filled', 'MarkerFaceAlpha', 0.7);
plot([-0.16 0.16], mean(hDet)*[1 1], 'k-', 'LineWidth', 2.5);
plot([ 0.84 1.16], mean(hNo) *[1 1], 'k-', 'LineWidth', 2.5);
errorbar(0, mean(hDet), std(hDet)/sqrt(numel(hDet)), 'k', 'LineStyle','none', 'LineWidth', 1.5, 'CapSize', 8);
errorbar(1, mean(hNo),  std(hNo) /sqrt(numel(hNo)),  'k', 'LineStyle','none', 'LineWidth', 1.5, 'CapSize', 8);
plot([-0.6 1.6], [0 0], 'k-', 'LineWidth', 0.5);
xlim([-0.6 1.6]);
set(gca, 'XTick', [0 1], 'XTickLabel', {sprintf('Plume detected\\newline(n = %d)', numel(hDet)), ...
                                        sprintf('No plume\\newline(n = %d)', numel(hNo))});
ylabel('Meters from mean transect height');
title(sprintf('(b) Same-day transect height (%d days), p = %.0e', nDaysB, pB), 'FontWeight','normal');
set(gca,'FontSize',11);

% (c)
subplot(1,3,3); hold on; box on;
for b = 1:nb
    w = 0.9*(edges(b+1)-edges(b));
    patch(cen(b) + [-1 1 1 -1]*w/2, [0 0 frac(b) frac(b)], lightblue, 'EdgeColor','k');
    text(cen(b), 0.03, sprintf('n=%d', nBin(b)), 'HorizontalAlignment','center', 'FontSize', 9);
end
errorbar(cen, frac, frac-lo, hi-frac, 'k.', 'LineStyle','none', 'LineWidth', 1.2, 'CapSize', 6);
zz = linspace(0, 1.5, 200);
hL = plot(zz, 1./(1+exp(-(bLog(1) + bLog(2)*zz))), 'r-', 'LineWidth', 1.5);
xlim([0 1.5]); ylim([0 1]);
xlabel('Transect height / ABL height');
ylabel('Fraction of transects with detectable plume');
title(sprintf('(c) Plume detection vs. height in ABL (all %d transects)', nTr), 'FontWeight','normal');
legend(hL, sprintf('Logistic fit, p = %.2f', pC), 'Location','northeast', 'Box','on');
set(gca,'FontSize',11);

%% ------------------------------------------------------------ save
set(fig, 'PaperPositionMode', 'auto');
print(fig, [outName '.png'], '-dpng', '-r300');
print(fig, [outName '.eps'], '-depsc');
if ~isempty(which('exportgraphics'))
    exportgraphics(fig, [outName '.pdf'], 'ContentType','vector');
end
fprintf('Saved %s.png / .eps\n', outName);

%% ============================================================ local functions
function tf = isDateCell(x)
    tf = false;
    if isa(x,'datetime'), tf = true; return; end
    if ischar(x) || isstring(x)
        tf = ~isempty(regexp(strtrim(char(x)), '^\d{4}/\d{2}/\d{2}$', 'once'));
    end
end

function s = dateString(x)
    if isa(x,'datetime'), s = datestr(x, 'yyyy/mm/dd'); else, s = strtrim(char(x)); end
end

function tf = isTransectCell(x)
    tf = (ischar(x) || isstring(x)) && ~isempty(regexpi(strtrim(char(x)), '^Transect', 'once'));
end

function v = toNum(x)
    % numeric cell -> value; text ('N/A', ' ', '') or missing -> NaN
    if isnumeric(x) && isscalar(x), v = double(x);
    elseif ischar(x) || isstring(x), v = str2double(strtrim(char(x)));
    else, v = NaN;
    end
end

function [b1, se1, p, b0] = olsFit(x, y)
    x = x(:); y = y(:); n = numel(x);
    X = [ones(n,1) x];
    beta = X \ y;
    res = y - X*beta;
    s2 = sum(res.^2) / (n-2);
    cov = s2 * inv(X'*X);               %#ok<MINV>
    b0 = beta(1); b1 = beta(2); se1 = sqrt(cov(2,2));
    t = b1/se1; df = n-2;
    p = betainc(df/(df + t^2), df/2, 0.5);   % two-sided t-test, base MATLAB
end

function p = ranksumP(a, b)
    % Wilcoxon rank-sum (Mann-Whitney) two-sided p, normal approximation with tie correction
    a = a(:); b = b(:); na = numel(a); nb = numel(b); n = na + nb;
    v = [a; b];
    [vs, idx] = sort(v); rk = zeros(n,1);
    i = 1;
    while i <= n
        j = i; while j < n && vs(j+1) == vs(i), j = j + 1; end
        rk(idx(i:j)) = (i + j)/2;  i = j + 1;
    end
    Ra = sum(rk(1:na));
    U  = Ra - na*(na+1)/2;
    mu = na*nb/2;
    [~, ~, ic] = unique(v); tcount = accumarray(ic, 1);
    tie = sum(tcount.^3 - tcount) / (n*(n-1));
    sig = sqrt(na*nb/12 * ((n+1) - tie));
    z = (U - mu) / sig;
    p = erfc(abs(z)/sqrt(2));
end

function [f, lo, hi] = wilson(k, n)
    z = 1.96; p = k/n;
    c = (p + z^2/(2*n)) / (1 + z^2/n);
    h = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n);
    f = p; lo = c - h; hi = c + h;
end

function [b, se, p] = logisticFit(x, y)
    % logistic regression y ~ b0 + b1*x by Newton-Raphson; p for b1 (Wald)
    X = [ones(size(x)) x]; b = zeros(2,1);
    for it = 1:50
        pr = 1./(1+exp(-X*b)); W = pr.*(1-pr);
        H = X' * (W .* X);  g = X' * (y - pr);
        step = H \ g; b = b + step;
        if max(abs(step)) < 1e-10, break; end
    end
    pr = 1./(1+exp(-X*b)); W = pr.*(1-pr);
    cov = inv(X' * (W .* X));            %#ok<MINV>
    se = sqrt(diag(cov));
    z = b(2)/se(2);
    p = erfc(abs(z)/sqrt(2));
end

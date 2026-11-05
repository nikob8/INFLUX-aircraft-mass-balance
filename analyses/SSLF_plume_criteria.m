%% SSLF_plume_criteria.m
%  Objective criteria for a "well-defined" South Side Landfill (SSLF) plume
%  (Reviewer 2, Line 326) and their comparison with the manual selection of
%  Section 2.4.  Applied to all downwind transects of the 37 flight days.
%
%  Definition.  On each downwind transect the CH4 series (15 s running mean,
%  as in the landfill mass balance) within WIN_KM of the projected SSLF
%  position is fitted by least squares with a Gaussian spike on a smooth
%  polynomial background (order NPOLY, which represents the broad city plume
%  and the regional gradient):
%        CH4(s) = p(s) + A exp(-(s - mu)^2 / (2 sigma^2))
%  The spike centre mu is searched within W_KM of the projected position and
%  sigma over SIG_KM; the fit with the smallest residual is kept.  The
%  transect has a well-defined landfill plume when
%        A  >= PMIN_PPB          minimum peak height above the local background
%        A  >= SNRMIN * s_res    signal to noise: s_res is the root mean square
%                                residual of the fit, i.e. the variability of
%                                the transect that neither the background nor
%                                the spike explains
%        FWHM = 2.355 sigma <= FWMAX_KM   the feature is narrow (a point
%                                source plume, not the city plume itself)
%        P  >= PMIN_PPB          the measured CH4 between the edges also rises
%                                at least PMIN above the edge background
%  The plume edges used for the emission rate are mu +/- KEDGE*sigma (KEDGE =
%  2.5 encloses 99% of the area of the fitted spike); as in
%  the manual method the background is the mean of the CH4 values at the two
%  edges, E is the mean enhancement between them and D their distance (edges
%  beyond the end of the transect are moved to the end, as an analyst would).
%  Only positive spikes are considered (fits with A < 0 are discarded).  The
%  projected SSLF position is the crossing of the mean-wind line through the
%  landfill with the transect (same construction as the county corners in the
%  Downwind Edge method).  The thresholds (10 ppb, 7, 10 km) sit at the 10th
%  (amplitude, signal to noise) and 90th (width) percentiles of the plumes the
%  analysts quantified, i.e. they describe the manual selection; the
%  sensitivity section prints the result for other thresholds, including the
%  permissive set RELAXED (Text S3 and Table S2 of the Supporting Information).
%
%  Output: SSLF_plume_criteria.csv (one row per transect: fit parameters,
%  classification by the criteria and by the analysts, both emission rates,
%  position of the manual plume), a printout with the agreement between the
%  two classifications, the emission rates with both edge definitions, the
%  daily SSLF means and sensitivity tables over the thresholds,
%  SSLF_criteria.png, and two inspection figures with the CH4 series and the
%  fit of every transect: SSLF_criteria_disagree.png (transects on which the
%  criteria and the analysts differ) and SSLF_criteria_agree.png.
%
%  Inputs: FlightList.mat, LegAndProfileList/<Date>.mat, INFLUX_<Date>,
%  InterX.m, <Date>/Solution-MethodA.mat (mean wind, density, ABL, cos theta
%  per transect), Landfill/<Date>_Loop<k>.mat (manual plumes: CH4LF series,
%  used to locate the manual edges on the transect) and
%  Landfill_transects.csv (manual emission rates).

clear; clc; close all

%% ------------------------------------------------------------ settings
MB_ROOT  = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';   % FlightList.mat, LegAndProfileList/
OBS_DIR  = '/storage/work/zrb5027/WRFGUI_V3/Observations/INFLUX';               % INFLUX_<yyyymmdd>
CODE_DIR = '/storage/work/zrb5027/WRFGUI_V3/Code';                              % InterX.m
SOL_ROOT = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB';   % <Date>/Solution-MethodA.mat
LF_DIR   = fullfile(pwd, 'Landfill');                 % Loop files written by MassBalanceLandfill.m
LF_CSV   = fullfile(pwd, 'Landfill_transects.csv');   % from Collect_Landfill.m (manual emission rates)
OUT_DIR  = pwd;
SSLF_LATLON = [39.714691162109375, -86.205429077148438];
prm.W     = 8;              % km: the spike centre is searched within W of the projected SSLF position
prm.WIN   = 16;             % km: half-width of the segment used for the fit
prm.NPOLY = 4;              % order of the polynomial background
prm.SIG   = 0.5:0.25:7;     % km: candidate sigma of the spike
prm.PMIN  = 10;             % ppb: minimum amplitude
prm.SNRMIN = 7;             % minimum amplitude / rms residual of the fit
prm.FWMAX = 10;             % km: maximum full width at half maximum
prm.KEDGE = 2.5;            % edges at mu +/- KEDGE*sigma
SENS_SNR = [4 5 7 8];  SENS_FWMAX = [8 10 15];  SENS_PMIN = [5 10 15];  SENS_KEDGE = [2 2.5 3];
RELAXED  = [5 4 15];        % permissive thresholds [PMIN SNRMIN FWMAX] reported in Text S3

addpath(MB_ROOT); addpath(OBS_DIR); addpath(CODE_DIR);
assert(~isempty(which('InterX')), 'InterX.m not on the path (check CODE_DIR)');
FL = load(fullfile(MB_ROOT, 'FlightList.mat')); FlightList = FL.FlightList;
man = readtable(LF_CSV, 'TextType', 'string', 'DatetimeType', 'text');

%% ------------------------------------------------------------ loop over all downwind transects
R = []; CUR = {};  t0 = tic;
for k = 1:numel(FlightList)
    Date = FlightList(k).name; Date = Date(1:8);
    L = load(fullfile(MB_ROOT, 'LegAndProfileList', Date));
    x = load(['INFLUX_' Date]); x(x < -998) = NaN;
    Lat = x(:,2); Lon = x(:,3);
    CH4 = movmean(x(:,6), 15, 'omitnan');            % ppm, as in MassBalanceLandfill.m
    if isfield(L, 'CH4Flagged'), CH4 = L.CH4Flagged(:); end
    SA = load(fullfile(SOL_ROOT, Date, 'Solution-MethodA.mat'));
    MeanU = -SA.MeanWSPD*sind(SA.MeanWDIR); MeanV = -SA.MeanWSPD*cosd(SA.MeanWDIR);
    PC = [SSLF_LATLON(1)+MeanV, SSLF_LATLON(2)+MeanU*(111/86)];
    PA = [SSLF_LATLON(1)-MeanV, SSLF_LATLON(2)-MeanU*(111/86)];
    LineLF = [PA(2) SSLF_LATLON(2) PC(2); PA(1) SSLF_LATLON(1) PC(1)];     % [lon; lat] of the mean-wind line through the SSLF
    dstr = string(sprintf('%s-%s-%s', Date(1:4), Date(5:6), Date(7:8)));
    for j = 1:size(L.LegList, 1)
        Leg = L.LegList(j,1):L.LegList(j,2);
        lat = Lat(Leg); lon = Lon(Leg); c = CH4(Leg)*1000;               % ppb
        ok = ~isnan(lat) & ~isnan(lon) & ~isnan(c);
        s = NaN(size(lat)); la0 = mean(lat, 'omitnan');
        dx = [0; diff(lon(ok))*111.32*cosd(la0)]; dy = [0; diff(lat(ok))*110.57];
        s(ok) = cumsum(hypot(dx, dy));                                    % along-track distance (km)
        r = struct('Date', dstr, 'Transect', j, 's0_km', NaN, 'mu_km', NaN, 'sigma_km', NaN, 'fwhm_km', NaN, ...
            'A_ppb', NaN, 'sres_ppb', NaN, 'SNR', NaN, 'P_ppb', NaN, 'E_ppb', NaN, 'D_km', NaN, 'sL_km', NaN, 'sR_km', NaN, ...
            'edges_in', false, 'det_alg', false, 'reason', "", 'det_manual', false, 'Flux_alg_mols', NaN, 'Flux_alg_k20_mols', NaN, 'Flux_alg_k30_mols', NaN, ...
            'Flux_manual_mols', NaN, 'sManStart_km', NaN, 'sManEnd_km', NaN, 'sManPeak_km', NaN, 'cos_theta', SA.CosTheta(j));
        % ---- manual plume (if any): locate the CH4LF series of the Loop file on the transect
        mi = find(man.Date == r.Date & man.Transect == j, 1);
        if ~isempty(mi)
            r.det_manual = true; r.Flux_manual_mols = man.Flux_Landfill_mols(mi);
            lf = fullfile(LF_DIR, sprintf('%s_Loop%d.mat', Date, j));
            if isfile(lf)
                Q = load(lf);
                if isfield(Q, 'CH4LF')
                    v = Q.CH4LF(:)*1000; n = numel(v); i0 = NaN;
                    for i = 1:numel(c)-n+1
                        if abs(c(i) - v(1)) < 1e-6 && max(abs(c(i:i+n-1) - v)) < 1e-6, i0 = i; break; end
                    end
                    if ~isnan(i0)
                        r.sManStart_km = s(i0); r.sManEnd_km = s(i0+n-1);
                        [~, ip] = max(v); r.sManPeak_km = s(i0+ip-1);
                    end
                end
            end
        end
        % ---- projected SSLF position on the transect
        X = InterX(LineLF, [lon(ok)'; lat(ok)']);
        if isempty(X), r.reason = "projection outside the transect"; R = [R; r]; CUR{end+1} = struct('s', s, 'c', c, 'sf', [], 'fit', [], 'bg', []); continue; end %#ok<AGROW>
        [~, i0] = min(hypot((lon - X(1,1))*111.32*cosd(la0), (lat - X(2,1))*110.57));
        r.s0_km = s(i0);
        % ---- fit and criteria
        f = detect_plume(s, c, r.s0_km, prm);
        r.mu_km = f.mu_km; r.sigma_km = f.sigma_km; r.fwhm_km = f.fwhm_km; r.A_ppb = f.A_ppb; r.sres_ppb = f.sres_ppb; r.SNR = f.SNR;
        r.det_alg = f.det; r.reason = f.reason; r.edges_in = f.edges_in;
        CUR{end+1} = struct('s', s, 'c', c, 'sf', f.ss, 'fit', f.fit, 'bg', f.bg); %#ok<SAGROW>
        if ~isnan(f.mu_km)
            q = edge_quantities(s, c, f.mu_km, f.sigma_km, prm.KEDGE);
            r.P_ppb = q.P; r.E_ppb = q.E; r.D_km = q.D; r.sL_km = q.sL; r.sR_km = q.sR;
            if r.det_alg && q.P < prm.PMIN, r.det_alg = false; r.reason = "no bump in the data"; end
            flx = @(qq) SA.Density(j)*SA.MeanWSPD*(qq.E/1e9)*(qq.D*1000)*SA.DownwindABL(j)*SA.CosTheta(j);   % Equation 1 with cos theta, as Table S1
            r.Flux_alg_mols = flx(q);
            r.Flux_alg_k20_mols = flx(edge_quantities(s, c, f.mu_km, f.sigma_km, 2));
            r.Flux_alg_k30_mols = flx(edge_quantities(s, c, f.mu_km, f.sigma_km, 3));
        end
        R = [R; r]; %#ok<AGROW>
    end
    fprintf('%2d %s: %d transects, manual plumes %d, objective %d (%.0f s)\n', k, Date, size(L.LegList,1), nnz([R.det_manual] & [R.Date] == dstr), nnz([R.det_alg] & [R.Date] == dstr), toc(t0));
end
T = struct2table(R);
writetable(T, fullfile(OUT_DIR, 'SSLF_plume_criteria.csv'));

%% ------------------------------------------------------------ agreement with the manual selection
n = height(T); a = T.det_alg; m = T.det_manual;
fprintf('\n%d transects: manual plumes %d, objective plumes %d\n', n, nnz(m), nnz(a));
fprintf('agreement %d of %d (%.0f%%): both %d, neither %d, manual only %d, objective only %d\n', nnz(a == m), n, 100*mean(a == m), nnz(a & m), nnz(~a & ~m), nnz(m & ~a), nnz(a & ~m));
mo = m & ~a;
if any(mo)
    fprintf('manual only (reason the criteria fail):\n');
    for i = find(mo)', fprintf('   %s T%d: %s (A %.1f ppb, SNR %.1f, FWHM %.1f km, manual peak - projection %.1f km, manual %.1f mol/s)\n', T.Date(i), T.Transect(i), T.reason(i), T.A_ppb(i), T.SNR(i), T.fwhm_km(i), T.sManPeak_km(i) - T.s0_km(i), T.Flux_manual_mols(i)); end
end
ao = a & ~m;
if any(ao)
    fprintf('objective only:\n');
    for i = find(ao)', fprintf('   %s T%d: A %.1f ppb, SNR %.1f, FWHM %.1f km, mu - projection %.1f km, %.1f mol/s\n', T.Date(i), T.Transect(i), T.A_ppb(i), T.SNR(i), T.fwhm_km(i), T.mu_km(i) - T.s0_km(i), T.Flux_alg_mols(i)); end
end
both = a & m;
if any(both)
    d = T.Flux_alg_mols(both) - T.Flux_manual_mols(both); rat = T.Flux_alg_mols(both)./T.Flux_manual_mols(both);
    fprintf('both detected (%d): objective minus manual emission rate: median %.1f mol/s (IQR %.1f to %.1f), median ratio %.2f, r = %.2f\n', nnz(both), median(d), prctile(d,25), prctile(d,75), median(rat), corr(T.Flux_alg_mols(both), T.Flux_manual_mols(both)));
    dp = T.mu_km(both) - T.sManPeak_km(both); dw = T.D_km(both) - (T.sManEnd_km(both) - T.sManStart_km(both));
    fprintf('spike centre minus manual peak: median %.1f km, 90th percentile of |d| %.1f km; width objective minus manual: median %.1f km (IQR %.1f to %.1f)\n', median(dp,'omitnan'), prctile(abs(dp),90), median(dw,'omitnan'), prctile(dw,25), prctile(dw,75));
    for kk = [2 3]
        if kk == 2, fa = T.Flux_alg_k20_mols(both); else, fa = T.Flux_alg_k30_mols(both); end
        fprintf('   edges at %.1f sigma: median ratio to manual %.2f\n', kk, median(fa./T.Flux_manual_mols(both)));
    end
end
days = unique(T.Date);
dm = NaN(numel(days),1); da = dm;
for i = 1:numel(days)
    v = T.Flux_manual_mols(T.Date == days(i) & m); if ~isempty(v), dm(i) = mean(v); end
    v = T.Flux_alg_mols(T.Date == days(i) & a);    if ~isempty(v), da(i) = mean(v); end
end
fprintf('daily SSLF means: manual %d days, mean %.1f, median %.1f, 2 SE %.1f; objective %d days, mean %.1f, median %.1f, 2 SE %.1f\n', ...
    nnz(~isnan(dm)), mean(dm,'omitnan'), median(dm,'omitnan'), 2*std(dm,'omitnan')/sqrt(nnz(~isnan(dm))), ...
    nnz(~isnan(da)), mean(da,'omitnan'), median(da,'omitnan'), 2*std(da,'omitnan')/sqrt(nnz(~isnan(da))));
fprintf('manual plumes: A median %.1f ppb (10th pct %.1f), SNR median %.1f (10th pct %.1f), FWHM median %.1f km (90th pct %.1f); other transects: A median %.1f ppb, SNR median %.1f, FWHM median %.1f km\n', ...
    median(T.A_ppb(m),'omitnan'), prctile(T.A_ppb(m),10), median(T.SNR(m),'omitnan'), prctile(T.SNR(m),10), median(T.fwhm_km(m),'omitnan'), prctile(T.fwhm_km(m),90), median(T.A_ppb(~m),'omitnan'), median(T.SNR(~m),'omitnan'), median(T.fwhm_km(~m),'omitnan'));

%% ------------------------------------------------------------ sensitivity to the thresholds
dayMean = @(sel, col) mean(arrayfun(@(i) mean(col(T.Date == days(i) & sel)), 1:numel(days)), 'omitnan');
crit = @(p, sn, fw) T.A_ppb >= p & T.SNR >= sn & T.fwhm_km <= fw & T.P_ppb >= p;
fprintf('\nsensitivity (SNRMIN x FWMAX_KM, PMIN = %g ppb): n detected / agreement with manual / daily mean SSLF (mol/s)\n', prm.PMIN);
fprintf('%14s', ''); for fw = SENS_FWMAX, fprintf('%24s', sprintf('FWMAX = %g km', fw)); end; fprintf('\n');
for sn = SENS_SNR
    fprintf('SNRMIN = %3g   ', sn);
    for fw = SENS_FWMAX
        aa = crit(prm.PMIN, sn, fw);
        fprintf('%9d /%4.0f%% /%7.1f', nnz(aa), 100*mean(aa == m), dayMean(aa, T.Flux_alg_mols));
    end
    fprintf('\n');
end
fprintf('PMIN (SNRMIN = %g, FWMAX = %g):', prm.SNRMIN, prm.FWMAX);
for p = SENS_PMIN, aa = crit(p, prm.SNRMIN, prm.FWMAX); fprintf('   %g ppb: %d / %.0f%% / %.1f', p, nnz(aa), 100*mean(aa == m), dayMean(aa, T.Flux_alg_mols)); end; fprintf('\n');
aa = crit(RELAXED(1), RELAXED(2), RELAXED(3));
fprintf('relaxed thresholds (A >= %g ppb, SNR >= %g, FWHM <= %g km): %d detected, agreement %.0f%% (both %d, neither %d, manual only %d, objective only %d), daily mean SSLF %.1f mol/s (%d days); the additional detections have A median %.1f ppb, SNR median %.1f, emission rate median %.1f mol/s\n', ...
    RELAXED, nnz(aa), 100*mean(aa == m), nnz(aa & m), nnz(~aa & ~m), nnz(m & ~aa), nnz(aa & ~m), dayMean(aa, T.Flux_alg_mols), nnz(arrayfun(@(i) any(aa & T.Date == days(i)), 1:numel(days))), ...
    median(T.A_ppb(aa & ~m)), median(T.SNR(aa & ~m)), median(T.Flux_alg_mols(aa & ~m)));
fprintf('edges at KEDGE sigma (default criteria): daily mean SSLF %.1f (2), %.1f (2.5), %.1f (3) mol/s; on the transects detected by both: median ratio to manual %.2f (2), %.2f (2.5), %.2f (3)\n', ...
    dayMean(a, T.Flux_alg_k20_mols), mean(da, 'omitnan'), dayMean(a, T.Flux_alg_k30_mols), ...
    median(T.Flux_alg_k20_mols(both)./T.Flux_manual_mols(both)), median(T.Flux_alg_mols(both)./T.Flux_manual_mols(both)), median(T.Flux_alg_k30_mols(both)./T.Flux_manual_mols(both)));

%% ------------------------------------------------------------ figure
FN = 'Arial';
fig = figure('Color', 'w', 'Position', [100 100 1250 520]);
ax1 = subplot(1,2,1); hold(ax1, 'on'); box(ax1, 'on');
okf = ~isnan(T.A_ppb); wide = okf & T.fwhm_km > prm.FWMAX;
h1 = scatter(ax1, T.SNR(okf & ~m), max(T.A_ppb(okf & ~m), 2), 45, [0.8 0.3 0.3], 'o', 'LineWidth', 1.2);
h2 = scatter(ax1, T.SNR(okf & m), max(T.A_ppb(okf & m), 2), 45, [0.2 0.4 0.75], 'filled');
h4 = scatter(ax1, T.SNR(wide), max(T.A_ppb(wide), 2), 30, 'k', 'x', 'LineWidth', 1.2);
h3 = plot(ax1, [prm.SNRMIN prm.SNRMIN], [2 1200], 'k--', 'LineWidth', 1.2); plot(ax1, [1 50], [prm.PMIN prm.PMIN], 'k--', 'LineWidth', 1.2);
set(ax1, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12, 'FontName', FN, 'TickDir', 'out'); xlim(ax1, [1 50]); ylim(ax1, [2 1200]);
set(ax1, 'XTick', [1 2 5 10 20 50], 'YTick', [2 5 10 20 50 100 200 500 1000]);
xlabel(ax1, 'Amplitude of the fitted spike / rms residual of the fit', 'FontSize', 13, 'FontName', FN); ylabel(ax1, 'Amplitude of the fitted spike (ppb)', 'FontSize', 13, 'FontName', FN);
legend(ax1, [h2 h1 h4 h3], {sprintf('plume quantified by the analysts (%d)', nnz(okf & m)), sprintf('no plume quantified (%d)', nnz(okf & ~m)), ...
    sprintf('width above %g km (%d)', prm.FWMAX, nnz(wide)), sprintf('thresholds: %g ppb, %g times the residual', prm.PMIN, prm.SNRMIN)}, ...
    'Location', 'northwest', 'FontSize', 10, 'FontName', FN);    % upper left of the log-log panel is empty
title(ax1, '(a) Best-fit spike on every transect', 'FontSize', 13, 'FontName', FN);
ax2 = subplot(1,2,2); hold(ax2, 'on'); box(ax2, 'on');
lim = [0 max([T.Flux_alg_mols(both); T.Flux_manual_mols(both)])*1.05];
plot(ax2, lim, lim, 'k--', 'LineWidth', 1);
scatter(ax2, T.Flux_manual_mols(both), T.Flux_alg_mols(both), 45, [0.2 0.4 0.75], 'filled');
set(ax2, 'FontSize', 12, 'FontName', FN, 'TickDir', 'out'); xlim(ax2, lim); ylim(ax2, lim); axis(ax2, 'square');
xlabel(ax2, 'SSLF emission rate, edges chosen by the analysts (mol s^{-1})', 'FontSize', 13, 'FontName', FN); ylabel(ax2, 'SSLF emission rate, fitted edges (mol s^{-1})', 'FontSize', 13, 'FontName', FN);
title(ax2, sprintf('(b) Transects selected by both (n = %d)', nnz(both)), 'FontSize', 13, 'FontName', FN);
print(fig, fullfile(OUT_DIR, 'SSLF_criteria.png'), '-dpng', '-r300');

%% ------------------------------------------------------------ inspection figures: every transect with its fit
groups = {find(a ~= m), 'SSLF_criteria_disagree.png', 'criteria and analysts differ'; find(a & m), 'SSLF_criteria_agree.png', 'plume by both'};
for g = 1:size(groups, 1)
    idx = groups{g, 1}; if isempty(idx), continue; end
    nc = 6; nr = ceil(numel(idx)/nc);
    fg = figure('Color', 'w', 'Position', [50 50 2200 320*nr], 'Visible', 'off');
    for q = 1:numel(idx)
        i = idx(q); ax = subplot(nr, nc, q); hold(ax, 'on'); box(ax, 'on');
        u = CUR{i};
        plot(ax, u.s, u.c, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.8);
        if ~isempty(u.fit)
            plot(ax, u.sf, u.bg, '-', 'Color', [0.1 0.6 0.2], 'LineWidth', 1.4);
            plot(ax, u.sf, u.fit, '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.4);
        end
        yl = [min(u.c, [], 'omitnan') max(u.c, [], 'omitnan')]; yl = yl + [-0.05 0.08]*diff(yl);
        plot(ax, [T.s0_km(i) T.s0_km(i)], yl, 'k-', 'LineWidth', 1);
        if ~isnan(T.sL_km(i)), plot(ax, [T.sL_km(i) T.sL_km(i)], yl, 'r--', 'LineWidth', 1); plot(ax, [T.sR_km(i) T.sR_km(i)], yl, 'r--', 'LineWidth', 1); end
        if ~isnan(T.sManStart_km(i)), plot(ax, [T.sManStart_km(i) T.sManStart_km(i)], yl, 'b:', 'LineWidth', 1.4); plot(ax, [T.sManEnd_km(i) T.sManEnd_km(i)], yl, 'b:', 'LineWidth', 1.4); end
        ylim(ax, yl); xlim(ax, [max(0, T.s0_km(i) - 25) T.s0_km(i) + 25]);
        tag = ''; if T.det_manual(i), tag = [tag 'manual ' sprintf('%.0f', T.Flux_manual_mols(i)) ' | ']; end
        if T.det_alg(i), tag = [tag 'objective ' sprintf('%.0f', T.Flux_alg_mols(i)) ' mol/s']; else, tag = [tag char(T.reason(i))]; end
        title(ax, sprintf('%s T%d  A %.0f ppb, SNR %.1f, FWHM %.1f km\n%s', T.Date(i), T.Transect(i), T.A_ppb(i), T.SNR(i), T.fwhm_km(i), tag), 'FontSize', 9, 'FontWeight', 'normal');
        set(ax, 'FontSize', 8);
    end
    print(fg, fullfile(OUT_DIR, groups{g, 2}), '-dpng', '-r110'); close(fg);
    fprintf('wrote %s (%d transects, %s): grey CH4, green background, red fit, black projected SSLF, red dashed objective edges, blue dotted manual edges\n', groups{g, 2}, numel(idx), groups{g, 3});
end
fprintf('wrote SSLF_plume_criteria.csv and SSLF_criteria.png in %s\n', OUT_DIR);

%% ------------------------------------------------------------ local functions
function r = detect_plume(s, c, s0, prm)
% Gaussian spike on a polynomial background, fitted within prm.WIN of s0.
    r = struct('mu_km', NaN, 'sigma_km', NaN, 'A_ppb', NaN, 'sres_ppb', NaN, 'SNR', NaN, 'fwhm_km', NaN, 'edges_in', false, 'det', false, 'reason', "", 'ss', [], 'fit', [], 'bg', []);
    ok = ~isnan(s) & ~isnan(c) & abs(s - s0) <= prm.WIN;
    if nnz(ok) < 20, r.reason = "projection outside the transect"; return; end
    ss = s(ok); cc = c(ok); n = numel(ss);
    X0 = ((ss - s0)/prm.WIN).^(0:prm.NPOLY);                 % polynomial background basis
    best = struct('rss', Inf);
    for mu = (s0 - prm.W):0.25:(s0 + prm.W)
        for sg = prm.SIG
            X = [X0, exp(-0.5*((ss - mu)/sg).^2)];
            beta = X \ cc; if beta(end) <= 0, continue; end   % positive spikes only
            res = cc - X*beta; rss = res'*res;
            if rss < best.rss, best = struct('rss', rss, 'mu', mu, 'sg', sg, 'beta', beta, 'X', X); end
        end
    end
    if isinf(best.rss), r.reason = "no positive spike"; return; end
    A = best.beta(end);
    r.ss = ss; r.fit = best.X*best.beta; r.bg = best.X(:,1:end-1)*best.beta(1:end-1);
    dof = max(n - size(best.X, 2), 1); sres = sqrt(best.rss/dof);        % rms residual of the fit
    r.mu_km = best.mu; r.sigma_km = best.sg; r.A_ppb = A; r.sres_ppb = sres; r.SNR = A/sres; r.fwhm_km = 2.3548*best.sg;
    sL = best.mu - prm.KEDGE*best.sg; sR = best.mu + prm.KEDGE*best.sg; r.edges_in = sL >= ss(1) && sR <= ss(end);
    if A < prm.PMIN, r.reason = "amplitude below PMIN";
    elseif r.SNR < prm.SNRMIN, r.reason = "signal to noise below SNRMIN";
    elseif r.fwhm_km > prm.FWMAX, r.reason = "wider than FWMAX";
    else, r.det = true; r.reason = "plume";
    end
end

function q = edge_quantities(s, c, mu, sg, kedge)
% Background, enhancement and width with the edges at mu +/- kedge*sigma (manual convention).
    ok = ~isnan(s) & ~isnan(c); ss = s(ok); cc = c(ok);
    [~, iL] = min(abs(ss - (mu - kedge*sg))); [~, iR] = min(abs(ss - (mu + kedge*sg)));
    if iR <= iL, iR = min(iL + 1, numel(ss)); end
    B = (cc(iL) + cc(iR))/2;
    q.P = max(cc(iL:iR)) - B; q.E = mean(cc(iL:iR)) - B; q.D = ss(iR) - ss(iL); q.sL = ss(iL); q.sR = ss(iR);
end

%% plot_Figure8_altitude.m
%  Figure 8: flux deviation vs. transect height deviation for the four
%  background methods, with the cos(theta) fluxes of All_Solutions.mat and the
%  method names of the revision (Downwind Edge).  Each panel shows the per
%  transect deviations from the daily mean, a least-squares line, its slope
%  with 1 SE and the p-value of the slope (R1.10).  The script also prints the
%  correlation between height deviation and time deviation (quoted as r = 0.55
%  in Section 3.1.1).
%
%  Days with a single downwind transect contribute no deviation and are skipped;
%  Methods B to D use the 32 days with an upwind transect.

clear; close all; clc;

%% ------------------------------------------------------------ settings
SOL_FILE = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB/All_Solutions.mat';
outName  = 'Fig8_flux_vs_height';
FN = 'Arial';
names  = {'Downwind Edge', 'Upwind', 'Adjusted ABL Upwind', 'Upwind Edge'};
fields = {'SolA', 'SolB', 'SolC', 'SolD'};

%% ------------------------------------------------------------ data
assert(isfile(SOL_FILE), 'All_Solutions.mat not found at %s', SOL_FILE);
S = load(SOL_FILE);
fig = figure('Color','w','Units','pixels','Position',[80 80 1150 820]);
for m = 1:4
    Sol = S.(fields{m});
    dF = []; dZ = []; dT = [];
    for i = 1:numel(Sol)
        f = Sol(i).Flux(:); z = Sol(i).LegHeight(:); t = Sol(i).DownwindTime(:);
        ok = ~isnan(f) & ~isnan(z);
        if nnz(ok) < 2, continue; end
        dF = [dF; f(ok) - mean(f(ok))]; dZ = [dZ; z(ok) - mean(z(ok))]; dT = [dT; t(ok) - mean(t(ok))]; %#ok<AGROW>
    end
    lm = fitlm(dZ, dF);
    slope = lm.Coefficients.Estimate(2); se = lm.Coefficients.SE(2); pv = lm.Coefficients.pValue(2);
    r = corrcoef(dZ, dT); rZT = r(1,2);
    fprintf('%-20s n = %3d transects: slope %+.4f +- %.4f mol/s per m, p = %.3f; corr(height dev, time dev) = %.2f\n', ...
        names{m}, numel(dF), slope, se, pv, rZT);

    ax = subplot(2,2,m); hold(ax,'on'); box(ax,'on');
    xl = [-600 600];
    scatter(ax, dZ, dF, 28, [0 0.45 0.74], 'filled', 'MarkerEdgeColor', [0 0.25 0.5], 'LineWidth', 0.3);
    xf = linspace(xl(1), xl(2), 50);
    plot(ax, xf, lm.Coefficients.Estimate(1) + slope*xf, '-r', 'LineWidth', 1.6);
    yl = max(150, 50*ceil(max(abs(dF))/50));
    xlim(ax, xl); ylim(ax, [-yl yl]);
    grid(ax, 'on'); ax.GridAlpha = 0.12;
    xlabel(ax, 'Meters from mean transect height', 'FontSize', 12, 'FontName', FN);
    ylabel(ax, 'Flux deviation (mol s^{-1})', 'FontSize', 12, 'FontName', FN);
    title(ax, names{m}, 'FontSize', 13, 'FontName', FN);
    text(ax, 0.97, 0.95, sprintf('slope = %.3f \\pm %.3f mol s^{-1} m^{-1}\np = %.2f  (n = %d)', slope, se, pv, numel(dF)), ...
        'Units','normalized', 'HorizontalAlignment','right', 'VerticalAlignment','top', 'FontSize', 10, 'FontName', FN, ...
        'BackgroundColor','w', 'EdgeColor','k', 'Margin', 3);
    set(ax, 'FontSize', 11, 'FontName', FN, 'TickDir', 'out');
end
set(fig, 'PaperPositionMode', 'auto');
print(fig, [outName '.png'], '-dpng', '-r300');
print(fig, [outName '.eps'], '-depsc');
fprintf('Saved %s.png / .eps\n', outName);

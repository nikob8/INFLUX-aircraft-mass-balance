%% plot_Figure9_SSLF.m
%  Revised Figure 9 (old Figure 10): SSLF CH4 emission rate per flight day.
%  Reviewer 1 (R1.11): the date labels are replaced by a chronological flight
%  index (dates are listed in Table S1).
%
%  Same data and logic as Zach's ReadLandfillsDaily.m: every Landfill/<Date>_Loop<k>.mat
%  with Flag > 0 counts; the day value is the mean of its transects; the error
%  bar is std(...,1) of the transects; red = only one transect that day.
%  New: with APPLY_COS = true each transect flux is multiplied by cos(theta) of
%  its transect (from All_Solutions.mat), so the figure matches Table S1 v3.

clear; close all; clc;

%% ------------------------------------------------------------ settings
LF_DIR    = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance/Landfill';
SOL_FILE  = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB/All_Solutions.mat';   % from Collect_Solutions.m
APPLY_COS = true;
outName   = 'Fig9_SSLF_daily';

%% ------------------------------------------------------------ cos(theta) lookup
useCos = APPLY_COS;
if useCos
    assert(isfile(SOL_FILE), 'All_Solutions.mat not found at %s (fix SOL_FILE above)', SOL_FILE);
    S = load(SOL_FILE, 'SolA'); SolA = S.SolA;
end

%% ------------------------------------------------------------ read the Loop files (as ReadLandfillsDaily.m)
first = dir(fullfile(LF_DIR, '*Loop1*'));
dates = sort(unique(cellfun(@(n) n(1:8), {first.name}, 'UniformOutput', false)));
dayMean = []; daySD = []; dayN = []; dayDate = {};
for i = 1:numel(dates)
    D = dates{i};
    files = dir(fullfile(LF_DIR, [D '_Loop*.mat']));
    fl = [];
    for j = 1:numel(files)
        L = load(fullfile(files(j).folder, files(j).name));
        if L.Flag > 0
            v = L.Flux;
            if useCos
                k  = sscanf(regexprep(files(j).name, '.*_Loop', ''), '%d');   % transect number
                kk = find(strcmp({SolA.Date}, D), 1);
                assert(~isempty(kk) && numel(SolA(kk).CosTheta) >= k, 'No CosTheta for %s transect %d', D, k);
                v = v * SolA(kk).CosTheta(k);
            end
            fl(end+1) = v; %#ok<AGROW>
        end
    end
    if ~isempty(fl)
        dayMean(end+1) = mean(fl);  daySD(end+1) = std(fl, 1); %#ok<AGROW>
        dayN(end+1) = numel(fl);    dayDate{end+1} = D;         %#ok<AGROW>
    end
end
n = numel(dayMean); idx = 1:n;
good = dayN > 1; one = dayN == 1;
fprintf('%d days with a quantified plume (%d with one transect only); cos(theta) applied: %d\n', n, nnz(one), useCos);
fprintf('Mean of daily means %.1f mol/s, SD %.1f, 2SE %.1f, median %.1f\n', mean(dayMean), std(dayMean), ...
    2*std(dayMean)/sqrt(n), median(dayMean));
for i = 1:n, fprintf('  %2d  %s  %5.1f ± %4.1f  (n=%d)\n', i, dayDate{i}, dayMean(i), daySD(i), dayN(i)); end

%% ------------------------------------------------------------ figure
FN = 'Helvetica';
fig = figure('Color','w','Units','pixels','Position',[100 100 1150 500]);
ax = axes(fig); hold(ax,'on'); box(ax,'on');
mu = mean(dayMean);
h0 = plot(ax, [0 n+1], [mu mu], '--', 'Color',[0.45 0.45 0.45], 'LineWidth',1.4);
h1 = errorbar(ax, idx(good), dayMean(good), daySD(good), 'o', 'LineWidth',1.6, 'Color',[0 0 0.8], ...
    'MarkerFaceColor',[0 0 0.8], 'MarkerSize',7, 'CapSize',5);
h2 = errorbar(ax, idx(one), dayMean(one), daySD(one), 'o', 'LineWidth',1.6, 'Color',[0.8 0 0], ...
    'MarkerFaceColor',[0.8 0 0], 'MarkerSize',7, 'CapSize',5);
xlim(ax, [0 n+1]); xticks(ax, 1:n); xtickangle(ax, 0);
ylim(ax, [0 10*ceil(max(dayMean+daySD)/10 + 0.3)]);         % room above the tallest error bar
set(ax, 'FontSize',12, 'FontName',FN, 'TickDir','out', 'YGrid','on', 'GridAlpha',0.12);
xlabel(ax, 'Flight day (chronological order; dates in Table S1)', 'FontName',FN, 'FontSize',13);
ylabel(ax, 'Emission rate (mol s^{-1})', 'FontName',FN, 'FontSize',13);
lg = legend(ax, [h1 h2 h0], {'two or more transects (mean \pm 1 SD)', 'one transect', ...
    sprintf('mean of daily values (%.0f mol s^{-1})', mu)}, ...
    'Location','northoutside', 'Orientation','horizontal', 'FontSize',11, 'FontName',FN, 'Box','off');

set(fig,'PaperPositionMode','auto');
print(fig, [outName '.png'], '-dpng', '-r300');
print(fig, [outName '.eps'], '-depsc');
fprintf('Saved %s.png / .eps\n', outName);

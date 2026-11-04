%% Collect_Landfill.m
%  Gathers the South Side Landfill (SSLF) mass balance results of all flight
%  days into two tables (Section 2.4 and Section 3.1.3 of the paper):
%     Landfill_transects.csv   one row per downwind transect with a quantified
%                              landfill plume: emission rate, plume width
%                              (distance), CH4 enhancement and peak
%     Landfill_daily.csv       one row per flight day: mean and SD of the
%                              transects, number of transects
%
%  Input: Landfill/<yyyymmdd>_Loop<k>.mat, one file per downwind transect k,
%  written by MassBalanceLandfill.m, with the variables
%     Flag          -1 = no landfill plume could be identified on the transect;
%                   > 0 = plume quantified (these transects are used in the
%                   paper; the original ReadLandfills.m also kept a stricter
%                   list with Flag > 1)
%     Start, End    indices of the plume end points, chosen by eye, inside the
%                   downwind segment of the transect (MassBalanceLandfill.m
%                   recomputes the results from them without interaction)
%     Flux          emission rate of the landfill plume (mol/s), before the
%                   cos(theta) factor of Equation 1
%     Distance      horizontal extent of the plume along the transect (km)
%     Enhancement   mean CH4 enhancement of the plume above its local
%                   background (ppm; written in ppb to the csv)
%     Peak          peak CH4 enhancement of the plume (ppm; ppb in the csv)
%  and All_Solutions.mat (from Collect_Solutions.m) for cos(theta) of each
%  transect.  The cos(theta)-corrected emission rates are the landfill values
%  of Table S1 (column Flux_Landfill_mols) and of Figures 9, 10 and S3.
%
%  This script replaces ReadLandfills.m and ReadLandfillsDaily.m (same data
%  and the same selection, Flag > 0; the daily plot of ReadLandfillsDaily.m is
%  now plot_Figure9_SSLF.m).

clear; clc;

%% ------------------------------------------------------------ settings
LF_DIR   = fullfile(pwd, 'Landfill');                                                      % Loop files written by MassBalanceLandfill.m
LF_REF   = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance/Landfill';    % optional: original results to compare with ('' to skip)
SOL_FILE = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB/All_Solutions.mat';
APPLY_COS = true;                       % multiply each transect by cos(theta) of the transect (Equation 1)
OUT_DIR   = pwd;

%% ------------------------------------------------------------ read every Loop file
if APPLY_COS
    assert(isfile(SOL_FILE), 'All_Solutions.mat not found at %s (fix SOL_FILE above)', SOL_FILE);
    S = load(SOL_FILE, 'SolA'); SolA = S.SolA;
end
files = dir(fullfile(LF_DIR, '2*_Loop*.mat'));
assert(~isempty(files), 'no <yyyymmdd>_Loop<k>.mat files in %s', LF_DIR);
T = table('Size', [0 11], 'VariableTypes', {'string','double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Date','Transect','Flag','Start','End','cos_theta','Flux_raw_mols','Flux_Landfill_mols','Distance_km','Enhancement_ppb','Peak_ppb'});
for j = 1:numel(files)
    tok = regexp(files(j).name, '^(\d{8})_Loop(\d+)\.mat$', 'tokens', 'once');
    if isempty(tok), continue; end
    D = tok{1}; k = str2double(tok{2});
    L = load(fullfile(files(j).folder, files(j).name));
    if ~isfield(L, 'Flag') || L.Flag(1) <= 0, continue; end       % no quantified plume on this transect
    ct = NaN;
    if APPLY_COS
        kk = find(strcmp({SolA.Date}, D), 1);
        assert(~isempty(kk) && numel(SolA(kk).CosTheta) >= k, 'no CosTheta for %s transect %d', D, k);
        ct = SolA(kk).CosTheta(k);
    end
    fl = double(L.Flux(1)); if APPLY_COS, fl = fl*ct; end
    T(end+1,:) = {string(sprintf('%s-%s-%s', D(1:4), D(5:6), D(7:8))), k, double(L.Flag(1)), getfield_or(L,'Start'), getfield_or(L,'End'), ct, double(L.Flux(1)), fl, ...
        getfield_or(L,'Distance'), 1000*getfield_or(L,'Enhancement'), 1000*getfield_or(L,'Peak')}; %#ok<AGROW>   % Enhancement and Peak are in ppm in the Loop files
end
T = sortrows(T, {'Date','Transect'});
if ~APPLY_COS, T.cos_theta = []; T.Flux_Landfill_mols = []; T = renamevars(T, 'Flux_raw_mols', 'Flux_Landfill_mols'); end
writetable(T, fullfile(OUT_DIR, 'Landfill_transects.csv'));

%% ------------------------------------------------------------ daily table
[dates, ~, g] = unique(T.Date);
Dy = table(dates, 'VariableNames', {'Date'});
Dy.N_transects = accumarray(g, 1);
Dy.Flux_Landfill_mean_mols = accumarray(g, T.Flux_Landfill_mols, [], @mean);
Dy.Flux_Landfill_SD_mols   = accumarray(g, T.Flux_Landfill_mols, [], @(v) std(v, 1));   % population SD, as in ReadLandfillsDaily.m
Dy.Distance_mean_km        = accumarray(g, T.Distance_km, [], @(v) mean(v, 'omitnan'));
Dy.Enhancement_mean_ppb    = accumarray(g, T.Enhancement_ppb, [], @(v) mean(v, 'omitnan'));
Dy.Peak_max_ppb            = accumarray(g, T.Peak_ppb, [], @(v) max(v, [], 'omitnan'));
writetable(Dy, fullfile(OUT_DIR, 'Landfill_daily.csv'));

%% ------------------------------------------------------------ comparison with the original results (optional)
if ~isempty(LF_REF) && isfolder(LF_REF)
    dmax = 0; nref = 0; missing = {};
    if ismember('Flux_raw_mols', T.Properties.VariableNames), raw = T.Flux_raw_mols; else, raw = T.Flux_Landfill_mols; end
    for j = 1:height(T)
        D = strrep(char(T.Date(j)), '-', '');
        f = fullfile(LF_REF, sprintf('%s_Loop%d.mat', D, T.Transect(j)));
        if ~isfile(f), missing{end+1} = sprintf('%s T%d', D, T.Transect(j)); continue; end %#ok<AGROW>
        R = load(f); nref = nref + 1;
        if isfield(R, 'Flux'), dmax = max(dmax, abs(double(R.Flux(1)) - raw(j))); end
    end
    fprintf('comparison with %s: %d transects, largest difference in Flux (before cos theta) %.3f mol/s', LF_REF, nref, dmax);
    if ~isempty(missing), fprintf('; not in the reference: %s', strjoin(missing, ', ')); end
    fprintf('\n');
end

%% ------------------------------------------------------------ summary
fprintf('%d transects with a quantified landfill plume on %d flight days (cos theta applied: %d)\n', height(T), height(Dy), APPLY_COS);
fprintf('transects: mean %.1f mol/s, median %.1f, SD %.1f; Flag > 1 on %d transects\n', mean(T.Flux_Landfill_mols), median(T.Flux_Landfill_mols), std(T.Flux_Landfill_mols), nnz(T.Flag > 1));
fprintf('daily means: mean %.1f mol/s, median %.1f, SD %.1f, 2 SE %.1f (%d days with one transect only)\n', mean(Dy.Flux_Landfill_mean_mols), ...
    median(Dy.Flux_Landfill_mean_mols), std(Dy.Flux_Landfill_mean_mols), 2*std(Dy.Flux_Landfill_mean_mols)/sqrt(height(Dy)), nnz(Dy.N_transects == 1));
fprintf('wrote Landfill_transects.csv and Landfill_daily.csv in %s\n', OUT_DIR);

function v = getfield_or(S, f)
    if isfield(S, f) && ~isempty(S.(f)), v = double(S.(f)(1)); else, v = NaN; end
end

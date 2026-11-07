%% Check_PBLH_20121004.m  (optional, for R2.9)
%  Modeled boundary layer height at the Marion County centre on 10/04/2012,
%  hour by hour, from the WRF-Chem run used for the diagnostic ratio.  Uses the
%  same wrfout lookup as plot_Figure7_WRF.m.  Prints PBLH (WRF's own diagnostic) and,
%  as a second opinion, the height where the modeled potential temperature
%  first exceeds its 100 m value by 1 K.

clear; clc;
WRF_ROOT = '/storage/group/zrb5027/default/INFLUX/nvb5011/';
CASE     = '20121004'; DOMAIN = 2; HOURS = 14:21;
LATC = 39.7806; LONC = -86.1401;

yyyy = CASE(1:4); mm = CASE(5:6); dd = CASE(7:8);
fprintf('WRF boundary layer height at the Marion County centre, %s/%s/%s\n', mm, dd, yyyy);
fprintf('  UTC   PBLH (m)   z(theta - theta_100m > 1 K) (m)\n');
for hh = HOURS
    fn = find_wrfout(WRF_ROOT, DOMAIN, yyyy, mm, dd, hh);
    if isempty(fn), fprintf('  %02d    no wrfout found\n', hh); continue; end
    XLAT = double(ncread(fn,'XLAT')); XLAT = XLAT(:,:,1);
    XLON = double(ncread(fn,'XLONG')); XLON = XLON(:,:,1);
    [~, idx] = min((XLAT(:)-LATC).^2 + ((XLON(:)-LONC)*cosd(LATC)).^2);
    [i, j] = ind2sub(size(XLAT), idx);
    pblh = double(ncread(fn, 'PBLH', [i j 1], [1 1 1]));
    nk = 30;
    T  = double(ncread(fn, 'T',   [i j 1 1], [1 1 nk 1])) + 300;          % potential temperature (K)
    ph = double(ncread(fn, 'PH',  [i j 1 1], [1 1 nk+1 1])); phb = double(ncread(fn, 'PHB', [i j 1 1], [1 1 nk+1 1]));
    zw = squeeze(ph + phb)/9.81; hgt = double(ncread(fn, 'HGT', [i j 1], [1 1 1]));
    z  = 0.5*(zw(1:end-1) + zw(2:end)) - hgt;                               % m AGL at mass levels
    th = squeeze(T);
    th100 = interp1(z, th, 100);
    k = find(th - th100 > 1, 1);
    ztop = NaN; if ~isempty(k) && k > 1, ztop = interp1(th(k-1:k) - th100, z(k-1:k), 1); end
    fprintf('  %02d    %6.0f     %6.0f\n', hh, pblh, ztop);
end

function fn = find_wrfout(wrfRoot, dom, yyyy, mm, dd, hh)
    key = sprintf('%s-%s-%s_%02d', yyyy, mm, dd, hh);
    pat = sprintf('wrfout_d0%d_(\\d{4}-\\d{2}-\\d{2})_(\\d{2})', dom);
    dirs = dir(fullfile(wrfRoot,'WRF_*')); dirs = dirs([dirs.isdir]);
    fn = ''; cands = {};
    for i = 1:numel(dirs)
        L = dir(fullfile(wrfRoot, dirs(i).name, sprintf('wrfout_d0%d_*', dom)));
        for j = 1:numel(L)
            tok = regexp(L(j).name, pat, 'tokens', 'once');
            if ~isempty(tok) && strcmp([tok{1} '_' tok{2}], key)
                cands{end+1} = fullfile(wrfRoot, dirs(i).name, L(j).name); %#ok<AGROW>
            end
        end
    end
    if isempty(cands), return; end
    isdi = endsWith(cands, '_di');
    if any(isdi), fn = cands{find(isdi,1)}; else, fn = cands{1}; end
end

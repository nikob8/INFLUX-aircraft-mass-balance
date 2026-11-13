%% Table_S1_statistics.m
%  Statistics quoted in the paper that follow directly from Table S1
%  (Indy_Mass_Balance_Table_v5.csv, folder data/ of the repository; one row per
%  downwind transect) and from SSLF_plume_criteria.csv (Data Set S4).
%
%   1. Time of day of the downwind transects, in UTC and in local time (Section 2.2).
%   2. Order in which the transect altitudes were flown on each day (Sections 2.3 and 3.1.1).
%   3. South Side Landfill (SSLF) plume by wind direction (Sections 2.4 and 3.1.3, Text S3):
%      path of the plume within Marion County along the mean wind, number of transects
%      with a quantified plume in the lower and upper half of the boundary layer, plume
%      width, and daily SSLF emission rates, for northerly/easterly and for
%      southerly/westerly winds, with a Fisher exact test and a rank-sum test.
%
%  Wind sectors are the four quadrants centred on N, E, S and W (N = 315 to 45 deg),
%  from the mean (mass balance) wind direction of the flight day, as in
%  plot_FigS4_flux_vs_meteorology.py.
%  Local time: Eastern Standard Time (UTC - 5 h); Eastern Daylight Time (UTC - 4 h) from
%  the second Sunday of March to the first Sunday of November.
%
%  Base MATLAB only (no toolboxes, no readtable); also runs in GNU Octave.
%  Expected output: see the end of this file.

clear; clc;
csvFile  = 'Indy_Mass_Balance_Table_v5.csv';
critFile = 'SSLF_plume_criteria.csv';            % optional (plume widths, automated selection)
SSLF     = [39.714691162109375, -86.205429077148438];                 % lat, lon
% Marion County corners (lat, lon): bottom right, bottom left, top left, top right
County   = [39.63760423971629, -85.95148383144084; 39.63337385377353, -86.32570560518701; ...
            39.92446141560437, -86.32734193809195; 39.9271330073512,  -85.95582253907553];
KM_PER_DEG = 111.195;

%% ------------------------------------------------------------ read Table S1
fid = fopen(csvFile, 'r'); assert(fid > 0, 'cannot open %s', csvFile);
hdr = regexp(strtrim(fgetl(fid)), ',', 'split'); R = {};
while true
    l = fgetl(fid); if ~ischar(l), break; end
    l = strtrim(l); if isempty(l), continue; end
    f = regexp(l, ',', 'split'); f(end+1:numel(hdr)) = {''}; R(end+1,:) = f(1:numel(hdr)); %#ok<SAGROW>
end
fclose(fid);
col = @(name) R(:, strcmp(hdr, name));
num = @(name) str2double(col(name));
dateStr = col('Date'); n = numel(dateStr);
ymd = zeros(n,3); for i = 1:n, ymd(i,:) = sscanf(dateStr{i}, '%d-%d-%d')'; end
hm = zeros(n,2); tt = col('Time_UTC'); for i = 1:n, hm(i,:) = sscanf(tt{i}, '%d:%d')'; end
tUTC = hm(:,1)*60 + hm(:,2);                      % start time of the transect, minutes UTC
z = num('Height_AGL_m'); zi = num('ABL_height_m');
wd = num('Mean_wind_direction_deg'); FL = num('Flux_Landfill_mols');
[days, iFirst] = unique(dateStr, 'stable'); nD = numel(days);     % flight days in chronological order
[~, dayIdx] = ismember(dateStr, days);                            % flight day of each transect
fprintf('%d downwind transects on %d flight days\n\n', n, nD);

%% ------------------------------------------------------------ 1. time of day
dn = datenum(ymd(:,1), ymd(:,2), ymd(:,3)); isDST = false(n,1);
for i = 1:n
    s = datenum(ymd(i,1), 3, 8):datenum(ymd(i,1), 3, 14);  s = s(weekday(s) == 1);   % second Sunday of March
    e = datenum(ymd(i,1), 11, 1):datenum(ymd(i,1), 11, 7); e = e(weekday(e) == 1);   % first Sunday of November
    isDST(i) = dn(i) >= s && dn(i) < e;
end
tLoc = tUTC - 300 + 60*isDST;
hhmm = @(t) sprintf('%02d:%02d', floor(t/60), round(mod(t,60)));
fprintf('1. Start times of the downwind transects\n');
fprintf('   UTC:   %s to %s, median %s; %d of %d (%.0f%%) between 17 and 20 UTC\n', hhmm(min(tUTC)), hhmm(max(tUTC)), ...
    hhmm(median(tUTC)), sum(tUTC >= 17*60 & tUTC < 20*60), n, 100*mean(tUTC >= 17*60 & tUTC < 20*60));
fprintf('   local: %s to %s, median %s (%d flight days on daylight saving time, %d on standard time)\n\n', ...
    hhmm(min(tLoc)), hhmm(max(tLoc)), hhmm(median(tLoc)), sum(isDST(iFirst)), sum(~isDST(iFirst)));

%% ------------------------------------------------------------ 2. order of the transect altitudes
nUp = 0; nDown = 0; nMulti = 0; nAsc = 0; nDesc = 0; nFirstLowest = 0; rho = [];
for d = 1:nD
    k = find(dayIdx == d); [~, o] = sort(tUTC(k)); k = k(o);
    if numel(k) < 2, continue; end
    nMulti = nMulti + 1; dz = diff(z(k));
    nUp = nUp + sum(dz > 0); nDown = nDown + sum(dz < 0);
    nAsc = nAsc + all(dz > 0); nDesc = nDesc + all(dz < 0);
    nFirstLowest = nFirstLowest + (z(k(1)) == min(z(k)));
    if numel(k) >= 3
        [~, o] = sort(z(k)); rz = zeros(numel(k),1); rz(o) = (1:numel(k))';     % rank of the altitude
        c = corrcoef((1:numel(k))', rz); rho(end+1) = c(1,2); %#ok<SAGROW>
    end
end
fprintf('2. Order of the transect altitudes (%d days with two or more downwind transects)\n', nMulti);
fprintf('   consecutive pairs of transects: %d; later transect higher in %d, lower in %d\n', nUp + nDown, nUp, nDown);
fprintf('   days with strictly ascending altitudes: %d; strictly descending: %d; first transect the lowest: %d\n', nAsc, nDesc, nFirstLowest);
fprintf('   rank correlation of altitude with flight order (days with three or more transects): median %.2f, positive on %d of %d days\n\n', ...
    median(rho), sum(rho > 0), numel(rho));

%% ------------------------------------------------------------ 3. SSLF plume by wind direction
wdDay = wd(iFirst); sec = mod(floor(mod(wdDay + 45, 360)/90), 4) + 1;            % 1 N, 2 E, 3 S, 4 W
names = 'NESW';
fprintf('3. SSLF plume by wind direction\n   flight days by sector:');
for s = 1:4, fprintf(' %s %d', names(s), sum(sec == s)); end
fprintf('\n');
isNE = sec <= 2;                                                                  % northerly or easterly
% path of the plume inside the county: from the SSLF along the mean wind to the county boundary
P = [(County(:,2) - SSLF(2))*cosd(SSLF(1)), County(:,1) - SSLF(1)]*KM_PER_DEG;   % corners, km east and north of the SSLF
pathKm = NaN(nD,1);
for d = 1:nD
    u = [sind(wdDay(d) + 180), cosd(wdDay(d) + 180)];                             % direction in which the air moves
    for c = 1:4
        a = P(c,:); b = P(mod(c,4) + 1,:); e = b - a; den = u(1)*e(2) - u(2)*e(1);
        if abs(den) < 1e-12, continue; end
        t = (a(1)*e(2) - a(2)*e(1))/den; w = (a(1)*u(2) - a(2)*u(1))/den;
        if t > 0 && w >= 0 && w <= 1, pathKm(d) = min([pathKm(d) t]); end
    end
end
fprintf('   SSLF to the county boundaries: %.1f km (south), %.1f km (west), %.1f km (north), %.1f km (east)\n', ...
    -mean(P([1 2],2)), -mean(P([2 3],1)), mean(P([3 4],2)), mean(P([1 4],1)));
fprintf('   path of the plume inside the county along the mean wind: median %.1f km on the %d N/E days, %.1f km on the %d S/W days\n', ...
    median(pathKm(isNE)), sum(isNE), median(pathKm(~isNE)), sum(~isNE));

quant = ~isnan(FL); low = z./zi < 0.5; ne = isNE(dayIdx);
cnt = @(m) sprintf('%d of %d', sum(quant & m), sum(m));
fprintf('   transects with a quantified plume: N/E %s (lower half of the ABL %s, upper half %s)\n', cnt(ne), cnt(ne & low), cnt(ne & ~low));
fprintf('                                      S/W %s (lower half of the ABL %s, upper half %s)\n', cnt(~ne), cnt(~ne & low), cnt(~ne & ~low));

% Fisher exact test (two-sided: tables no more probable than the observed one)
lnC = @(nn, kk) gammaln(nn + 1) - gammaln(kk + 1) - gammaln(nn - kk + 1);
tab = {'all transects', true(n,1); 'lower half of the ABL', low; 'upper half of the ABL', ~low};
for q = 1:3
    m = tab{q,2}; a = sum(quant & ne & m); b = sum(~quant & ne & m); c = sum(quant & ~ne & m); d2 = sum(~quant & ~ne & m);
    N = a + b + c + d2; r1 = a + b; c1 = a + c; x = max(0, c1 - (N - r1)):min(r1, c1);
    px = exp(lnC(r1, x) + lnC(N - r1, c1 - x) - lnC(N, c1)); p0 = px(x == a);
    fprintf('   Fisher exact test, N/E against S/W, %s: p = %.2f\n', tab{q,1}, sum(px(px <= p0*(1 + 1e-9))));
end

% daily SSLF emission rates (mean of the quantified transects of the day)
daily = NaN(nD,1); for d = 1:nD, v = FL(dayIdx == d); v = v(~isnan(v)); if ~isempty(v), daily(d) = mean(v); end, end
ok = ~isnan(daily);
fprintf('   daily SSLF emission rate: %d days, mean %.1f +/- %.1f mol/s (2 standard errors), median %.1f mol/s\n', ...
    sum(ok), mean(daily(ok)), 2*std(daily(ok))/sqrt(sum(ok)), median(daily(ok)));
sets = {'selection of the analysts', daily};
if exist(critFile, 'file')                                                        % plume widths and automated selection
    fid = fopen(critFile, 'r'); h2 = regexp(strtrim(fgetl(fid)), ',', 'split'); R2 = {};
    while true
        l = fgetl(fid); if ~ischar(l), break; end
        l = strtrim(l); if isempty(l), continue; end
        f = regexp(l, ',', 'split'); f(end+1:numel(h2)) = {''}; R2(end+1,:) = f(1:numel(h2)); %#ok<SAGROW>
    end
    fclose(fid);
    c2 = @(name) str2double(R2(:, strcmp(h2, name)));
    assert(isequal(R2(:, strcmp(h2, 'Date')), dateStr) && isequal(c2('Transect'), num('Transect')), 'the two files list different transects');
    fw = c2('fwhm_km'); man = c2('det_manual') == 1; alg = c2('det_alg') == 1; Falg = c2('Flux_alg_mols');
    assert(isequal(man, quant), 'det_manual differs from the landfill column of Table S1');
    fprintf('   median width (FWHM) of the quantified plumes: N/E %.1f km, S/W %.1f km\n', median(fw(man & ne)), median(fw(man & ~ne)));
    dAlg = NaN(nD,1); for d = 1:nD, v = Falg(dayIdx == d & alg); v = v(~isnan(v)); if ~isempty(v), dAlg(d) = mean(v); end, end
    sets(2,:) = {'automated selection of Text S3', dAlg};
end
for q = 1:size(sets,1)
    v = sets{q,2}; A = v(isNE & ~isnan(v)); B = v(~isNE & ~isnan(v)); nA = numel(A); nB = numel(B);
    [~, o] = sort([A; B]); rk = zeros(nA + nB, 1); rk(o) = (1:nA + nB)';          % ranks (no ties in these data)
    U = sum(rk(1:nA)) - nA*(nA + 1)/2; mu = nA*nB/2;
    zs = (abs(U - mu) - 0.5)/sqrt(nA*nB*(nA + nB + 1)/12); pNorm = erfc(zs/sqrt(2));
    pExact = NaN;
    if exp(lnC(nA + nB, nA)) < 2e6                                                % exact distribution of the rank sum
        cmb = nchoosek(1:nA + nB, nA); Uall = sum(cmb, 2) - nA*(nA + 1)/2;
        pExact = mean(abs(Uall - mu) >= abs(U - mu) - 1e-9);
    end
    fprintf('   %s: N/E %.1f mol/s (%d days, %.1f to %.1f), S/W %.1f mol/s (%d days); rank-sum test p = %.2f (exact), %.2f (normal approximation)\n', ...
        sets{q,1}, mean(A), nA, min(A), max(A), mean(B), nB, pExact, pNorm);
end

%  Expected output (Indy_Mass_Balance_Table_v5.csv and SSLF_plume_criteria.csv of the repository):
%   147 downwind transects on 37 flight days
%   1. UTC 15:09 to 20:56, median 18:27; 131 of 147 (89%) between 17 and 20 UTC
%      local 11:09 to 15:56, median 14:04 (20 flight days on daylight saving time, 17 on standard time)
%   2. 33 days with two or more transects; 110 consecutive pairs, later transect higher in 86, lower in 24;
%      14 days strictly ascending, 0 strictly descending, first transect the lowest on 20;
%      median rank correlation 0.71, positive on 24 of 30 days
%   3. sectors N 5, E 1, S 9, W 22; path inside the county 10 km (N/E) and 23 km (S/W);
%      quantified plumes N/E 14 of 25 (10 of 14 and 4 of 11), S/W 72 of 122 (34 of 49 and 38 of 73);
%      Fisher p = 0.83, 1.00, 0.52; FWHM 4.1 km in both groups;
%      daily SSLF emission rates N/E 16.8 mol/s (5 days), S/W 27.6 mol/s (26 days), p = 0.06;
%      automated selection 17.3 mol/s (5 days) and 25.0 mol/s (27 days), p = 0.22

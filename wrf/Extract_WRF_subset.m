%% Extract_WRF_subset.m
%  Writes, for every flight day, the part of the WRF-Chem output that the
%  analyses of the paper use into one MATLAB file (WRF_subset_<yyyymmdd>.mat,
%  MATLAB v7.3 format, single precision), so that the model fields can be
%  released without the full wrfout files.
%
%  What each file contains
%    grid     XLAT, XLONG, HGT, COSALPHA, SINALPHA of the sub-domain (nx x ny),
%             the indices of the sub-domain in the full d02 grid, the eta levels
%             (ZNU, ZNW) and the global attributes of the run (namelist settings)
%    hours    UTC hours of the stored fields
%    U, V     earth-relative wind on mass points (m/s), nx x ny x KMAX x nHours
%    TH       potential temperature (K), same size
%    Z        height of the mass levels above sea level (m), same size
%    P        pressure (Pa), same size
%    PBLH     boundary layer height of the MYNN scheme (m), nx x ny x nHours
%    TRACER   CH4 enhancement of each of the eight tracers (ppb), stored for
%             the hours tracerHours only: nx x ny x KMAX x nTracerHours x 8
%             (1 oil and gas, 2 coal, 3 agriculture, 4 landfills, 5 other
%              regional sources, 6 to 8 Marion County).  WRF-Chem writes the
%             tracers with a 300 offset and emissions scaled by 1e6; the
%             enhancement in ppb is (tracer - 300)/1e6, as in all the analyses
%    info     text description of every variable
%
%  Sub-domain: the d02 (3 km) grid points within BOX_KM of the Marion County
%  centre (the whole d02 grid for the days listed in FULL_DAYS, which is needed
%  for the regional map of Figure 12).  Hours: from HOURS_BEFORE before the
%  first transect (the trajectory analysis of Text S2 needs 8 h) to HOURS_AFTER
%  after the last one, always including DIAG_HOURS (16 to 20 UTC, the window of
%  the diagnostic ratio).  Tracers are stored from one hour before the first
%  transect to one hour after the last one and for DIAG_HOURS.
%
%  Inputs: FlightList.mat, LegAndProfileList/<Date>.mat, INFLUX_<Date> (for
%  the transect times) and the wrfout_d02 files, found with the same lookup as
%  in the other scripts.  With BOX_KM = 160 and KMAX = 25 one day is about
%  150 MB before compression (about 600 MB for a FULL_DAYS day).

clear; clc;

%% ------------------------------------------------------------ settings
MB_ROOT  = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';
LEG_DIR  = fullfile(MB_ROOT, 'LegAndProfileList');
OBS_DIR  = '/storage/work/zrb5027/WRFGUI_V3/Observations/INFLUX';
WRF_ROOT = '/storage/group/zrb5027/default/INFLUX/nvb5011/';
OUT_DIR  = fullfile(pwd, 'WRF_subsets');
DOMAIN   = 2;          % 3 km domain
KMAX     = 25;         % lowest model levels kept (about 3 km above ground)
BOX_KM   = 160;        % half-width of the sub-domain around the county centre (as in Trajectory_check_R2_7.m)
FULL_DAYS = {'20141119'};            % days stored on the whole d02 grid (Figure 12: coal mine plumes from southern Illinois)
HOURS_BEFORE = 8;      % hours before the first transect (back trajectories)
HOURS_AFTER  = 1;      % hours after the last transect
DIAG_HOURS = 16:20;    % UTC window of the diagnostic ratio (Diagnostic_ratio.m); always stored, with tracers
LATC = 39.7806; LONC = -86.1401;     % centre of Marion County
START_AT = 1;  SKIP_DONE = true;

addpath(MB_ROOT); addpath(LEG_DIR); addpath(OBS_DIR);
if ~isfolder(OUT_DIR), mkdir(OUT_DIR); end
FL = load('FlightList'); FlightList = FL.FlightList;

for k = START_AT:numel(FlightList)
    Date = FlightList(k).name; Date = Date(1:8);
    outFile = fullfile(OUT_DIR, ['WRF_subset_' Date '.mat']);
    if SKIP_DONE && isfile(outFile), fprintf('%2d %s exists, skipped\n', k, Date); continue; end
    t0 = tic;
    try
        % ---- transect times
        L = load(fullfile(LEG_DIR, Date));
        x = load(['INFLUX_' Date]); x(x < -998) = NaN; Time = x(:,1); clear x
        idx = [];
        for i = 1:size(L.LegList,1), idx = [idx, L.LegList(i,1):L.LegList(i,2)]; end %#ok<AGROW>
        if isfield(L, 'UpwindLeg'), idx = [idx, L.UpwindLeg(1):L.UpwindLeg(2)]; end
        tLeg = Time(idx); tLeg = tLeg(~isnan(tLeg));
        hFirst = floor(min(tLeg)/3600); hLast = ceil(max(tLeg)/3600);
        hours = max(0, min(hFirst - HOURS_BEFORE, DIAG_HOURS(1))):min(23, max(hLast + HOURS_AFTER, DIAG_HOURS(end)));
        hasTracer = (hours >= hFirst - 1 & hours <= hLast + 1) | ismember(hours, DIAG_HOURS);

        % ---- files
        yyyy = Date(1:4); mm = Date(5:6); dd = Date(7:8);
        fns = arrayfun(@(h) find_wrfout(WRF_ROOT, DOMAIN, yyyy, mm, dd, h), hours, 'UniformOutput', false);
        have = ~cellfun(@isempty, fns);
        if ~any(have), error('no wrfout found for %s', Date); end
        fn0 = fns{find(have,1)};

        % ---- grid and sub-domain
        XLAT = double(ncread(fn0,'XLAT')); XLAT = XLAT(:,:,1);
        XLON = double(ncread(fn0,'XLONG')); XLON = XLON(:,:,1);
        xkm = (XLON - LONC)*111.32*cosd(LATC); ykm = (XLAT - LATC)*110.57;
        in = abs(xkm) <= BOX_KM & abs(ykm) <= BOX_KM;
        if ismember(Date, FULL_DAYS), in(:) = true; end
        [ii, jj] = find(in); i1 = min(ii); i2 = max(ii); j1 = min(jj); j2 = max(jj);
        ni = i2 - i1 + 1; nj = j2 - j1 + 1;
        grid.XLAT = single(XLAT(i1:i2, j1:j2)); grid.XLONG = single(XLON(i1:i2, j1:j2));
        hgt = double(ncread(fn0,'HGT')); grid.HGT = single(hgt(i1:i2, j1:j2, 1));
        try
            ca = double(ncread(fn0,'COSALPHA')); sa = double(ncread(fn0,'SINALPHA'));
            ca = ca(i1:i2, j1:j2, 1); sa = sa(i1:i2, j1:j2, 1);
        catch
            ca = ones(ni,nj); sa = zeros(ni,nj);
        end
        grid.COSALPHA = single(ca); grid.SINALPHA = single(sa);
        grid.i_range = [i1 i2]; grid.j_range = [j1 j2];
        grid.ZNU = double(ncread(fn0,'ZNU')); grid.ZNU = grid.ZNU(1:KMAX,1);
        grid.ZNW = double(ncread(fn0,'ZNW')); grid.ZNW = grid.ZNW(1:KMAX+1,1);
        grid.global_attributes = global_attributes(fn0);
        [~, runDir] = fileparts(fileparts(fn0)); grid.run_directory = runDir;
        grid.flight_hours = [hFirst hLast];                 % UTC hours of the first and last transect

        % ---- hourly fields
        nH = numel(hours); nF = nnz(hasTracer);
        U = NaN(ni,nj,KMAX,nH,'single'); V = U; TH = U; Z = U; P = U; PBLH = NaN(ni,nj,nH,'single');
        TRACER = NaN(ni,nj,KMAX,nF,8,'single'); fi = 0; tracerHours = hours(hasTracer);
        for q = 1:nH
            if ~have(q), if hasTracer(q), fi = fi + 1; end, continue; end
            fn = fns{q};
            u = double(ncread(fn,'U',[i1 j1 1 1],[ni+1 nj KMAX 1])); u = 0.5*(u(1:end-1,:,:) + u(2:end,:,:));
            v = double(ncread(fn,'V',[i1 j1 1 1],[ni nj+1 KMAX 1])); v = 0.5*(v(:,1:end-1,:) + v(:,2:end,:));
            U(:,:,:,q) = single(u.*ca - v.*sa); V(:,:,:,q) = single(v.*ca + u.*sa);
            TH(:,:,:,q) = single(double(ncread(fn,'T',[i1 j1 1 1],[ni nj KMAX 1])) + 300);
            zw = (double(ncread(fn,'PH',[i1 j1 1 1],[ni nj KMAX+1 1])) + double(ncread(fn,'PHB',[i1 j1 1 1],[ni nj KMAX+1 1])))/9.81;
            Z(:,:,:,q) = single(0.5*(zw(:,:,1:end-1) + zw(:,:,2:end)));
            P(:,:,:,q) = single(double(ncread(fn,'P',[i1 j1 1 1],[ni nj KMAX 1])) + double(ncread(fn,'PB',[i1 j1 1 1],[ni nj KMAX 1])));
            PBLH(:,:,q) = single(ncread(fn,'PBLH',[i1 j1 1],[ni nj 1]));
            if hasTracer(q)
                fi = fi + 1;
                for tr = 1:8        % (tracer - 300)/1e6 = enhancement in ppb (WRF-Chem offset and emission scaling)
                    TRACER(:,:,:,fi,tr) = single((double(ncread(fn, sprintf('tracer_%d',tr), [i1 j1 1 1], [ni nj KMAX 1])) - 300)/1e6);
                end
            end
        end
        info = struct();
        info.description = 'Subset of the WRF-Chem v3.6.1 output (3 km domain) used in Balashov et al. (2026), JGR Atmospheres. Fields are on mass points of the sub-domain grid.XLAT/XLONG; levels 1:KMAX from the surface up.';
        info.hours = 'UTC hours of the stored fields (one file per flight day; date in the file name); grid.flight_hours = UTC hours of the first and last aircraft transect';
        grid.z_agl_mean = squeeze(mean(mean(double(Z(:,:,:,find(have,1))) - double(grid.HGT), 1, 'omitnan'), 2, 'omitnan'))';   % mean height above ground of each level (m)
        info.U_V = 'earth-relative horizontal wind on mass points (m/s)';
        info.TH = 'potential temperature (K)'; info.Z = 'height of the mass levels above sea level (m); subtract grid.HGT for height above ground';
        info.P = 'pressure (Pa)'; info.PBLH = 'boundary layer height of the MYNN scheme (m)';
        info.TRACER = 'CH4 enhancement (ppb) of each tracer at the hours tracerHours: 1 oil and gas, 2 coal, 3 agriculture, 4 landfills, 5 other regional sources (EPA gridded inventory, Maasakkers et al. 2016), 6 to 8 Marion County (Lamb et al. 2016). Regional background enhancement = sum(TRACER(:,:,:,:,1:5),5); city enhancement = sum(TRACER(:,:,:,:,6:8),5).';
        info.grid = 'XLAT, XLONG, HGT (terrain, m), COSALPHA/SINALPHA (map rotation, already applied to U and V), i_range/j_range (indices of the sub-domain in the full 3 km grid), ZNU/ZNW (eta levels), z_agl_mean (mean height above ground of each level, m), global_attributes (namelist settings), run_directory, flight_hours';
        info.missing_hours = hours(~have);
        save(outFile, 'grid', 'hours', 'hasTracer', 'tracerHours', 'U', 'V', 'TH', 'Z', 'P', 'PBLH', 'TRACER', 'info', '-v7.3');
        d = dir(outFile);
        fprintf('%2d %s: hours %02d-%02d (%d files, %d missing), %dx%dx%d, %.0f MB, %.0f s\n', k, Date, hours(1), hours(end), nnz(have), nnz(~have), ni, nj, KMAX, d.bytes/1e6, toc(t0));
    catch err
        fprintf('%2d %s: ERROR %s (line %d)\n', k, Date, err.message, err.stack(1).line);
    end
end

function A = global_attributes(fn)
    I = ncinfo(fn); A = struct();
    for i = 1:numel(I.Attributes)
        name = matlab.lang.makeValidName(I.Attributes(i).Name);
        A.(name) = I.Attributes(i).Value;
    end
end

function fn = find_wrfout(wrfRoot, dom, yyyy, mm, dd, hh)
    persistent keys names root
    if isempty(keys) || ~strcmp(root, wrfRoot)
        root = wrfRoot; keys = {}; names = {};
        pat = sprintf('wrfout_d0%d_(\\d{4}-\\d{2}-\\d{2})_(\\d{2})', dom);
        dirs = dir(fullfile(wrfRoot,'WRF_*')); dirs = dirs([dirs.isdir]);
        for i = 1:numel(dirs)
            L = dir(fullfile(wrfRoot, dirs(i).name, sprintf('wrfout_d0%d_*', dom)));
            for j = 1:numel(L)
                tok = regexp(L(j).name, pat, 'tokens','once');
                if ~isempty(tok), keys{end+1} = [tok{1} '_' tok{2}]; names{end+1} = fullfile(wrfRoot, dirs(i).name, L(j).name); end %#ok<AGROW>
            end
        end
    end
    key = sprintf('%s-%s-%s_%02d', yyyy, mm, dd, hh);
    cands = names(strcmp(keys, key)); fn = '';
    if isempty(cands), return; end
    isdi = endsWith(cands, '_di');
    if any(isdi), fn = cands{find(isdi,1)}; else, fn = cands{1}; end
end

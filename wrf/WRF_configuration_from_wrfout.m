%% WRF_configuration_from_wrfout.m
%  Documents the WRF-Chem configuration from the model output itself.  Every
%  wrfout file carries the namelist settings as global attributes (physics
%  options, grid, time step, nesting, projection, start date), so the
%  configuration of each simulation can be recovered even without the
%  namelist.input files.  The script reads one wrfout file of each domain from
%  every run directory and writes
%     WRF_configuration.txt   the settings of the 3 km domain and of the 9 km
%                             parent domain, the eta levels, and a table of
%                             the attributes that differ between runs (start
%                             dates, spin-up), with a check that all physics
%                             and dynamics options are identical across runs
%  Run it once; the text file goes into the repository (wrf/ folder).

clear; clc;
WRF_ROOT = '/storage/group/zrb5027/default/INFLUX/nvb5011/';
OUT_FILE = fullfile(pwd, 'WRF_configuration.txt');

dirs = dir(fullfile(WRF_ROOT, 'WRF_*')); dirs = dirs([dirs.isdir]);
runs = struct('name', {}, 'd01', {}, 'd02', {}, 'nfiles', {}, 'first', {}, 'last', {});
for i = 1:numel(dirs)
    r.name = dirs(i).name; r.d01 = []; r.d02 = []; r.nfiles = 0; r.first = ''; r.last = '';
    for dom = 1:2
        L = dir(fullfile(WRF_ROOT, dirs(i).name, sprintf('wrfout_d0%d_*', dom)));
        if isempty(L), continue; end
        names = sort({L.name});
        A = attrs(fullfile(WRF_ROOT, dirs(i).name, names{1}));
        if dom == 1, r.d01 = A; else, r.d02 = A; r.nfiles = numel(names); r.first = names{1}; r.last = names{end}; end
    end
    if ~isempty(r.d02), runs(end+1) = r; end %#ok<SAGROW>
end
fprintf('%d run directories with d02 output\n', numel(runs));
assert(~isempty(runs), 'no wrfout_d02 files found under %s', WRF_ROOT);

fid = fopen(OUT_FILE, 'w');
fprintf(fid, 'WRF-Chem configuration recovered from the global attributes of the wrfout files\n');
fprintf(fid, '%s\n\n', datestr(now));
ref = runs(1);
fprintf(fid, '==== Settings common to all runs (from %s) ====\n\n', ref.name);
for dom = 1:2
    if dom == 1, A = ref.d01; else, A = ref.d02; end
    if isempty(A), continue; end
    fprintf(fid, '--- domain d0%d ---\n', dom);
    f = fieldnames(A);
    for i = 1:numel(f)
        v = A.(f{i});
        if ischar(v), fprintf(fid, '%-36s = %s\n', f{i}, v);
        elseif isnumeric(v) && numel(v) == 1, fprintf(fid, '%-36s = %g\n', f{i}, v);
        elseif isnumeric(v), fprintf(fid, '%-36s = [%s]\n', f{i}, num2str(v(:)', '%g '));
        end
    end
    fprintf(fid, '\n');
end
% eta levels of d02
fn2 = fullfile(WRF_ROOT, ref.name, ref.first);
try
    znw = ncread(fn2, 'ZNW'); znw = znw(:,1);
    fprintf(fid, '--- eta levels (ZNW, %d levels) of d02 ---\n%s\n\n', numel(znw), num2str(znw', '%.4f '));
catch, end

% attributes that differ between runs
fprintf(fid, '==== Per-run values of the attributes that differ between runs (d02) ====\n\n');
f = fieldnames(ref.d02); differing = {};
for i = 1:numel(f)
    vals = cell(numel(runs),1);
    for q = 1:numel(runs), if isfield(runs(q).d02, f{i}), vals{q} = runs(q).d02.(f{i}); else, vals{q} = '(absent)'; end, end
    same = all(cellfun(@(v) isequal(v, vals{1}), vals));
    if ~same, differing{end+1} = f{i}; end %#ok<SAGROW>
end
fprintf(fid, '%-14s %-22s %-22s %-8s %-30s %-30s\n', 'run', 'SIMULATION_START', 'START_DATE', 'nfiles', 'first d02 file', 'last d02 file');
for q = 1:numel(runs)
    A = runs(q).d02;
    fprintf(fid, '%-14s %-22s %-22s %-8d %-30s %-30s\n', runs(q).name, getf(A,'SIMULATION_START_DATE'), getf(A,'START_DATE'), runs(q).nfiles, runs(q).first, runs(q).last);
end
fprintf(fid, '\nAttributes that differ between runs: %s\n', strjoin(differing, ', '));
phys = {'MP_PHYSICS','RA_LW_PHYSICS','RA_SW_PHYSICS','SF_SFCLAY_PHYSICS','SF_SURFACE_PHYSICS','BL_PBL_PHYSICS','CU_PHYSICS','DX','DY','DT','BOTTOM_TOP_GRID_DIMENSION','WEST_EAST_GRID_DIMENSION','SOUTH_NORTH_GRID_DIMENSION','DIFF_OPT','KM_OPT','DAMP_OPT','GRID_FDDA','TRACER_PBLMIX'};
bad = intersect(differing, phys);
if isempty(bad), fprintf(fid, 'All physics, dynamics and grid options are identical across the runs.\n');
else
    fprintf(fid, 'NOTE: these physics/grid options differ between runs: %s\n', strjoin(bad, ', '));
    for i = 1:numel(bad)                                   % per-run values of each differing option (d01 and d02)
        fprintf(fid, '\n%s per run (d01 / d02):\n', bad{i});
        for q = 1:numel(runs)
            fprintf(fid, '  %-14s %s / %s\n', runs(q).name, getf(runs(q).d01, bad{i}), getf(runs(q).d02, bad{i}));
        end
    end
end
fclose(fid);
fprintf('Wrote %s\n', OUT_FILE);
type(OUT_FILE);

function A = attrs(fn)
    I = ncinfo(fn); A = struct();
    for i = 1:numel(I.Attributes)
        A.(matlab.lang.makeValidName(I.Attributes(i).Name)) = I.Attributes(i).Value;
    end
end
function s = getf(A, f)
    if isfield(A, f), s = A.(f); if ~ischar(s), s = num2str(s); end, else, s = '?'; end
end

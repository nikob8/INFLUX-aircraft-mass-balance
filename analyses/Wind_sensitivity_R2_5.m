%% Wind_sensitivity_R2_5.m
%  Reviewer 2, Equation 1: does it matter that U cos(theta) is a single daily
%  value outside the integral instead of the aircraft-measured wind along the
%  transect inside it?
%
%  For every downwind transect the Edge enhancement DX(x) is computed exactly
%  as in MassBalanceMethodA.m and integrated two ways:
%     F_mean  = n_air zi  Ubar cos(theta_bar)  INT DX dx          (what the paper does)
%     F_local = n_air zi  INT U(x) cos(theta(x)) DX dx           (wind inside the integral)
%  U(x), theta(x) from the 1 Hz aircraft wind along the transect, smoothed with
%  a 60 s moving mean; theta(x) is the angle between the local wind and the
%  normal to the transect.  Both integrals use the same trapezoid rule on the
%  along-track distance, so the ratio isolates the wind treatment.
%
%  Output: Wind_sensitivity_R2_5.csv (one row per transect) and a summary in
%  the Command Window.  Needs INFLUX_<Date>, LegAndProfileList, InterX.m,
%  All_Solutions.mat (for n_air, zi).

clear; close all; clc;

MB_ROOT  = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';
LEG_DIR  = fullfile(MB_ROOT, 'LegAndProfileList');
OBS_DIR  = '/storage/work/zrb5027/WRFGUI_V3/Observations/INFLUX';
CODE_DIR = '/storage/work/zrb5027/WRFGUI_V3/Code';
SOL_FILE = '/storage/work/nvb5011/WRF_GUI_CT_2/WRFGUI_V2.0/revision_codes/revis_oct26/MassB/All_Solutions.mat';
SMOOTH_S = 60;                                % s, smoothing of the along-track wind
addpath(MB_ROOT); addpath(LEG_DIR); addpath(OBS_DIR); addpath(CODE_DIR);
%BottomRight BottomLeft TopLeft TopRight
CountyCorners=[39.63760423971629, -85.95148383144084; 39.63337385377353, -86.32570560518701; 39.92446141560437, -86.32734193809195; 39.9271330073512, -85.95582253907553];

S = load(SOL_FILE, 'SolA'); SolA = S.SolA;
out = {};
for iDay = 1:numel(SolA)
    Date = SolA(iDay).Date;
    try
        L = load(fullfile(LEG_DIR, Date)); LegList = L.LegList;
        if isfield(L,'UpwindLeg'), UpwindLeg = L.UpwindLeg(1):L.UpwindLeg(2); AllLeg = UpwindLeg; else, AllLeg = []; end
        for k = 1:size(LegList,1), AllLeg = [AllLeg, LegList(k,1):LegList(k,2)]; end %#ok<AGROW>
        x = load(['INFLUX_',Date]); x(x<-998) = NaN;
        Lat=x(:,2); Lon=x(:,3); WSPD=x(:,11); WDIR=x(:,12);
        U=-1*WSPD.*sind(WDIR); V=-1*WSPD.*cosd(WDIR);
        CH4=movmean(x(:,6),30,'omitnan'); if isfield(L,'CH4Flagged'), CH4 = L.CH4Flagged; end
        MeanU=mean(U(AllLeg),'omitnan'); MeanV=mean(V(AllLeg),'omitnan');
        MeanWSPD=hypot(MeanU,MeanV);
        for k = 1:size(LegList,1)
            Leg = LegList(k,1):LegList(k,2);
            [bg, Downwind] = edge_background(Lat, Lon, CH4, Leg, MeanU, MeanV, CountyCorners);
            isDown = ismember(Leg, Downwind);
            % along-track axis (km) and the unit vector of the leg
            la0 = mean(Lat(Leg),'omitnan');
            xk = (Lon(Leg)-mean(Lon(Leg),'omitnan'))*111.32*cosd(la0); yk = (Lat(Leg)-mean(Lat(Leg),'omitnan'))*110.57;
            ok = ~isnan(xk) & ~isnan(yk);
            [Vv, Dd] = eig(cov([xk(ok) yk(ok)])); [~, m] = max(diag(Dd)); e = Vv(:,m)';
            s = xk*e(1) + yk*e(2);                                    % km along the leg
            % winds along the leg (smoothed), and the daily mean wind
            uL = movmean(U(Leg), SMOOTH_S, 'omitnan'); vL = movmean(V(Leg), SMOOTH_S, 'omitnan');
            spdL = hypot(uL, vL);
            cosL = abs(uL*e(2) - vL*e(1)) ./ spdL;                   % |wind . normal| / |wind|
            cosBar = abs(MeanU*e(2) - MeanV*e(1)) / MeanWSPD;
            DX = (CH4(Leg) - bg(:)) .* isDown(:);                      % enhancement, zero outside the downwind segment
            g = ~isnan(DX) & ~isnan(spdL) & ~isnan(s);
            [ss, io] = sort(s(g)); DXs = DX(g); DXs = DXs(io); wL = spdL(g).*cosL(g); wL = wL(io);
            Imean  = MeanWSPD*cosBar * trapz(ss*1000, DXs);
            Ilocal = trapz(ss*1000, wL.*DXs);
            n_air = SolA(iDay).Density(k); zi = SolA(iDay).DownwindABL(k);
            Fmean = n_air*zi*Imean/1e6*1;   % CH4 in ppm -> mol/mol handled by /1e6 as in the Method scripts
            Flocal = n_air*zi*Ilocal/1e6;
            out(end+1,:) = {Date, k, Fmean, Flocal, Flocal/Fmean, MeanWSPD, mean(spdL(isDown),'omitnan'), std(spdL(isDown),'omitnan'), ...
                            cosBar, mean(cosL(isDown),'omitnan'), SolA(iDay).Flux(k)}; %#ok<AGROW>
        end
    catch err
        fprintf('%s: ERROR %s\n', Date, err.message);
    end
end
T = cell2table(out, 'VariableNames', {'Date','Transect','Flux_meanWind','Flux_localWind','Ratio','Ubar','U_leg_mean','U_leg_SD','cos_bar','cos_leg_mean','Flux_TableS1'});
writetable(T, 'Wind_sensitivity_R2_5.csv');

%% ------------------------------------------------------------ summary
r = T.Ratio; good = isfinite(r) & abs(T.Flux_meanWind) > 5;        % ratios are meaningless for near-zero fluxes
fprintf('\n%d transects; %d with |flux| > 5 mol/s used for the ratio statistics\n', height(T), nnz(good));
fprintf('F_local / F_mean: median %.3f, IQR %.3f to %.3f, 5th to 95th percentile %.3f to %.3f\n', ...
    median(r(good)), prctile(r(good),25), prctile(r(good),75), prctile(r(good),5), prctile(r(good),95));
fprintf('|ratio - 1| > 0.10 on %d transects, > 0.20 on %d\n', nnz(abs(r(good)-1)>0.1), nnz(abs(r(good)-1)>0.2));
fprintf('check: F_mean vs Table S1 flux, median |difference| %.1f mol/s (should be small)\n', median(abs(T.Flux_meanWind - T.Flux_TableS1)));
d = unique(T.Date); dm = NaN(numel(d),1); dl = dm;
for i = 1:numel(d), q = strcmp(T.Date, d{i}); dm(i) = mean(T.Flux_meanWind(q)); dl(i) = mean(T.Flux_localWind(q)); end
fprintf('daily means: mean wind %.1f ± %.1f mol/s, local wind %.1f ± %.1f mol/s; median |change| %.1f, max %.1f\n', ...
    mean(dm), std(dm), mean(dl), std(dl), median(abs(dl-dm)), max(abs(dl-dm)));
fprintf('along-leg wind speed SD: median %.2f m/s; along-leg wind speed vs daily mean: median ratio %.3f\n', ...
    median(T.U_leg_SD,'omitnan'), median(T.U_leg_mean./T.Ubar,'omitnan'));

%% ============================================================ Edge background (Method A, verbatim)
function [BackgroundLine, Downwind] = edge_background(Lat, Lon, CH4, Leg, MeanU, MeanV, CountyCorners)
if (MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0) %#ok<*OR2,*AND2>
    PointB1=CountyCorners(3,:); PointB2=CountyCorners(1,:);
else
    PointB1=CountyCorners(2,:); PointB2=CountyCorners(4,:);
end
PointC=[PointB1(1)+MeanV, PointB1(2)+MeanU*(111/86)]; PointA=[PointB1(1)-MeanV, PointB1(2)-MeanU*(111/86)];
Line1a=[PointA(2) PointB1(2) PointC(2);PointA(1) PointB1(1) PointC(1)];
PointC=[PointB2(1)+MeanV, PointB2(2)+MeanU*(111/86)]; PointA=[PointB2(1)-MeanV, PointB2(2)-MeanU*(111/86)];
Line1b=[PointA(2) PointB2(2) PointC(2);PointA(1) PointB2(1) PointC(1)];
Line2=[Lon(Leg)';Lat(Leg)'];
Intersection1=flipud(InterX(Line1a,Line2))'; Intersection2=flipud(InterX(Line1b,Line2))';
if ((MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0))
    if isempty(Intersection1), Intersection1=[max(Lat(Leg)) min(Lon(Leg))]; end
    if isempty(Intersection2), Intersection2=[min(Lat(Leg)) max(Lon(Leg))]; end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05,     Upwind1=find(Lon(Leg)<Intersection1(2));
    elseif abs(Intersection1(2)-Intersection2(2)) < 0.05, Upwind1=find(Lat(Leg)>Intersection1(1));
    else,                                                 Upwind1=find(Lat(Leg)>Intersection1(1) | Lon(Leg)<Intersection1(2)); end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05,     Upwind2=find(Lon(Leg)>Intersection2(2));
    elseif abs(Intersection2(2)-Intersection2(2)) < 0.05, Upwind2=find(Lat(Leg)<Intersection2(1));
    else,                                                 Upwind2=find(Lat(Leg)<Intersection2(1) | Lon(Leg)>Intersection2(2)); end
else
    if isempty(Intersection1), Intersection1=[min(Lat(Leg)) min(Lon(Leg))]; end
    if isempty(Intersection2), Intersection2=[max(Lat(Leg)) max(Lon(Leg))]; end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05,     Upwind1=find(Lon(Leg)<Intersection1(2));
    elseif abs(Intersection1(2)-Intersection2(2)) < 0.05, Upwind1=find(Lat(Leg)<Intersection1(1));
    else,                                                 Upwind1=find(Lat(Leg)<Intersection1(1) | Lon(Leg)<Intersection1(2)); end
    if abs(Intersection1(1)-Intersection2(1)) < 0.05,     Upwind2=find(Lon(Leg)>Intersection2(2));
    elseif abs(Intersection2(2)-Intersection2(2)) < 0.05, Upwind2=find(Lat(Leg)>Intersection2(1));
    else,                                                 Upwind2=find(Lat(Leg)>Intersection2(1) | Lon(Leg)>Intersection2(2)); end
end
Upwind1=Leg(Upwind1); Upwind2=Leg(Upwind2);
Downwind=Leg(~ismember(Leg,[Upwind1 Upwind2]));
Upwind1Background=median(CH4(Upwind1),'omitnan'); Upwind2Background=median(CH4(Upwind2),'omitnan');
mTemp=round(length(Leg)/10);
if isempty(Upwind1)
    if Leg(1)==Upwind2(1), Upwind1Background=median(CH4(Leg(end-mTemp:end)),'omitnan'); else, Upwind1Background=median(CH4(Leg(1:mTemp)),'omitnan'); end
elseif isempty(Upwind2)
    if Leg(1)==Upwind1(1), Upwind2Background=median(CH4(Leg(end-mTemp:end)),'omitnan'); else, Upwind2Background=median(CH4(Leg(1:mTemp)),'omitnan'); end
end
LatDistance=distance(min(Lat(Leg)),mean(Lon(Leg),'omitnan'),max(Lat(Leg)),mean(Lon(Leg),'omitnan'));
LonDistance=distance(mean(Lat(Leg),'omitnan'),min(Lon(Leg)),mean(Lat(Leg),'omitnan'),max(Lon(Leg)));
if LonDistance>LatDistance
    if Lon(Leg(1))>Lon(Leg(end)), [Upwind1Background,Upwind2Background]=deal(Upwind2Background,Upwind1Background); end
    LineCoeff=polyfit([Lon(Leg(1)),Lon(Leg(end))],[Upwind1Background,Upwind2Background],1);
    BackgroundLine=LineCoeff(1).*Lon(Leg)+LineCoeff(2);
else
    if (MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0), [Upwind1Background,Upwind2Background]=deal(Upwind2Background,Upwind1Background); end
    if Lat(Leg(1))>Lat(Leg(end)), [Upwind1Background,Upwind2Background]=deal(Upwind2Background,Upwind1Background); end
    LineCoeff=polyfit([Lat(Leg(1)),Lat(Leg(end))],[Upwind1Background,Upwind2Background],1);
    BackgroundLine=LineCoeff(1).*Lat(Leg)+LineCoeff(2);
end
end

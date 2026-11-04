%% Run_MethodsCD.m
%  Methods C (Adjusted ABL Upwind) and D (Upwind Edge) for all flight days in
%  one run.  The calculations are copied line by line from MassBalanceMethodC.m
%  and MassBalanceMethodD.m (including cos(theta) and their small differences:
%  C fills the bottom of the profile with mean() and averages the ABL with
%  nanmean(); D does the opposite).  What changed:
%
%    * each flight day / method runs inside its own function, so all memory is
%      released when it returns (no clearvars, nothing piles up between days)
%    * MAKE_FIGS = false skips all plotting (the hidden full-screen map figures
%      were the memory hog).  The fluxes do not depend on the figures.  Set it
%      to true if you also want the T#-MethodC/D.png pages; they are closed
%      after saving.
%    * a day that throws an error is logged and skipped instead of stopping the
%      run (Run_MethodsCD_log.txt)
%    * all paths are set below; results go to <pwd>/<Date>/Solution-MethodC.mat
%      and Solution-MethodD.mat with the same variable names as before, so
%      Collect_Solutions.m works unchanged
%
%  Run from your MassB folder.  Afterwards run Collect_Solutions.m.

clear; close all force; clc;

%% ------------------------------------------------------------ settings
METHODS   = {'C','D'};       % which methods to run
START_AT  = 1;               % first flight index in FlightList (resume)
SKIP_DONE = false;           % true: skip days that already have Solution-MethodX.mat
MAKE_FIGS = false;           % true: also save profile and transect PNGs

MB_ROOT  = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';
LEG_DIR  = fullfile(MB_ROOT, 'LegAndProfileList');
OBS_DIR  = '/storage/work/zrb5027/WRFGUI_V3/Observations/INFLUX';      % INFLUX_<Date>
CODE_DIR = '/storage/work/zrb5027/WRFGUI_V3/Code';                     % InterX.m, PlotINFLUX.m
OUT_ROOT = pwd;

addpath(MB_ROOT); addpath(LEG_DIR); addpath(OBS_DIR); addpath(CODE_DIR);
need = {'FlightList.mat','ElevationData.nc','InterX'};
for q = 1:numel(need)
    assert(~isempty(which(need{q})), 'Cannot find %s on the path. Add its folder in the settings block.', need{q});
end
if MAKE_FIGS && isempty(which('PlotINFLUX')), warning('PlotINFLUX.m not found: maps will have no roads.'); end

%% ------------------------------------------------------------ shared inputs (read once)
FL = load('FlightList'); FlightList = FL.FlightList;
E.Elev = ncread('ElevationData.nc','data'); E.Elev(E.Elev<0) = 0;
E.Lat  = ncread('ElevationData.nc','lat');
E.Lon  = ncread('ElevationData.nc','lon'); E.Lon(E.Lon>180) = E.Lon(E.Lon>180) - 360;

P = struct('LEG_DIR',LEG_DIR, 'OUT_ROOT',OUT_ROOT, 'MAKE_FIGS',MAKE_FIGS);
logFile = fullfile(OUT_ROOT, 'Run_MethodsCD_log.txt');
fid = fopen(logFile, 'a'); fprintf(fid, '\n==== %s ====\n', datestr(now)); fclose(fid);

%% ------------------------------------------------------------ loop
for k = START_AT:numel(FlightList)
    Date = FlightList(k).name; Date = Date(1:8);
    if ~isfolder(fullfile(OUT_ROOT,Date)), mkdir(fullfile(OUT_ROOT,Date)); end
    for m = 1:numel(METHODS)
        M = METHODS{m};
        outFile = fullfile(OUT_ROOT, Date, ['Solution-Method' M '.mat']);
        if SKIP_DONE && isfile(outFile), fprintf('%2d %s %s: exists, skipped\n', k, Date, M); continue; end
        try
            Sol = solve_flight(Date, M, P, E);
            if isempty(Sol)
                msg = 'no upwind transect in LegAndProfileList: nothing to compute (as in the original script)';
            else
                save(outFile, '-struct', 'Sol');
                msg = sprintf('flux %s mol/s, CosTheta %.3f-%.3f', mat2str(round(Sol.Flux)), min(Sol.CosTheta), max(Sol.CosTheta));
            end
        catch err
            msg = ['ERROR: ' err.message];
            if ~isempty(err.stack), msg = sprintf('%s (line %d of %s)', msg, err.stack(1).line, err.stack(1).name); end
        end
        txt = sprintf('%2d %s Method %s: %s', k, Date, M, msg);
        disp(txt);
        fid = fopen(logFile, 'a'); fprintf(fid, '%s\n', txt); fclose(fid);
        close all force;
    end
end
fprintf('Done. Log: %s\nNext: run Collect_Solutions.m\n', logFile);


%% ============================================================ one flight, one method
function Sol = solve_flight(Date, Method, P, E)
Sol = [];
L = load(fullfile(P.LEG_DIR, Date));
if ~isfield(L,'UpwindLeg'), return; end             % original: whole calculation inside "if exist('UpwindLeg')"
LegList = L.LegList; ProfileList = L.ProfileList;
UpwindLeg = L.UpwindLeg(1):L.UpwindLeg(2);
%BottomRight BottomLeft TopLeft TopRight
CountyCorners=[39.63760423971629, -85.95148383144084; 39.63337385377353, -86.32570560518701; 39.92446141560437, -86.32734193809195; 39.9271330073512, -85.95582253907553];

AllLeg=UpwindLeg;
for LegLoop=1:size(LegList,1)
    AllLeg=[AllLeg,LegList(LegLoop,1):LegList(LegLoop,2)]; %#ok<AGROW>
end

x=load(['INFLUX_',Date]);
x(x<-998)=NaN;
Time=x(:,1); Lat=x(:,2); Lon=x(:,3); Alt=x(:,4);

% terrain height (same nearest-within-0.2-deg lookup as the original loop)
Elevation=NaN(length(Alt),1);
for i=1:length(Alt)
    Row=find(abs(E.Lon-Lon(i))<0.2);
    Col=find(abs(E.Lat-Lat(i))<0.2);
    if ~isempty(Row) && ~isempty(Col)
        Elevation(i,1)=E.Elev(Row(1),Col(1));
    end
end
AltAGL=Alt-Elevation;

WSPD=x(:,11); WDIR=x(:,12);
U=-1*WSPD.*sind(WDIR);
V=-1*WSPD.*cosd(WDIR);
CH4=movmean(x(:,6),30,'omitnan');
if isfield(L,'CH4Flagged'), CH4=L.CH4Flagged; end
T=x(:,10); P_=x(:,9);
clear x
MeanU=nanmean(U(AllLeg));
MeanV=nanmean(V(AllLeg));
MeanWSPD=(MeanU.^2+MeanV.^2).^.5;
MeanWDIR=atan2d(MeanU,MeanV)+180;

%% ---- profiles: ABL height, mean ABL T and p, CH4 profile
nProf=size(ProfileList,1);
ProfileTime=NaN(nProf,1); ABLHeight=NaN(nProf,1); ABLTheta=NaN(nProf,1); ABLPressure=NaN(nProf,1); ABLT=NaN(nProf,1);
for ProfileLoop=1:nProf
    Profile=ProfileList(ProfileLoop,1):ProfileList(ProfileLoop,2);
    ProfileTime(ProfileLoop,1)=nanmean(Time(Profile));
    Heights=-25:50:5025;
    AltAGLX=AltAGL(Profile);
    Theta=T(Profile).*(1000./P_(Profile)).^(2/7);
    CH4X=CH4(Profile);
    PressureX=P_(Profile);
    CH4Averaged=NaN(1,length(Heights)-1); ThetaAveraged=CH4Averaged; PressureAveraged=CH4Averaged;
    for HeightLoop=1:length(Heights)-1
        GoodSpot=find(AltAGLX>=Heights(HeightLoop) & AltAGLX<Heights(HeightLoop+1));
        CH4Averaged(HeightLoop)=nanmean(CH4X(GoodSpot));
        ThetaAveraged(HeightLoop)=nanmean(Theta(GoodSpot));
        PressureAveraged(HeightLoop)=nanmean(PressureX(GoodSpot));
    end
    Heights=0:50:5000;
    ThetaDerivative=NaN(1,length(Heights)-2); CH4Derivative=ThetaDerivative;
    for Slope=1:length(Heights)-2
        ThetaDerivative(Slope)=(ThetaAveraged(Slope)-ThetaAveraged(Slope+2))/100;
        CH4Derivative(Slope)=(CH4Averaged(Slope)-CH4Averaged(Slope+2))/100;
    end
    ABLTop=find(CH4Derivative>2e-5 & ThetaDerivative<-5e-3);
    if isfield(L,'ManualABLHeights')
        ABLTopSpot=find(Heights==L.ManualABLHeights(ProfileLoop));
    else
        ABLTopSpot=ABLTop(1)+1;
    end
    ABLTop=Heights(ABLTopSpot);
    ABLHeight(ProfileLoop,1)=ABLTop;
    ABLTheta(ProfileLoop,1)=nanmean(ThetaAveraged(1:ABLTopSpot));
    ABLPressure(ProfileLoop,1)=nanmean(PressureAveraged(1:ABLTopSpot));
    ABLT(ProfileLoop,1)=ABLTheta(ProfileLoop)./((1000./ABLPressure(ProfileLoop)).^(2/7));
    if ProfileLoop==1
        CH4ProfileOG=CH4Averaged;
        FirstPoint=min(find(isnan(CH4ProfileOG)==0)); %#ok<MXFND>
        if FirstPoint~=1
            if strcmp(Method,'C')
                CH4ProfileOG(1:FirstPoint-1)=mean(CH4ProfileOG(FirstPoint:FirstPoint+3));
            else
                CH4ProfileOG(1:FirstPoint-1)=nanmean(CH4ProfileOG(FirstPoint:FirstPoint+3));
            end
        end
        LastPoint=max(find(isnan(CH4ProfileOG)==0)); %#ok<MXFND>
        if LastPoint~=length(Heights)
            CH4ProfileOG(LastPoint+1:end)=CH4ProfileOG(LastPoint);
        end
        Heights2=0:10:5000;
        CH4ProfileOG=interp1(Heights,CH4ProfileOG,Heights2);
    end
    if P.MAKE_FIGS
        plot_profile(Date, ProfileLoop, ProfileTime(ProfileLoop), CH4Averaged, ThetaAveraged, Heights, ABLTop, P.OUT_ROOT);
    end
end

%% ---- ABL height, T, p and density at each downwind transect
nLeg=size(LegList,1);
DownwindTime=NaN(1,nLeg);
for LegLoop=1:nLeg
    Leg=LegList(LegLoop,1):LegList(LegLoop,2);
    DownwindTime(LegLoop)=nanmean(Time(Leg));
end
UpwindHeight=nanmean(AltAGL(UpwindLeg));
if length(ProfileTime)==1
    if ProfileTime<min(DownwindTime)
        ProfileTime(2,1)=max(DownwindTime);
        ABLHeight(2,1)=ABLHeight*1.2;
    elseif ProfileTime>max(DownwindTime)
        error('Lacking early morning profile')
    end
    ABLT(2,1)=ABLT;
    ABLPressure(2,1)=ABLPressure;
end
DownwindTime(DownwindTime>max(ProfileTime))=max(ProfileTime);
DownwindTime(DownwindTime<min(ProfileTime))=min(ProfileTime);
ABLT(isnan(ABLT)==1)=nanmean(ABLT);
ABLPressure(isnan(ABLPressure)==1)=nanmean(ABLPressure);
DownwindABL=interp1(ProfileTime,ABLHeight,DownwindTime);
DownwindT=interp1(ProfileTime,ABLT,DownwindTime);
DownwindPressure=interp1(ProfileTime,ABLPressure,DownwindTime);
Density=DownwindPressure*100./8.314./DownwindT;

%% ---- ABL dilution factor (Methods C and D)
CorrectionFactor=NaN(1,length(DownwindABL));
for DownwindLoop=1:length(DownwindABL)
    Temp=abs(Heights2-DownwindABL(DownwindLoop));
    GoodSpot=find(Temp==min(Temp));
    if strcmp(Method,'C')
        MeanABLCH4=nanmean(CH4ProfileOG(1:GoodSpot));
    else
        MeanABLCH4=mean(CH4ProfileOG(1:GoodSpot));
    end
    Temp=abs(Heights2-UpwindHeight);
    GoodSpot=find(Temp==min(Temp));
    ReasonableUpwindCH4=nanmean(CH4ProfileOG(GoodSpot-7:GoodSpot+7));
    CorrectionFactor(DownwindLoop)=MeanABLCH4/ReasonableUpwindCH4;
end

%% ---- county-corner wind lines
if (MeanU>0 && MeanV>0) || (MeanU<0 && MeanV<0)
    PointB1=CountyCorners(3,:); PointB2=CountyCorners(1,:);
else
    PointB1=CountyCorners(2,:); PointB2=CountyCorners(4,:);
end
PointC=[PointB1(1)+MeanV, PointB1(2)+MeanU*(111/86)];
PointA=[PointB1(1)-MeanV, PointB1(2)-MeanU*(111/86)];
Line1a=[PointA(2) PointB1(2) PointC(2);PointA(1) PointB1(1) PointC(1)];
PointC=[PointB2(1)+MeanV, PointB2(2)+MeanU*(111/86)];
PointA=[PointB2(1)-MeanV, PointB2(2)-MeanU*(111/86)];
Line1b=[PointA(2) PointB2(2) PointC(2);PointA(1) PointB2(1) PointC(1)];

LegHeight=NaN(1,nLeg); CosTheta=NaN(1,nLeg); CorrectionFactor2=NaN(1,nLeg);
UpwindMatrix=NaN(nLeg,1); DownwindMatrix=NaN(nLeg,1); Distance=NaN(nLeg,1);
Enhancement=NaN(1,nLeg); Flux=NaN(1,nLeg);
UpwindLine=[Lon(UpwindLeg), Lat(UpwindLeg)];

for LegLoop=1:nLeg
    Leg=LegList(LegLoop,1):LegList(LegLoop,2);
    LegHeight(LegLoop)=nanmean(AltAGL(Leg));
    % cos(theta): angle between the mean wind and the normal to this leg (Eq. 1)
    LatLeg=Lat(Leg); LonLeg=Lon(Leg); okLeg=~isnan(LatLeg) & ~isnan(LonLeg);
    xKm=(LonLeg(okLeg)-mean(LonLeg(okLeg)))*111.32*cosd(mean(LatLeg(okLeg)));
    yKm=(LatLeg(okLeg)-mean(LatLeg(okLeg)))*110.57;
    [EigV,EigD]=eig(cov([xKm yKm])); [~,mAx]=max(diag(EigD)); eLeg=EigV(:,mAx)';
    uFlow=[MeanU MeanV]/MeanWSPD;
    CosTheta(LegLoop)=sqrt(max(0,1-(uFlow*eLeg')^2));

    % ---- edges vs. directly downwind (verbatim)
    Line2=[Lon(Leg)';Lat(Leg)'];
    Intersection1=flipud(InterX(Line1a,Line2))';
    Intersection2=flipud(InterX(Line1b,Line2))';
    if ((MeanU>0 && MeanV>0) || (MeanU<0 && MeanV<0))
        if isempty(Intersection1), Intersection1=[max(Lat(Leg)) min(Lon(Leg))]; end
        if isempty(Intersection2), Intersection2=[min(Lat(Leg)) max(Lon(Leg))]; end
        if abs(Intersection1(1)-Intersection2(1)) < 0.01
            Upwind1=find(Lon(Leg)<Intersection1(2));
        elseif abs(Intersection1(2)-Intersection2(2)) < 0.01
            Upwind1=find(Lat(Leg)>Intersection1(1));
        else
            Upwind1=find(Lat(Leg)>Intersection1(1) | Lon(Leg)<Intersection1(2));
        end
        if abs(Intersection1(1)-Intersection2(1)) < 0.01
            Upwind2=find(Lon(Leg)>Intersection2(2));
        elseif abs(Intersection2(2)-Intersection2(2)) < 0.01        % (sic, as in the original)
            Upwind2=find(Lat(Leg)<Intersection2(1));
        else
            Upwind2=find(Lat(Leg)<Intersection2(1) | Lon(Leg)>Intersection2(2));
        end
    else
        if isempty(Intersection1), Intersection1=[min(Lat(Leg)) min(Lon(Leg))]; end
        if isempty(Intersection2), Intersection2=[max(Lat(Leg)) max(Lon(Leg))]; end
        if abs(Intersection1(1)-Intersection2(1)) < 0.01
            Upwind1=find(Lon(Leg)<Intersection1(2));
        elseif abs(Intersection1(2)-Intersection2(2)) < 0.01
            Upwind1=find(Lat(Leg)<Intersection1(1));
        else
            Upwind1=find(Lat(Leg)<Intersection1(1) | Lon(Leg)<Intersection1(2));
        end
        if abs(Intersection1(1)-Intersection2(1)) < 0.01
            Upwind2=find(Lon(Leg)>Intersection2(2));
        elseif abs(Intersection2(2)-Intersection2(2)) < 0.01        % (sic)
            Upwind2=find(Lat(Leg)>Intersection2(1));
        else
            Upwind2=find(Lat(Leg)>Intersection2(1) | Lon(Leg)>Intersection2(2));
        end
    end
    Intersection=[Intersection1;Intersection2];
    Upwind1=Leg(Upwind1);
    Upwind2=Leg(Upwind2);
    Downwind=Leg(~ismember(Leg,[Upwind1 Upwind2]));

    % ---- upwind background traced along the wind, times the dilution factor
    BackgroundLine=NaN(1,length(Leg));
    for UpwindLoop=1:length(Leg)
        PointB=[Lat(Leg(UpwindLoop)), Lon(Leg(UpwindLoop))];
        PointC=[PointB(1)+MeanV, PointB(2)+MeanU*(111/86)];
        PointA=[PointB(1)-MeanV, PointB(2)-MeanU*(111/86)];
        Line1A=[PointA(2) PointB(2) PointC(2);PointA(1) PointB(1) PointC(1)];
        UpwindIntersection=flipud(InterX(Line1A,UpwindLine'))';
        if ~isempty(UpwindIntersection)
            TempDistance=distance(UpwindIntersection(1),UpwindIntersection(2),UpwindLine(:,2),UpwindLine(:,1));
            GoodSpot=find(TempDistance==min(TempDistance));
            UpwindAligned=UpwindLeg(GoodSpot);
            BackgroundLine(UpwindLoop)=CH4(UpwindAligned)*CorrectionFactor(LegLoop);
        end
    end
    BackgroundLine=fill_gaps(BackgroundLine);
    DownwindLine=CH4(Leg);
    R=NaN(1,101); ShiftCounter=0;
    for Shift=-50:50
        ShiftCounter=ShiftCounter+1;
        TempBackground=NaN(length(BackgroundLine)+100,1);
        TempDownwind=NaN(length(BackgroundLine)+100,1);
        TempDownwind(51:end-50)=DownwindLine;
        TempBackground(51+Shift:end-50+Shift)=BackgroundLine;
        TempCorr=corrcoef(TempBackground,TempDownwind,'Rows','complete');
        R(ShiftCounter)=TempCorr(2);
    end
    Shift=find(R==max(R),1)-51;
    BackgroundLineDoppel=NaN(length(BackgroundLine)+100,1);
    BackgroundLineDoppel(51+Shift:end-50+Shift)=BackgroundLine;
    BackgroundLine=BackgroundLineDoppel(51:end-50)';
    BackgroundGoodSpot=find(~isnan(BackgroundLine));
    BackgroundLine=fill_gaps(BackgroundLine);
    BackgroundLine(isnan(BackgroundLine))=nanmedian(BackgroundLine);

    % ---- Method D: scale so the low (10th percentile) values match
    if strcmp(Method,'D')
        BadSpot=find(isnan(DownwindLine)==0 & isnan(BackgroundLine')==0);
        TempSpot=round(length(BadSpot)*0.1);
        TempSort1=sort(BackgroundLine(BadSpot));
        TempSort2=sort(DownwindLine(BadSpot));
        CorrectionFactor2(LegLoop)=TempSort2(TempSpot)/TempSort1(TempSpot);
        BackgroundLine=BackgroundLine*CorrectionFactor2(LegLoop);
    end

    % ---- mean enhancement over the directly-downwind segment (verbatim)
    LatDistance=distance(nanmin(Lat(Leg)),nanmean(Lon(Leg)),nanmax(Lat(Leg)),nanmean(Lon(Leg)));
    LonDistance=distance(nanmean(Lat(Leg)),nanmin(Lon(Leg)),nanmean(Lat(Leg)),nanmax(Lon(Leg)));
    CH4Downwind=CH4(Downwind);
    BackgroundLineDownwind=BackgroundLine(ismember(Leg,Downwind))';
    if LonDistance>LatDistance
        xx=Lon(Downwind); col=2;
    else
        xx=Lat(Downwind); col=1;
    end
    y1=CH4Downwind;
    y1((y1-BackgroundLineDownwind)>0)=BackgroundLineDownwind((y1-BackgroundLineDownwind)>0);
    GoodSpot=find(~isnan(y1));
    Intp=min(Intersection(:,col)):0.001:max(Intersection(:,col));
    TempMatrix=[xx(GoodSpot),CH4Downwind(GoodSpot), BackgroundLineDownwind(GoodSpot)];
    [~,UniqueRows]=unique(TempMatrix(:,1));
    TempMatrix=sortrows(TempMatrix(UniqueRows,:),1);
    DownwindCH4Intp=nanmean(interp1(TempMatrix(:,1),TempMatrix(:,2),Intp));
    UpwindCH4Intp=nanmean(interp1(TempMatrix(:,1),TempMatrix(:,3),Intp));

    UpwindMatrix(LegLoop,1)=UpwindCH4Intp;
    DownwindMatrix(LegLoop,1)=DownwindCH4Intp;
    Distance(LegLoop,1)=distance(Intersection(1,:),Intersection(2,:))*111;
    Enhancement(LegLoop)=DownwindCH4Intp-UpwindCH4Intp;
    Flux(LegLoop)=Density(LegLoop).*MeanWSPD.*CosTheta(LegLoop).*(Enhancement(LegLoop)/1e6).*(Distance(LegLoop)*1000).*DownwindABL(LegLoop);

    if P.MAKE_FIGS
        cfMap=CorrectionFactor(LegLoop); if strcmp(Method,'D'), cfMap=cfMap*CorrectionFactor2(LegLoop); end
        plot_transect(Date, Method, LegLoop, DownwindTime(LegLoop), Lat, Lon, CH4, Leg, Upwind1, Upwind2, Downwind, UpwindLeg, ...
            cfMap, CountyCorners, PointB1, PointB2, Intersection, xx, CH4Downwind, BackgroundLineDownwind, ...
            BackgroundLine, BackgroundGoodSpot, LonDistance>LatDistance, Flux(LegLoop), P.OUT_ROOT);
    end
end

Sol = struct('Date',Date, 'UpwindMatrix',UpwindMatrix, 'DownwindMatrix',DownwindMatrix, ...
    'CorrectionFactor',CorrectionFactor, 'DownwindABL',DownwindABL, 'LegHeight',LegHeight, ...
    'Distance',Distance, 'Enhancement',Enhancement, 'DownwindTime',DownwindTime, 'ABLHeight',ABLHeight, ...
    'ABLTheta',ABLTheta, 'ABLT',ABLT, 'ABLPressure',ABLPressure, 'ProfileTime',ProfileTime, ...
    'MeanWSPD',MeanWSPD, 'MeanWDIR',MeanWDIR, 'Density',Density, 'Flux',Flux, 'CosTheta',CosTheta);
if strcmp(Method,'D'), Sol.CorrectionFactor2 = CorrectionFactor2; end
end

%% ============================================================ helpers
function y = fill_gaps(y)
% three-pass moving-mean gap filling, as in the original scripts
    for w = [120 240 360]
        Temp = movmean(y, w, 'omitnan');
        Bad = isnan(y);
        y(Bad) = Temp(Bad);
    end
end

function plot_profile(Date, ProfileLoop, tP, CH4Averaged, ThetaAveraged, Heights, ABLTop, OUT_ROOT)
    f=figure('Color',[1 1 1],'visible','off');
    t=tiledlayout(f,1,1);
    ax1=axes(t); plot(ax1,CH4Averaged,Heights,'-k'); xlabel(ax1,'CH_4 (ppm)');
    ax2=axes(t); plot(ax2,ThetaAveraged,Heights,'-r'); hold(ax2,'on');
    ax2.XAxisLocation='top'; ax2.Color='none'; ax1.Box='off'; ax2.Box='off'; ax2.XColor='r';
    xlabel(ax2,'Pot Temp (K)');
    plot(ax2,ThetaAveraged,ABLTop*ones(size(ThetaAveraged)),'-k','LineWidth',3);
    title(ax2,sprintf('%s/%s/%s  Profile %d  %02d:%02dZ   ABL Height = %g', Date(5:6),Date(7:8),Date(1:4), ProfileLoop, floor(tP/3600), round(mod(tP,3600)/60), ABLTop));
    saveas(f, fullfile(OUT_ROOT, Date, sprintf('Profile%d.png', ProfileLoop)));
    close(f);
end

function plot_transect(Date, Method, LegLoop, tD, Lat, Lon, CH4, Leg, Upwind1, Upwind2, Downwind, UpwindLeg, ...
        cfMap, CountyCorners, PointB1, PointB2, Intersection, xx, CH4Downwind, BgDown, BackgroundLine, BgGood, isLon, FluxVal, OUT_ROOT)
    f=figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1],'visible','off');
    subplot(1,2,1)
    worldmap([39.2 40.6],[-87 -85.2]); setm(gca,'FontSize',14); hold on
    scatterm(Lat(Upwind1),Lon(Upwind1),200,[0 0 0],'filled')
    scatterm(Lat(Downwind),Lon(Downwind),200,[1 0 0],'filled')
    scatterm(Lat(Upwind2),Lon(Upwind2),200,[0 0 0],'filled')
    scatterm(Lat(UpwindLeg),Lon(UpwindLeg),200,[0 0 0],'filled')
    scatterm(Lat(UpwindLeg),Lon(UpwindLeg),50,CH4(UpwindLeg)*cfMap,'filled')
    scatterm(Lat(Leg),Lon(Leg),50,CH4(Leg),'filled')
    if ~isempty(which('PlotINFLUX')), PlotINFLUX; end
    plotm([PointB1(1),Intersection(1,1)],[PointB1(2),Intersection(1,2)],'LineWidth',2,'Color',[0 0 0])
    plotm([PointB2(1),Intersection(2,1)],[PointB2(2),Intersection(2,2)],'LineWidth',2,'Color',[0 0 0])
    for k=1:4, k2=mod(k,4)+1; plotm(CountyCorners([k k2],1),CountyCorners([k k2],2),'LineWidth',4,'Color',[0 0 0]); end
    subplot(1,2,2); hold on
    y1=CH4Downwind; b=(y1-BgDown)<0; y1(b)=BgDown(b); g=~isnan(y1);
    patch([xx(g)' flip(xx(g))'],[BgDown(g)' flip(y1(g))'],'g')
    y1=CH4Downwind; b=(y1-BgDown)>0; y1(b)=BgDown(b); g=~isnan(y1);
    patch([xx(g)' flip(xx(g))'],[BgDown(g)' flip(y1(g))'],'r')
    if isLon, xa=Lon(Leg); xlabel('Longitude'); else, xa=Lat(Leg); xlabel('Latitude'); end
    plot(xa,CH4(Leg),'Color',[0 0 0],'LineWidth',2)
    plot(xa(BgGood),BackgroundLine(BgGood),'Color','blue','LineWidth',2)
    ylabel('CH_4 (ppm)'); set(gca,'FontSize',14)
    sgtitle(sprintf('%s/%s/%s  %02d:%02dZ  (T%d) - Method %s   Flux = %d mol/s', Date(5:6),Date(7:8),Date(1:4), ...
        floor(tD/3600), round(mod(tD,3600)/60), LegLoop, Method, round(FluxVal)),'FontSize',18)
    saveas(f, fullfile(OUT_ROOT, Date, sprintf('T%d-Method%s.png', LegLoop, Method)));
    close(f);
end

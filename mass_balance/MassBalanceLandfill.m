%% MassBalanceLandfill.m
%  Mass balance of the South Side Landfill (SSLF) plume on every downwind
%  transect (Section 2.4 of the paper).  Original script FindTheLandfill.m by
%  Zachary Barkley; the computations are unchanged.  This version adds the
%  settings block, the loop over all flight days, the header comments and the
%  REUSE_BOUNDS option.
%
%  For each flight day the script repeats the first steps of
%  MassBalanceMethodA.m (terrain-corrected altitude, mean wind of all
%  transects, boundary layer height from the profiles, the part of each
%  downwind transect that lies downwind of Marion County) and then, for each
%  downwind transect, uses the start and end points of the landfill plume in
%  the CH4 series of that transect segment:
%     Background   = mean of the CH4 values at the two end points of the plume
%     Enhancement  = mean(CH4 inside the plume) - Background            (ppm)
%     Peak         = max(CH4 inside the plume) - Background             (ppm)
%     Distance     = great-circle distance between the end points        (km)
%     Flux         = Density * MeanWSPD * (Enhancement/1e6) * Distance*1000
%                    * DownwindABL                                      (mol/s)
%  i.e. Equation 1 without the cos(theta) factor, which Collect_Landfill.m
%  applies afterwards.  Flag = 1 marks a quantified plume, Flag = -1 a transect
%  on which no landfill plume could be identified.
%
%  The plume end points were chosen by eye from the CH4 series of each transect
%  (interactive mode: the script plots the segment, stops with "keyboard", and
%  the user types A = <start index>, B = <end index>, optionally Flag = -1, then
%  dbcont).  The chosen indices are saved with the results (Start, End), so the
%  published values can be recomputed without interaction: with REUSE_BOUNDS =
%  true the script reads Start, End and Flag from the existing Loop files (or
%  from the released Landfill_transects.csv) and skips the prompt.
%
%  Output: Landfill/<yyyymmdd>_Loop<k>.mat with Flag, Distance, Background,
%  Enhancement, Peak, Flux, CH4LF (the CH4 series of the plume), Start, End.

clear
close all

%% ------------------------------------------------------------ settings
MB_ROOT  = '/storage/group/zrb5027/default/INFLUX/zrb5027/INFLUXMassBalance';   % FlightList.mat, LegAndProfileList/
OBS_DIR  = '/storage/work/zrb5027/WRFGUI_V3/Observations/INFLUX';               % aircraft files INFLUX_<yyyymmdd>
CODE_DIR = '/storage/work/zrb5027/WRFGUI_V3/Code';                              % InterX.m, ElevationData.nc
LF_IN    = fullfile(MB_ROOT, 'Landfill');   % existing results (read only): the stored plume end points when REUSE_BOUNDS is true
LF_OUT   = fullfile(pwd, 'Landfill');       % where the new Loop files are written
BOUNDS_CSV = '';                % alternative source of the end points: Landfill_transects.csv from Collect_Landfill.m
REUSE_BOUNDS = true;            % true: take Start/End/Flag from LF_IN (or BOUNDS_CSV), no interaction
START_AT = 1;                   % first flight day to process
SAVE_PROFILE_FIGURES = false;   % profile plots of each day into <yyyymmdd>/Profile<k>.png (as the original)

addpath(MB_ROOT); addpath(OBS_DIR); addpath(CODE_DIR);
assert(~isempty(which('InterX')) && ~isempty(which('ElevationData.nc')), 'InterX.m or ElevationData.nc not on the path (check CODE_DIR)');
if ~isfolder(LF_OUT), mkdir(LF_OUT); end
bounds = [];
if REUSE_BOUNDS && ~isempty(BOUNDS_CSV) && isfile(BOUNDS_CSV), bounds = readtable(BOUNDS_CSV, 'TextType', 'string', 'DatetimeType', 'text'); end
FL = load(fullfile(MB_ROOT, 'FlightList.mat')); FlightList = FL.FlightList;
for FlightDateLoop=START_AT:length(FlightList)
    clearvars -except FlightDateLoop FlightList MB_ROOT LF_IN LF_OUT REUSE_BOUNDS bounds SAVE_PROFILE_FIGURES
    Date=FlightList(FlightDateLoop).name;
    Date=Date(1:8);
    if SAVE_PROFILE_FIGURES && ~isfolder(Date), mkdir(Date); end
    fprintf('%2d %s\n', FlightDateLoop, Date);

SelfAssignedBackground=0;
Percentile=0.1;
load(fullfile(MB_ROOT,'LegAndProfileList',Date))
File=['INFLUX_',Date];
%BottomRight BottomLeft TopLeft TopRight
CountyCorners=[39.63760423971629, -85.95148383144084; 39.63337385377353, -86.32570560518701; 39.92446141560437, -86.32734193809195; 39.9271330073512, -85.95582253907553];

if exist('UpwindLeg','var')==1
    UpwindLeg=UpwindLeg(1):UpwindLeg(2);
    AllLeg=[UpwindLeg];
else
    AllLeg=[];
end
for LegLoop=1:size(LegList,1)
    AllLeg=[AllLeg,LegList(LegLoop,1):LegList(LegLoop,2);];
end

%% ------------------------------------------------------------ aircraft data and height above ground
x=load(File);
x(x<-998)=NaN;
Time=x(:,1);
Lat=x(:,2);
Lon=x(:,3);
Alt=x(:,4);

Elev=ncread('ElevationData.nc','data');
Elev(Elev<0)=0;
WorldLat=ncread('ElevationData.nc','lat');
WorldLon=ncread('ElevationData.nc','lon');
X=find(WorldLon>180);
WorldLon(X)=WorldLon(X)-360;
for i=1:length(Alt)
    NewLat=abs(WorldLat-Lat(i));
    NewLon=abs(WorldLon-Lon(i));
    Row=find(NewLon<0.2);
    Col=find(NewLat<0.2);
    if isempty(Row)==1
        Elevation(i,1)=NaN;
    else
        Elevation(i,1)=Elev(Row(1),Col(1));
    end
end
AltAGL=Alt-Elevation;

WSPD=x(:,11);
WDIR=x(:,12);
U=-1*WSPD.*sind(WDIR);
V=-1*WSPD.*cosd(WDIR);

CO2=x(:,5);
CH4=movmean(x(:,6),15,'omitnan');          % 15 s running mean of CH4 (ppm)
if exist('CH4Flagged','var')==1
    CH4=CH4Flagged;
end
H20=x(:,7);
CO=x(:,8);
T=x(:,10);
P=x(:,9);
MeanU=nanmean(U(AllLeg));                  % mean wind of all transects of the day
MeanV=nanmean(V(AllLeg));
MeanWSPD=(MeanU.^2+MeanV.^2).^.5;
MeanWDIR=atan2d(MeanU,MeanV)+180;

%% ------------------------------------------------------------ boundary layer height from the profiles (as MassBalanceMethodA.m)
for ProfileLoop=1:size(ProfileList,1)
    Profile=ProfileList(ProfileLoop,1):ProfileList(ProfileLoop,2);
    ProfileTime(ProfileLoop,1)=nanmean(Time(Profile));
    Hour=floor(ProfileTime(ProfileLoop)/3600);
    Minute=num2str(round((ProfileTime(ProfileLoop)-Hour*3600)/60));
    Hour=num2str(Hour);
    if length(Minute)<2
        Minute=['0',Minute];
    end

    Heights=-25:50:5025;
    AltAGLX=AltAGL(Profile);
    Theta=T(Profile).*(1000./P(Profile)).^(2/7);
    CH4X=CH4(Profile);
    PressureX=P(Profile);
    for HeightLoop=1:length(Heights)-1
        GoodSpot=find(AltAGLX>=Heights(HeightLoop) & AltAGLX<Heights(HeightLoop+1));
        CH4Averaged(HeightLoop)=nanmean(CH4X(GoodSpot));
        ThetaAveraged(HeightLoop)=nanmean(Theta(GoodSpot));
        PressureAveraged(HeightLoop)=nanmean(PressureX(GoodSpot));
    end
    Heights=0:50:5000;

    for Slope=1:length(Heights)-2
        ThetaDerivative(Slope)=(ThetaAveraged(Slope)-ThetaAveraged(Slope+2))/100;
        CH4Derivative(Slope)=(CH4Averaged(Slope)-CH4Averaged(Slope+2))/100;
    end
    ABLTop=find(CH4Derivative>2e-5 & ThetaDerivative<-5e-3);
    if exist('ManualABLHeights','var')==1
        ABLTopSpot=find(Heights==ManualABLHeights(ProfileLoop));
    else
        ABLTopSpot=ABLTop(1)+1;
    end
    ABLTop=Heights(ABLTopSpot);
    ABLHeight(ProfileLoop,1)=ABLTop;
    ABLTheta(ProfileLoop,1)=nanmean(ThetaAveraged(1:ABLTopSpot));
    ABLPressure(ProfileLoop,1)=nanmean(PressureAveraged(1:ABLTopSpot));
    ABLT(ProfileLoop,1)=ABLTheta(ProfileLoop)./((1000./ABLPressure(ProfileLoop)).^(2/7));

    if SAVE_PROFILE_FIGURES
        f=figure('Color',[1 1 1],'visible','off');
        t = tiledlayout(1,1);
        ax1 = axes(t);
        plot(ax1,CH4Averaged,Heights,'-k')
        ax1.XColor = 'k'; ax1.YColor = 'k';
        xlabel('CH_4 (ppm)')
        ax2 = axes(t);
        plot(ax2,ThetaAveraged,Heights,'-r')
        ax2.XAxisLocation = 'top'; ax2.YAxisLocation = 'left'; ax2.Color = 'none';
        ax1.Box = 'off'; ax2.Box = 'off';
        set(ax2,'xcolor','red')
        xlabel('Pot Temp (K)')
        title([Date(5:6),'/',Date(7:8),'/',Date(1:4),'  ','Profile ',num2str(ProfileLoop),'  ',Hour,':',Minute,'Z'])
        hold on
        ABLTopLine=[]; ABLTopLine(1:length(ThetaAveraged))=ABLTop;
        plot(ax2,ThetaAveraged,ABLTopLine,'-k','LineWidth',3)
        annotation('textbox', [0.4, 0.15, 0.5, 0.05], 'String', ['ABL Height = ', num2str(ABLTop)]);
        annotation('rectangle', [0.4, 0.15, 0.5, 0.05],'color',[1 1 1]);
        saveas(f,[Date,'/Profile',num2str(ProfileLoop),'.png'])
        close(f)
    end
end

for LegLoop=1:size(LegList,1)
    Leg=LegList(LegLoop,1):LegList(LegLoop,2);
    DownwindTime(LegLoop,1)=nanmean(Time(Leg));
end

if length(ProfileTime)==1                  % one profile only: assume a growing boundary layer
    if ProfileTime<min(DownwindTime)
        ProfileTime(2,1)=max(DownwindTime);
        ABLHeight(2,1)=ABLHeight*1.2;
        ProfileTime(2,1)=max(DownwindTime);
    elseif ProfileTime>max(DownwindTime)
        ABLHeight(2,1)=ABLHeight;
        ABLHeight(1,1)=ABLHeight(2,1)*.833;
        ProfileTime(2,1)=ProfileTime;
        ProfileTime(1,1)=min(DownwindTime);
    end
    ABLT(2,1)=ABLT;
    ABLPressure(2,1)=ABLPressure;
end
DownwindTime(DownwindTime>max(ProfileTime))=max(ProfileTime);
DownwindTime(DownwindTime<min(ProfileTime))=min(ProfileTime);

if length(ProfileTime)>1
    ABLT(isnan(ABLT)==1)=nanmean(ABLT);
    ABLPressure(isnan(ABLPressure)==1)=nanmean(ABLPressure);
    DownwindABL=interp1(ProfileTime,ABLHeight,DownwindTime);
    DownwindT=interp1(ProfileTime,ABLT,DownwindTime);
    DownwindPressure=interp1(ProfileTime,ABLPressure,DownwindTime);
end
Density=DownwindPressure*100./8.314./DownwindT;    % mol/m3

%% ------------------------------------------------------------ the part of each transect downwind of Marion County (as MassBalanceMethodA.m)
if (MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0)
    PointB1=CountyCorners(3,:);
    PointB2=CountyCorners(1,:);
elseif (MeanU>0 & MeanV<0) | (MeanU<0 & MeanV>0)
    PointB1=CountyCorners(2,:);
    PointB2=CountyCorners(4,:);
end
PointC=[PointB1(1)+MeanV, PointB1(2)+MeanU*(111/86)];
PointA=[PointB1(1)-MeanV, PointB1(2)-MeanU*(111/86)];
Line1a=[PointA(2) PointB1(2) PointC(2);PointA(1) PointB1(1) PointC(1)];
PointC=[PointB2(1)+MeanV, PointB2(2)+MeanU*(111/86)];
PointA=[PointB2(1)-MeanV, PointB2(2)-MeanU*(111/86)];
Line1b=[PointA(2) PointB2(2) PointC(2);PointA(1) PointB2(1) PointC(1)];
for LegLoop=1:size(LegList,1)
    Leg=LegList(LegLoop,1):LegList(LegLoop,2);
    LegHeight(LegLoop)=nanmean(AltAGL(Leg));
    Line2=[Lon(Leg)';Lat(Leg)'];
    Intersection1=flipud(InterX(Line1a,Line2))';
    Intersection2=flipud(InterX(Line1b,Line2))';
    if ((MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0))
        if isempty(Intersection1)==1
            Intersection1=[max(Lat(Leg)) min(Lon(Leg))];
        end
        if isempty(Intersection2)==1
            Intersection2=[min(Lat(Leg)) max(Lon(Leg))];
        end
        if abs(Intersection1(1)-Intersection2(1)) < 0.05
            Upwind1=find(Lon(Leg)<Intersection1(2));
        elseif abs(Intersection1(2)-Intersection2(2)) < 0.05
            Upwind1=find(Lat(Leg)>Intersection1(1));
        else
            Upwind1=find(Lat(Leg)>Intersection1(1) | Lon(Leg)<Intersection1(2));
        end
        if abs(Intersection1(1)-Intersection2(1)) < 0.05
            Upwind2=find(Lon(Leg)>Intersection2(2));
        elseif abs(Intersection2(2)-Intersection2(2)) < 0.05    % always true in the original script; kept so that the stored Start/End indices apply
            Upwind2=find(Lat(Leg)<Intersection2(1));
        else
            Upwind2=find(Lat(Leg)<Intersection2(1) | Lon(Leg)>Intersection2(2));
        end
    elseif((MeanU>0 & MeanV<0) | (MeanU<0 & MeanV>0))
        if isempty(Intersection1)==1
            Intersection1=[min(Lat(Leg)) min(Lon(Leg))];
        end
        if isempty(Intersection2)==1
            Intersection2=[max(Lat(Leg)) max(Lon(Leg))];
        end
        if abs(Intersection1(1)-Intersection2(1)) < 0.05
            Upwind1=find(Lon(Leg)<Intersection1(2));
        elseif abs(Intersection1(2)-Intersection2(2)) < 0.05
            Upwind1=find(Lat(Leg)<Intersection1(1));
        else
            Upwind1=find(Lat(Leg)<Intersection1(1) | Lon(Leg)<Intersection1(2));
        end
        if abs(Intersection1(1)-Intersection2(1)) < 0.05
            Upwind2=find(Lon(Leg)>Intersection2(2));
        elseif abs(Intersection2(2)-Intersection2(2)) < 0.05    % always true in the original script; kept so that the stored Start/End indices apply
            Upwind2=find(Lat(Leg)>Intersection2(1));
        else
            Upwind2=find(Lat(Leg)>Intersection2(1) | Lon(Leg)>Intersection2(2));
        end
    end
    Intersection=[Intersection1;Intersection2];
    Upwind1=Leg(Upwind1);
    Upwind2=Leg(Upwind2);
    Temp=ismember(Leg,[Upwind1 Upwind2]);
    Downwind=Leg(Temp==0);                 % indices of the transect points downwind of the county

    %% -------------------------------------------------------- landfill plume: start (A) and end (B) index inside Downwind
    outFile=fullfile(LF_OUT,[Date,'_Loop',num2str(LegLoop),'.mat']);
    refFile=fullfile(LF_IN,[Date,'_Loop',num2str(LegLoop),'.mat']);
    A=[]; B=[]; Flag=1;
    if REUSE_BOUNDS
        if isfile(refFile)
            old=load(refFile);
            if isfield(old,'Start'), A=old.Start; B=old.End; Flag=old.Flag; end
        elseif ~isempty(bounds)
            r=find(string(bounds.Date)==sprintf('%s-%s-%s',Date(1:4),Date(5:6),Date(7:8)) & bounds.Transect==LegLoop, 1);
            if ~isempty(r), A=bounds.Start(r); B=bounds.End(r); Flag=bounds.Flag(r); end
        end
        if isempty(A)
            fprintf('   transect %d: no stored plume bounds, skipped\n', LegLoop);
            continue
        end
    else
        f=figure('visible','on');
        plot(CH4(Downwind)); title(sprintf('%s transect %d: type A = start, B = end of the landfill plume, then dbcont', Date, LegLoop))
        Flag=1;
        keyboard
        close all
        %%Identify start (A) and end point (B), then continue. Enter Flag=-1 if
        %%you don't know what's going on.
    end
    if Flag<0
        save(outFile,'Flag'); continue
    end
    Start=A; End=B;
    Range=Downwind(A:B);
    CH4LF=CH4(Range);
    LatLF=Lat(Range);
    LonLF=Lon(Range);
    Distance=distance(LatLF(1),LonLF(1),LatLF(end),LonLF(end))*111;    % km (degrees of arc x 111 km)
    Background=mean([CH4LF(1),CH4LF(end)]);
    Enhancement=mean(CH4LF)-Background;
    Peak=max(CH4LF(:))-Background;
    Flux=Density(LegLoop).*MeanWSPD.*(Enhancement/1e6).*(Distance*1000).*DownwindABL(LegLoop);
    save(outFile,'Flag','Distance','Background','Enhancement','Peak','Flux','CH4LF','Start','End')
    fprintf('   transect %d: plume %d to %d, %.1f km, %.1f ppb, %.1f mol/s\n', LegLoop, Start, End, Distance, Enhancement*1000, Flux);
end

end

clear
close all
addpath('/storage/work/zrb5027/WRFGUI_V3/Observations/INFLUX');
addpath('/storage/work/zrb5027/WRFGUI_V3/Code')
load FlightList
START_AT=1;   % set to k to resume at flight k (e.g. 4 = 20110630)
for FlightDateLoop=START_AT:length(FlightList)
    clearvars -except FlightDateLoop
    load FlightList
    Date=FlightList(FlightDateLoop).name;
    Date=Date(1:8);
   % Date='20141117';
    mkdir(Date);

%Date='20141110';
SelfAssignedBackground=0;
BackgroundPercentile=0.2;
load(['LegAndProfileList/',Date])
File=['INFLUX_',Date];
%BottomRight BottomLeft TopLeft TopRight
CountyCorners=[39.63760423971629, -85.95148383144084; 39.63337385377353, -86.32570560518701; 39.92446141560437, -86.32734193809195; 39.9271330073512, -85.95582253907553];

if exist('UpwindLeg','var')==1
UpwindLeg=UpwindLeg(1):UpwindLeg(2);
AllLeg=[UpwindLeg];
for LegLoop=1:size(LegList,1)
    AllLeg=[AllLeg,LegList(LegLoop,1):LegList(LegLoop,2);];
end




%111/86   U multiplier for line

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
                    %keyboard
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
CH4=movmean(x(:,6),30,'omitnan');
if exist('CH4Flagged','var')==1
    CH4=CH4Flagged;
end
H20=x(:,7);
CO=x(:,8);
T=x(:,10);
P=x(:,9);
MeanU=nanmean(U(AllLeg));
MeanV=nanmean(V(AllLeg));
MeanWSPD=(MeanU.^2+MeanV.^2).^.5;
MeanWDIR=atan2d(MeanU,MeanV)+180;



for ProfileLoop=1:size(ProfileList,1)
    Profile=ProfileList(ProfileLoop,1):ProfileList(ProfileLoop,2);
    ProfileTime(ProfileLoop,1)=nanmean(Time(Profile))
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

    %CH4X=interp1(AltAGLX,CH4X,Heights);
    %ThetaX=interp1(AltAGLX,Theta,Heights);
    

    f=figure('Color',[1 1 1],'visible','off')
    
t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1,CH4Averaged,Heights,'-k')
ax1.XColor = 'k';
ax1.YColor = 'k';
xlabel('CH_4 (ppm)')

ax2 = axes(t);
plot(ax2,ThetaAveraged,Heights,'-r')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'left';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
set(ax2,'xcolor','red')
xlabel('Pot Temp (K)')
title([Date(5:6),'/',Date(7:8),'/',Date(1:4),'  ','Profile ',num2str(ProfileLoop),'  ',Hour,':',Minute,'Z'])
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
%if ProfileLoop==2
%    ABLTopSpot=30
%end
ABLTop=Heights(ABLTopSpot);
ABLTopLine=[];
ABLTopLine(1:length(ThetaAveraged))=ABLTop;
hold on
plot(ax2,ThetaAveraged,ABLTopLine,'-k','LineWidth',3)
ABLHeight(ProfileLoop,1)=ABLTop;
annotation('textbox', [0.4, 0.15, 0.5, 0.05], 'String', ['ABL Height = ', num2str(ABLTopLine(1))]);
annotation('rectangle', [0.4, 0.15, 0.5, 0.05],'color',[1 1 1]);
ABLTheta(ProfileLoop,1)=nanmean(ThetaAveraged(1:ABLTopSpot));
ABLPressure(ProfileLoop,1)=nanmean(PressureAveraged(1:ABLTopSpot));
ABLT(ProfileLoop,1)=ABLTheta(ProfileLoop)./((1000./ABLPressure(ProfileLoop)).^(2/7));
if ProfileLoop==1
    CH4ProfileOG=CH4Averaged;
    FirstPoint=min(find(isnan(CH4ProfileOG)==0));
    if FirstPoint~=1
        CH4ProfileOG(1:FirstPoint-1)=mean(CH4ProfileOG(FirstPoint:FirstPoint+3));
    end
    LastPoint=max(find(isnan(CH4ProfileOG)==0));
    if LastPoint~=length(Heights)
        CH4ProfileOG(LastPoint+1:end)=CH4ProfileOG(LastPoint);
    end
Heights2=0:10:5000;
CH4ProfileOG=interp1(Heights,CH4ProfileOG,Heights2);
end
saveas(f,[Date,'/Profile',num2str(ProfileLoop),'.png'])
close(f)   % free memory: hidden figures were never closed
end



%%Start of area of interest
for LegLoop=1:size(LegList,1)
    Leg=LegList(LegLoop,1):LegList(LegLoop,2);
    DownwindTime(LegLoop,1)=nanmean(Time(Leg))
end

UpwindHeight=nanmean(AltAGL(UpwindLeg));
UpwindTime=nanmean(Time(UpwindLeg));
if length(ProfileTime)==1
    if ProfileTime<min(DownwindTime);
        ProfileTime(2,1)=max(DownwindTime);
        ABLHeight(2,1)=ABLHeight*1.2;
        ProfileTime(2,1)=max(DownwindTime)
    elseif ProfileTime>max(DownwindTime);
        error('Lacking early morning profile')
        ABLHeight(2,1)=ABLHeight;
        ABLHeight(1,1)=ABLHeight(2,1)*.833;
        ProfileTime(2,1)=ProfileTime
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
Density=DownwindPressure*100./8.314./DownwindT;

for DownwindLoop=1:length(DownwindABL)
    %Calculate Adjustment factor
    Temp=abs(Heights2-DownwindABL(DownwindLoop));
    GoodSpot=find(Temp==min(Temp));
    %Integrate CH4 to GoodSpot
    MeanABLCH4=mean(CH4ProfileOG(1:GoodSpot));
    %Find height of upwind and estimate roughly what part of the CH4
    %profile it was measuring
    Temp=abs(Heights2-UpwindHeight);
    GoodSpot=find(Temp==min(Temp));
    ReasonableUpwindCH4=nanmean(CH4ProfileOG(GoodSpot-7:GoodSpot+7));
%    CorrectionFactor(DownwindLoop)=MeanABLCH4/ReasonableUpwindCH4;
    % pause(2)
end
Density=DownwindPressure*100./8.314./DownwindT;

%%END OF ABL Cals


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
    %--- cos(theta): angle between the mean wind and the normal to this leg (Eq. 1)
    LatLeg=Lat(Leg); LonLeg=Lon(Leg); okLeg=~isnan(LatLeg) & ~isnan(LonLeg);
    xKm=(LonLeg(okLeg)-mean(LonLeg(okLeg)))*111.32*cosd(mean(LatLeg(okLeg)));
    yKm=(LatLeg(okLeg)-mean(LatLeg(okLeg)))*110.57;
    [EigV,EigD]=eig(cov([xKm yKm])); [~,mAx]=max(diag(EigD)); eLeg=EigV(:,mAx)';
    uFlow=[MeanU MeanV]/MeanWSPD;                      % unit vector of the flow (toward)
    CosTheta(LegLoop)=sqrt(max(0,1-(uFlow*eLeg')^2));  % = cos of angle to the leg normal

    %Create Line
    Line2=[Lon(Leg)';Lat(Leg)'];

        Intersection1=flipud(InterX(Line1a,Line2))';
        Intersection2=flipud(InterX(Line1b,Line2))';
    if ((MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0))
        if isempty(Intersection1)==1
            Intersection1=[max(Lat(Leg)) min(Lon(Leg))]
        end
        if isempty(Intersection2)==1
            Intersection2=[min(Lat(Leg)) max(Lon(Leg))]
        end
        if abs(Intersection1(1)-Intersection2(1)) < 0.01
                   Upwind1=find(Lon(Leg)<Intersection1(2));
        elseif abs(Intersection1(2)-Intersection2(2)) < 0.01
                   Upwind1=find(Lat(Leg)>Intersection1(1));
        else
                   Upwind1=find(Lat(Leg)>Intersection1(1) | Lon(Leg)<Intersection1(2));
        end
    

        
        if abs(Intersection1(1)-Intersection2(1)) < 0.01
                   Upwind2=find(Lon(Leg)>Intersection2(2));
        elseif abs(Intersection2(2)-Intersection2(2)) < 0.01
                   Upwind2=find(Lat(Leg)<Intersection2(1));
        else
                   Upwind2=find(Lat(Leg)<Intersection2(1) | Lon(Leg)>Intersection2(2));
        end
        
    elseif((MeanU>0 & MeanV<0) | (MeanU<0 & MeanV>0))
        if isempty(Intersection1)==1
            Intersection1=[min(Lat(Leg)) min(Lon(Leg))]
        end
        if isempty(Intersection2)==1
            Intersection2=[max(Lat(Leg)) max(Lon(Leg))]
        end
        if abs(Intersection1(1)-Intersection2(1)) < 0.01
                   Upwind1=find(Lon(Leg)<Intersection1(2));
        elseif abs(Intersection1(2)-Intersection2(2)) < 0.01
                   Upwind1=find(Lat(Leg)<Intersection1(1));
        else
                   Upwind1=find(Lat(Leg)<Intersection1(1) | Lon(Leg)<Intersection1(2));
        end
        
 
          if abs(Intersection1(1)-Intersection2(1)) < 0.01
                   Upwind2=find(Lon(Leg)>Intersection2(2));
        elseif abs(Intersection2(2)-Intersection2(2)) < 0.01
                   Upwind2=find(Lat(Leg)>Intersection2(1));
        else
                   Upwind2=find(Lat(Leg)>Intersection2(1) | Lon(Leg)>Intersection2(2));
        end
  
    end
    Intersection=[Intersection1;Intersection2];
        Upwind1=Leg(Upwind1);
        Upwind2=Leg(Upwind2);
        Temp=ismember(Leg,[Upwind1 Upwind2]);
        Downwind=Leg(Temp==0);
        UpwindLine=[Lon(UpwindLeg), Lat(UpwindLeg)];
        BackgroundLine=[];
        for UpwindLoop=1:length(Leg);
            PointB=[Lat(Leg(UpwindLoop)), Lon(Leg(UpwindLoop))];
            PointC=[PointB(1)+MeanV, PointB(2)+MeanU*(111/86)];
            PointA=[PointB(1)-MeanV, PointB(2)-MeanU*(111/86)];
            Line1A=[PointA(2) PointB(2) PointC(2);PointA(1) PointB(1) PointC(1)];
            UpwindIntersection=flipud(InterX(Line1A,UpwindLine'))';
            if isempty(UpwindIntersection)==1
                BackgroundLine(UpwindLoop)=NaN;
            else
            TempDistance=distance(UpwindIntersection(1),UpwindIntersection(2),UpwindLine(:,2),UpwindLine(:,1));
            GoodSpot=find(TempDistance==min(TempDistance));
            UpwindAligned=UpwindLeg(GoodSpot);
            BackgroundLine(UpwindLoop)=CH4(UpwindAligned);%*CorrectionFactor(LegLoop);
            end
            
        end
        %BackgroundLine=movmean(BackgroundLine,60,'omitnan');
        Temp=movmean(BackgroundLine,120,'omitnan');   
        BadSpot=find(isnan(BackgroundLine)==1);
        BackgroundLine(BadSpot)=Temp(BadSpot);
        
        Temp=movmean(BackgroundLine,240,'omitnan');   
        BadSpot=find(isnan(BackgroundLine)==1);
        BackgroundLine(BadSpot)=Temp(BadSpot);        
        Temp=movmean(BackgroundLine,360,'omitnan');   
        BadSpot=find(isnan(BackgroundLine)==1);
        BackgroundLine(BadSpot)=Temp(BadSpot);                 
                
DownwindLine=CH4(Leg);
ShiftCounter=0;
 for Shift=-50:50
     ShiftCounter=ShiftCounter+1;
TempBackground=NaN(length(BackgroundLine)+100,1);
TempDownwind=NaN(length(BackgroundLine)+100,1);
TempDownwind(51:end-50)=DownwindLine;
TempBackground(51+Shift:end-50+Shift)=BackgroundLine;
TempCorr=corrcoef(TempBackground,TempDownwind,'Rows','complete');
R(ShiftCounter)=TempCorr(2);
end
 Shift=find(R==max(R))-51;
 BackgroundLineDoppel=NaN(length(BackgroundLine)+100,1);
 BackgroundLineDoppel(51+Shift:end-50+Shift)=BackgroundLine;

 BackgroundLine=BackgroundLineDoppel(51:end-50);
 BackgroundLine=BackgroundLine';

 
        BackgroundGoodSpot=find(isnan(BackgroundLine)==0);

        Temp=movmean(BackgroundLine,120,'omitnan');   
        BadSpot=find(isnan(BackgroundLine)==1);
        BackgroundLine(BadSpot)=Temp(BadSpot);
        
        Temp=movmean(BackgroundLine,240,'omitnan');   
        BadSpot=find(isnan(BackgroundLine)==1);
        BackgroundLine(BadSpot)=Temp(BadSpot);        
        Temp=movmean(BackgroundLine,360,'omitnan');   
        BadSpot=find(isnan(BackgroundLine)==1);
        BackgroundLine(BadSpot)=Temp(BadSpot);         

        
        
        BackgroundLine(isnan(BackgroundLine)==1)=nanmedian(BackgroundLine);
        %BadSpot=find((DownwindLine-BackgroundLine')<0 & isnan(DownwindLine)==0 & isnan(BackgroundLine')==0);
%                 BadSpot=find(isnan(DownwindLine)==0 & isnan(BackgroundLine')==0);
%         TempSpot=round(length(BadSpot)*0.1); %Percentile Adjustment here Zach!
%         TempSort1=sort(BackgroundLine(BadSpot));
%         TempSort2=sort(DownwindLine(BadSpot));
%         CorrectionFactor2(LegLoop)=TempSort2(TempSpot)/TempSort1(TempSpot);
%        % CorrectionFactor2(LegLoop)=1;
% 
%         %CorrectionFactor2(LegLoop)=nanmedian(DownwindLine(BadSpot)/nanmedian(BackgroundLine(BadSpot));  
%         BackgroundLine=BackgroundLine*CorrectionFactor2(LegLoop);

        
f=figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1],'visible','off')

subplot(1,2,1)
worldmap([39.2 40.6],[-87 -85.2])
setm(gca,'FontSize',14)

hold on
scatterm(Lat(Upwind1),Lon(Upwind1),200,[0 0 0],'filled')
scatterm(Lat(Downwind),Lon(Downwind),200,[1 0 0],'filled')
scatterm(Lat(Upwind2),Lon(Upwind2),200,[0 0 0],'filled')
scatterm(Lat(UpwindLeg),Lon(UpwindLeg),200,[0 0 0],'filled')
scatterm(Lat(UpwindLeg),Lon(UpwindLeg),50,CH4(UpwindLeg),'filled')

scatterm(Lat(Leg),Lon(Leg),50,DownwindLine,'filled')
c=colorbar
c.FontSize=14;
PlotINFLUX
%plotm([PointB1(1),Intersection(1,1)],[PointB1(2),Intersection(1,2)],'LineWidth',2,'Color',[0 0 0])
%plotm([PointB2(1),Intersection(2,1)],[PointB2(2),Intersection(2,2)],'LineWidth',2,'Color',[0 0 0])
plotm([CountyCorners(1,1),CountyCorners(2,1)],[CountyCorners(1,2),CountyCorners(2,2)],'LineWidth',4,'Color',[0 0 0])
plotm([CountyCorners(2,1),CountyCorners(3,1)],[CountyCorners(2,2),CountyCorners(3,2)],'LineWidth',4,'Color',[0 0 0])
plotm([CountyCorners(3,1),CountyCorners(4,1)],[CountyCorners(3,2),CountyCorners(4,2)],'LineWidth',4,'Color',[0 0 0])
plotm([CountyCorners(4,1),CountyCorners(1,1)],[CountyCorners(4,2),CountyCorners(1,2)],'LineWidth',4,'Color',[0 0 0])


subplot(1,2,2)

%Make Lineplot
LatDistance=distance(nanmin(Lat(Leg)),nanmean(Lon(Leg)),nanmax(Lat(Leg)),nanmean(Lon(Leg)));
LonDistance=distance(nanmean(Lat(Leg)),nanmin(Lon(Leg)),nanmean(Lat(Leg)),nanmax(Lon(Leg)));
if LonDistance>LatDistance
    CH4Downwind=CH4(Downwind);
    BackgroundLineDownwind=BackgroundLine(ismember(Leg,Downwind)==1)';
    x=Lon(Downwind);
    y1=CH4Downwind;
    Temp=y1-BackgroundLineDownwind;
    BadSpot=find(Temp<0)
    y1(BadSpot)=BackgroundLineDownwind(BadSpot);
 %   GS=find(CH4Downwind-BackgroundLine)>0;
 GoodSpot=find(isnan(y1)==0);
 Temp=[x(GoodSpot)' flip(x(GoodSpot))'];
 Temp2=[BackgroundLineDownwind(GoodSpot)',flip(y1(GoodSpot))'];
 patch(Temp,Temp2, 'g')
     y1=CH4Downwind;
    Temp=y1-BackgroundLineDownwind;
    BadSpot=find(Temp>0)
    y1(BadSpot)=BackgroundLineDownwind(BadSpot);
 %   GS=find(CH4Downwind-BackgroundLine)>0;
 GoodSpot=find(isnan(y1)==0);
 Temp=[x(GoodSpot)' flip(x(GoodSpot))'];
 Temp2=[BackgroundLineDownwind(GoodSpot)',flip(y1(GoodSpot))'];
 patch(Temp,Temp2, 'r')
 xlabel('Longitude')
 hold on
     plot(Lon(Leg),DownwindLine,'Color',[0 0 0],'LineWidth',2)    
    plot(Lon(Leg(BackgroundGoodSpot)),BackgroundLine(BackgroundGoodSpot),'Color','blue','LineWidth',2)
 
 
 LonIntp=min(Intersection(:,2)):0.001:max(Intersection(:,2))
 TempMatrix=[x(GoodSpot),CH4Downwind(GoodSpot), BackgroundLineDownwind(GoodSpot)];
 [junk,UniqueRows]=unique(TempMatrix(:,1));
 TempMatrix=TempMatrix(UniqueRows,:);
 TempMatrix=sortrows(TempMatrix,1);
DownwindReimagined=interp1(TempMatrix(:,1),TempMatrix(:,2),LonIntp);
UpwindReimagined=interp1(TempMatrix(:,1),TempMatrix(:,3),LonIntp);
DownwindCH4Intp=nanmean(DownwindReimagined);
UpwindCH4Intp=nanmean(UpwindReimagined);
elseif LatDistance>LonDistance
    
    CH4Downwind=CH4(Downwind);
    BackgroundLineDownwind=BackgroundLine(ismember(Leg,Downwind)==1)';
    x=Lat(Downwind);
    y1=CH4Downwind;
    Temp=y1-BackgroundLineDownwind;
    BadSpot=find(Temp<0)
    y1(BadSpot)=BackgroundLineDownwind(BadSpot);
 %   GS=find(CH4Downwind-BackgroundLine)>0;
 GoodSpot=find(isnan(y1)==0);
 Temp=[x(GoodSpot)' flip(x(GoodSpot))'];
 Temp2=[BackgroundLineDownwind(GoodSpot)',flip(y1(GoodSpot))'];
 patch(Temp,Temp2, 'g')
     y1=CH4Downwind;
    Temp=y1-BackgroundLineDownwind;
    BadSpot=find(Temp>0)
    y1(BadSpot)=BackgroundLineDownwind(BadSpot);
 %   GS=find(CH4Downwind-BackgroundLine)>0;
 GoodSpot=find(isnan(y1)==0);
 Temp=[x(GoodSpot)' flip(x(GoodSpot))'];
 Temp2=[BackgroundLineDownwind(GoodSpot)',flip(y1(GoodSpot))'];
 patch(Temp,Temp2, 'r')
 xlabel('Latitude')
 hold on
     plot(Lat(Leg),DownwindLine,'Color',[0 0 0],'LineWidth',2)    
    plot(Lat(Leg),BackgroundLine,'Color','blue','LineWidth',2)
 
 
 LatIntp=min(Intersection(:,1)):0.001:max(Intersection(:,1))
 TempMatrix=[x(GoodSpot),CH4Downwind(GoodSpot), BackgroundLineDownwind(GoodSpot)];
 [junk,UniqueRows]=unique(TempMatrix(:,1));
 TempMatrix=TempMatrix(UniqueRows,:);
 TempMatrix=sortrows(TempMatrix,1);
DownwindReimagined=interp1(TempMatrix(:,1),TempMatrix(:,2),LatIntp);
UpwindReimagined=interp1(TempMatrix(:,1),TempMatrix(:,3),LatIntp);
DownwindCH4Intp=nanmean(DownwindReimagined);
UpwindCH4Intp=nanmean(UpwindReimagined);
end
ylabel('CH_4 (ppm)')
set(gca,'FontSize',14) 
Hour=floor(DownwindTime(LegLoop,1)/3600);
    Minute=num2str(round((DownwindTime(LegLoop,1)-Hour*3600)/60));
    Hour=num2str(Hour);
    if length(Minute)<2
        Minute=['0',Minute];
    end


sgtitle([Date(5:6),'/',Date(7:8),'/',Date(1:4),'  ',Hour,':',Minute,'Z','  (T',num2str(LegLoop),') - Method B'],'FontSize',18)

%UpwindMatrix(LegLoop,1:2)=[nanmean(CH4(Upwind1)) nanmean(CH4(Upwind2))];
UpwindMatrix(LegLoop,1)=UpwindCH4Intp;
DownwindMatrix(LegLoop,1)=DownwindCH4Intp;
Distance(LegLoop,1)=distance(Intersection(1,:),Intersection(2,:))*111; 

Enhancement(LegLoop)=DownwindCH4Intp-UpwindCH4Intp;
Flux(LegLoop)=Density(LegLoop).*MeanWSPD.*CosTheta(LegLoop).*(Enhancement(LegLoop)/1e6).*(Distance(LegLoop)*1000).*DownwindABL(LegLoop);
annotation('textbox', [0.68, 0.9, 0.15, 0.05], 'String', ['Flux = ', num2str(round(Flux(LegLoop))),' mol/s']);
annotation('rectangle', [0.68, 0.9, 0.15, 0.05],'color',[1 1 1]);
saveas(f,[Date,'/','T',num2str(LegLoop),'-MethodB','.png'])
close(f)   % free memory: hidden figures were never closed
end


    


clearvars -except DownwindMatrix UpwindMatrix FlightDateLoop CorrectionFactor CorrectionFactor2 Date UpwindMatrix DownwindMatrix DownwindABL LegHeight Distance Enhancement UpwindTime DownwindTime ABLHeight ABLTheta ABLT ABLPressure ProfileTime MeanWSPD MeanWDIR Density Flux CosTheta
save([Date,'/','Solution-MethodB'])
end
end










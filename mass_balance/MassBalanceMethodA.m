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
 %   Date='20120606';
    mkdir(Date);


%Date='20121004';
SelfAssignedBackground=0;
Percentile=0.1;
load(['LegAndProfileList/',Date])
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




%111/86   U multiplier for line

x=load(File);
% keyboard   (removed for batch run)
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

% keyboard   (removed for batch run)
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
saveas(f,[Date,'/Profile',num2str(ProfileLoop),'.png'])
close(f)   % free memory: hidden figures were never closed
end

% pause(2)


%%Start of area of interest
for LegLoop=1:size(LegList,1)
    Leg=LegList(LegLoop,1):LegList(LegLoop,2);
    DownwindTime(LegLoop,1)=nanmean(Time(Leg))
end

if length(ProfileTime)==1
    if ProfileTime<min(DownwindTime);
        ProfileTime(2,1)=max(DownwindTime);
        ABLHeight(2,1)=ABLHeight*1.2;
        ProfileTime(2,1)=max(DownwindTime)
    elseif ProfileTime>max(DownwindTime);
%        error('Lacking early morning profile')
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
        if abs(Intersection1(1)-Intersection2(1)) < 0.05
                   Upwind1=find(Lon(Leg)<Intersection1(2));
        elseif abs(Intersection1(2)-Intersection2(2)) < 0.05
                   Upwind1=find(Lat(Leg)>Intersection1(1));
        else
                   Upwind1=find(Lat(Leg)>Intersection1(1) | Lon(Leg)<Intersection1(2));
        end
    

        
        if abs(Intersection1(1)-Intersection2(1)) < 0.05
                   Upwind2=find(Lon(Leg)>Intersection2(2));
        elseif abs(Intersection2(2)-Intersection2(2)) < 0.05
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
        if abs(Intersection1(1)-Intersection2(1)) < 0.05
                   Upwind1=find(Lon(Leg)<Intersection1(2));
        elseif abs(Intersection1(2)-Intersection2(2)) < 0.05
                   Upwind1=find(Lat(Leg)<Intersection1(1));
        else
                   Upwind1=find(Lat(Leg)<Intersection1(1) | Lon(Leg)<Intersection1(2));
        end

          if abs(Intersection1(1)-Intersection2(1)) < 0.05
                   Upwind2=find(Lon(Leg)>Intersection2(2));
        elseif abs(Intersection2(2)-Intersection2(2)) < 0.05
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
        
f=figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1],'visible','off')

subplot(1,2,1)
worldmap([39.2 40.6],[-87 -85.2])
setm(gca,'FontSize',14)

hold on
scatterm(Lat(Upwind1),Lon(Upwind1),200,[0 0 0],'filled')
scatterm(Lat(Downwind),Lon(Downwind),200,[1 0 0],'filled')
scatterm(Lat(Upwind2),Lon(Upwind2),200,[0 0 0],'filled')


scatterm(Lat(Leg),Lon(Leg),50,CH4(Leg),'filled')
c=colorbar
c.FontSize=14;
PlotINFLUX
%plotm([PointB1(1),Intersection(1,1)],[PointB1(2),Intersection(1,2)],'LineWidth',2,'Color',[0 0 0])
%plotm([PointB2(1),Intersection(2,1)],[PointB2(2),Intersection(2,2)],'LineWidth',2,'Color',[0 0 0])
plotm([CountyCorners(1,1),CountyCorners(2,1)],[CountyCorners(1,2),CountyCorners(2,2)],'LineWidth',4,'Color',[0 0 0])
plotm([CountyCorners(2,1),CountyCorners(3,1)],[CountyCorners(2,2),CountyCorners(3,2)],'LineWidth',4,'Color',[0 0 0])
plotm([CountyCorners(3,1),CountyCorners(4,1)],[CountyCorners(3,2),CountyCorners(4,2)],'LineWidth',4,'Color',[0 0 0])
plotm([CountyCorners(4,1),CountyCorners(1,1)],[CountyCorners(4,2),CountyCorners(1,2)],'LineWidth',4,'Color',[0 0 0])


if isempty(Upwind1)==0
    Upwind1Flag=1;
else
    Upwind1Flag=-1;
end
if isempty(Upwind2)==0
    Upwind2Flag=1;
else
    Upwind2Flag=-1
end


Upwind1Background=nanmedian(CH4(Upwind1));
Upwind2Background=nanmedian(CH4(Upwind2));
mTemp=round(length(Leg)/10);
if Upwind1Flag==-1;
    CheckSpot1=Leg(1)-Upwind2(1);
    CheckSpot2=Leg(end)-Upwind2(end);
    if CheckSpot1==0
        Upwind1Background=nanmedian(CH4(Leg(end-mTemp:end)));
    elseif CheckSpot2==0;
        Upwind1Background=nanmedian(CH4(Leg(1:mTemp)));
    else
        error('Check upwind edge methodology')
    end       
elseif Upwind2Flag==-1
    CheckSpot1=Leg(1)-Upwind1(1);
    CheckSpot2=Leg(end)-Upwind1(end);
    if CheckSpot1==0
        Upwind2Background=nanmedian(CH4(Leg(end-mTemp:end)));
    elseif CheckSpot2==0;
        Upwind2Background=nanmedian(CH4(Leg(1:mTemp)));
    else
        error('Check upwind edge methodology')
    end       
end
subplot(1,2,2)
%Make Lineplot
LatDistance=distance(nanmin(Lat(Leg)),nanmean(Lon(Leg)),nanmax(Lat(Leg)),nanmean(Lon(Leg)));
LonDistance=distance(nanmean(Lat(Leg)),nanmin(Lon(Leg)),nanmean(Lat(Leg)),nanmax(Lon(Leg)));
if LonDistance>LatDistance
    plot(Lon(Leg),CH4(Leg),'color',[0 0 0],'linewidth',2)
    LegLength=length(Leg);
    HalfwayPoint=round(LegLength/2);
    FirstLeg=Leg(1:HalfwayPoint);
    SecondLeg=Leg(HalfwayPoint+1:end)
    %FirstLegCH4=sort(CH4(FirstLeg));
    %FirstLegCH4(isnan(FirstLegCH4)==1)=[];
   % Upwind1BackgroundSpot=round(length(FirstLegCH4).*Percentile);
   % Upwind1Background=FirstLegCH4(Upwind1BackgroundSpot);

    
    %SecondLegCH4=sort(CH4(SecondLeg));
    %SecondLegCH4(isnan(SecondLegCH4)==1)=[];
    %Upwind2BackgroundSpot=round(length(SecondLegCH4).*Percentile);
    %Upwind2Background=SecondLegCH4(Upwind2BackgroundSpot);
% if MeanU>0 & MeanV>0
%         Temp1=Upwind1Background;
%         Temp2=Upwind2Background;
%         Upwind1Background=Temp2;
%         Upwind2Background=Temp1;
% end
    
    if Lon(Leg(1))>Lon(Leg(end))
        Temp1=Upwind1Background;
        Temp2=Upwind2Background;
        Upwind1Background=Temp2;
        Upwind2Background=Temp1;
    end
    
    
    hold on
    Temp=[];
        LineCoeff=polyfit([nanmin(Lon(Leg(1))),nanmax(Lon(Leg(end)))],[Upwind1Background,Upwind2Background],1)
        
    m=LineCoeff(1);
    b=LineCoeff(2);
    BackgroundLine=m.*Lon(Leg)+b;
    BackgroundLineLeg=BackgroundLine(ismember(Leg,Downwind)==1);
    plot(Lon(Leg),BackgroundLine,'Color','blue','LineWidth',2)
    CH4Downwind=CH4(Downwind);
    x=Lon(Downwind);
    y1=CH4Downwind;
    Temp=y1-BackgroundLineLeg;
    BadSpot=find(Temp<0)
    y1(BadSpot)=BackgroundLineLeg(BadSpot);
 %   GS=find(CH4Downwind-BackgroundLine)>0;
 GoodSpot=find(isnan(y1)==0);
 Temp=[x(GoodSpot)' flip(x(GoodSpot))'];
 Temp2=[BackgroundLineLeg(GoodSpot)',flip(y1(GoodSpot))'];
 patch(Temp,Temp2, 'g')
     y1=CH4Downwind;
    Temp=y1-BackgroundLineLeg;
    BadSpot=find(Temp>0)
    y1(BadSpot)=BackgroundLineLeg(BadSpot);
 %   GS=find(CH4Downwind-BackgroundLine)>0;
 GoodSpot=find(isnan(y1)==0);
 Temp=[x(GoodSpot)' flip(x(GoodSpot))'];
 Temp2=[BackgroundLineLeg(GoodSpot)',flip(y1(GoodSpot))'];
 patch(Temp,Temp2, 'r')
  xlabel('Longitude')
 LonIntp=min(Intersection(:,2)):0.001:max(Intersection(:,2));
  TempMatrix=[x(GoodSpot),CH4Downwind(GoodSpot)];
 [junk,UniqueRows]=unique(TempMatrix(:,1));
 TempMatrix=TempMatrix(UniqueRows,:);
 TempMatrix=sortrows(TempMatrix,1);
 
 
DownwindReimagined=interp1(TempMatrix(:,1),TempMatrix(:,2),LonIntp);
DownwindCH4Intp=nanmean(DownwindReimagined);

  TempMatrix=[x(GoodSpot),BackgroundLineLeg(GoodSpot)];
 [junk,UniqueRows]=unique(TempMatrix(:,1));
 TempMatrix=TempMatrix(UniqueRows,:);
 TempMatrix=sortrows(TempMatrix,1);
 UpwindReimagined=interp1(TempMatrix(:,1),TempMatrix(:,2),LonIntp);
UpwindCH4Intp=nanmean(UpwindReimagined);
 
elseif LatDistance>LonDistance
    plot(Lat(Leg),CH4(Leg),'color',[0 0 0],'linewidth',2)
    if (MeanU>0 & MeanV>0) | (MeanU<0 & MeanV<0)
            Temp1=Upwind1Background;
            Temp2=Upwind2Background;
            Upwind1Background=Temp2;
            Upwind2Background=Temp1;
    end
    
    if Lat(Leg(1))>Lat(Leg(end))
        Temp1=Upwind1Background;
        Temp2=Upwind2Background;
        Upwind1Background=Temp2;
        Upwind2Background=Temp1;
    end
    
    
    hold on
    Temp=[];
        LineCoeff=polyfit([nanmin(Lat(Leg(1))),nanmax(Lat(Leg(end)))],[Upwind1Background,Upwind2Background],1)
        

    m=LineCoeff(1);
    b=LineCoeff(2);
    BackgroundLine=m.*Lat(Leg)+b;
    BackgroundLineLeg=BackgroundLine(ismember(Leg,Downwind)==1);
    plot(Lat(Leg),BackgroundLine,'Color','blue','linewidth',2)
    CH4Downwind=CH4(Downwind);
    x=Lat(Downwind);
    y1=CH4Downwind;
    Temp=y1-BackgroundLineLeg;
    BadSpot=find(Temp<0)
    y1(BadSpot)=BackgroundLineLeg(BadSpot);
 %   GS=find(CH4Downwind-BackgroundLine)>0;
 GoodSpot=find(isnan(y1)==0);
 Temp=[x(GoodSpot)' flip(x(GoodSpot))'];
 Temp2=[BackgroundLineLeg(GoodSpot)',flip(y1(GoodSpot))'];
 patch(Temp,Temp2, 'g')
     y1=CH4Downwind;
    Temp=y1-BackgroundLineLeg;
    BadSpot=find(Temp>0)
    y1(BadSpot)=BackgroundLineLeg(BadSpot);
 %   GS=find(CH4Downwind-BackgroundLine)>0;
 GoodSpot=find(isnan(y1)==0);
 Temp=[x(GoodSpot)' flip(x(GoodSpot))'];
 Temp2=[BackgroundLineLeg(GoodSpot)',flip(y1(GoodSpot))'];
 patch(Temp,Temp2, 'r')
  xlabel('Longitude')
 LatIntp=min(Intersection(:,1)):0.001:max(Intersection(:,1));
  TempMatrix=[x(GoodSpot),CH4Downwind(GoodSpot)];
 [junk,UniqueRows]=unique(TempMatrix(:,1));
 TempMatrix=TempMatrix(UniqueRows,:);
 TempMatrix=sortrows(TempMatrix,1);
 
 
DownwindReimagined=interp1(TempMatrix(:,1),TempMatrix(:,2),LatIntp);
DownwindCH4Intp=nanmean(DownwindReimagined);

  TempMatrix=[x(GoodSpot),BackgroundLineLeg(GoodSpot)];
 [junk,UniqueRows]=unique(TempMatrix(:,1));
 TempMatrix=TempMatrix(UniqueRows,:);
 TempMatrix=sortrows(TempMatrix,1);
 UpwindReimagined=interp1(TempMatrix(:,1),TempMatrix(:,2),LatIntp);
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


sgtitle([Date(5:6),'/',Date(7:8),'/',Date(1:4),'  ',Hour,':',Minute,'Z','  (T',num2str(LegLoop),')'],'FontSize',18)
%UpwindMatrix(LegLoop,1:2)=[nanmean(CH4(Upwind1)) nanmean(CH4(Upwind2))];
UpwindMatrix(LegLoop,1)=UpwindCH4Intp;

DownwindMatrix(LegLoop,1)=DownwindCH4Intp;
Distance(LegLoop,1)=distance(Intersection(1,:),Intersection(2,:))*111; 

Enhancement(LegLoop,1)=DownwindMatrix(LegLoop)-UpwindMatrix(LegLoop);
Flux(LegLoop)=Density(LegLoop).*MeanWSPD.*CosTheta(LegLoop).*(Enhancement(LegLoop)/1e6).*(Distance(LegLoop)*1000).*DownwindABL(LegLoop);
annotation('textbox', [0.68, 0.9, 0.15, 0.05], 'String', ['Flux = ', num2str(round(Flux(LegLoop))),' mol/s']);
annotation('rectangle', [0.68, 0.9, 0.15, 0.05],'color',[1 1 1]);
saveas(f,[Date,'/','T',num2str(LegLoop),'-MethodA','.png'])
close(f)   % free memory: hidden figures were never closed
end
clearvars -except UpwindMatrix DownwindMatrix DownwindABL FlightDateLoop FlightList Date LegHeight UpwindMatrix DownwindMatrix Distance Enhancement DownwindTime ABLHeight ABLTheta ABLT ABLPressure ProfileTime MeanWSPD MeanWDIR Density Flux CosTheta
save([Date,'/','Solution-MethodA'])
end









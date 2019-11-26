%% The topography
% load topo
step=1; startLat=41; endLat=45.9; startLon=-109.6; endLon=-101;
startLon=((startLon+180)*60*4)+1;
endLon=((endLon+180)*60*4)+1;
numLon=endLon-startLon;
startLat=((startLat+90)*60*4)+1;
endLat=((endLat+90)*60*4)+1;
numLat=endLat-startLat;
basePath= 'E:\topoData\';
elon=ncread([basePath 'topo15_compressed.nc'],'lon',startLon,numLon,step);
elat=ncread([basePath 'topo15_compressed.nc'],'lat',startLat,numLat,step);
elev=ncread([basePath 'topo15_compressed.nc'],'z',[startLon startLat],[numLon numLat],[step step]);
elev=elev';

disp('topo read')


%% The tStar 
% load tStar
WYRes = csvread('WY.txt');

mLat=WYRes(:,2);
mLon=WYRes(:,1);
mTstar=WYRes(:,3);

% load MAGIC_Deep.txt
% 
% mLat=MAGIC_Deep(:,2);
% mLon=MAGIC_Deep(:,1);
% mTstar=MAGIC_Deep(:,3);

% load MAGIC_Deep_slopes.txt
% 
% mLat=MAGIC_Deep_slopes(:,2);
% mLon=MAGIC_Deep_slopes(:,1);
% mTstar=MAGIC_Deep_slopes(:,3);

% make interpolant
F=scatteredInterpolant(mLon,mLat,mTstar,'nearest','none');

disp('t* interpolant made')

%% Beautiful map
[eLon eLat]=meshgrid(elon,elat);
eLon=eLon; eLat=eLat;
eTS=F(eLon,eLat); %basin depth interpolated onto topo points
disp('done interpolating tStar onto topo')

load CMfine
m=size(cm,1); TSmin=min(eTS(:))-0.00000001; TSmax=max(eTS(:))+0.00000001;
eTSind = ceil((eTS-TSmin)/(TSmax-TSmin)*m);

C=ind2rgb(eTSind,cm); 
BdRck=isnan(eTS);
BdRck=repmat(BdRck,[1 1 3]);

C(BdRck)=0.6; %no obs. BdRck is an inherited varname

%below sea level
oce=elev<0;
oce=repmat(oce,[1 1 3]);
C(oce)=0.1;

%in great lake
% load conus
% gtlk=inpolygon(eLon,eLat,gtlakelon,gtlakelat);
% gtlk=repmat(gtlk,[1 1 3]);
% C(gtlk)=0.1;



% C=permute(C,[2 1 3]); %this is like a generalized transpose

disp('making figure')

figure
surface(elon,elat,elev/10000,C,'EdgeColor','none')
daspect([1 cosd(40) 1])
xlabel('longitude (°)')
ylabel('latitude (°)')
material dull
light('position',[-88 47 150])
axis tight
box on
set(gca,'fontsize',14)
hold on

colormap(cm);
hcb=colorbar('EastOutside'); caxis([TSmin TSmax])
hcb.Label.String='\Deltat* (s)';
hcb.Label.FontSize=14;

set(gcf,'position',[400 264 700 540])

hold on
scatter3(lon_sta, lat_sta, ones(size(lon_sta)).*max(max(elev))/10000, 10, 'blue');

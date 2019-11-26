%% Pretty plots cont'd from smoothTSinversion.m
for i=1:size(MImodel,1)
    for j=1:size(MImodel, 2)
        [latMat(i,j), lonMat(i,j)] = minvtran(mstruct,xMat(i,j),yMat(i,j));
    end
end

figure(1);
% contourf(xMat,yMat,MImodel, v);
contourf(lonMat, latMat, MImodel, v);
h = colorbar; %caxis([0 2])
h.Label.String = '\Deltat*'; %delta t*
load('CMfine');
colormap(cm)
caxis(clim)
hold on

% [dataX, dataY] = mfwdtran(mstruct,lat_sta,lon_sta);
scatter(lon_sta,lat_sta,10, MIsta,'filled', 'MarkerEdgeColor', 'k')
set(gca,'YDir','normal')
daspect([1 1 1])

% C=load('coast');
%  
% [C.x,C.y]=mfwdtran(mstruct,C.lat,C.long);
% C.x(C.x>maxX)=nan;
% C.x(C.x<minX)=nan;
% C.y(C.y>maxY)=nan;
% C.y(C.y<minY)=nan;
% plot(C.lat,C.long,'k-','LineWidth',2)

xlim([ minX maxX ]);
ylim([ minY maxY ]);

ylabel('lat');
xlabel('lon');

figure(2)
histogram(MIsta)
xlabel('Size of the station terms');

save(name, 'MImodel', 'xMat', 'yMat', 'mstruct', 'MIevt', 'MIsta', 'dataE', 'dataI', 'dataV', 'dataX', 'dataY', 'uSta', 'lon_sta', 'lat_sta');
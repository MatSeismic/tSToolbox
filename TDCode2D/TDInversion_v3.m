%fit a smooth surface to all t* observations through an inversion scheme.
clear, close all

name = 'CIELOTimeDomain.pdf';

TD_parameters.sig_b        = 0.01;%in s
TD_parameters.sig_sig      = 0.01;%in s
TD_parameters.sig_flag     = 1;%1 is uniform, 2 is station, 3 is event
TD_parameters.sig_r        = 0.05;%in km
TD_parameters.max_nodes    = 100;
TD_parameters.min_error    = 0.001;
TD_parameters.interp_style = 'nearest';
    TD_parameters.burn_in        = 5e3;
    TD_parameters.keep_each      = 100;

TD_parameters.ydamp          = 0.01;%in km
TD_parameters.range          = [ -0.5 0.5 ];
TD_parameters.n_iter         = 1e5;%for now, just iterate to a max
TD_parameters.max_cells      = 100;%for starting
TD_parameters.min_cells      = 10;%CANNOT be below 4
TD_parameters.keep_bestN     = 0.05;%fraction

nodeSpacing    = 0.1;%for plotting
buffr          = 0.5;%in km
rotation       = 0;%in degrees
n_chains       = 10;
noise          = 0.075;%if a synthetic test

pdf_x = -0.4:0.025:0.4;

allLats   = [];
allLons   = [];
allTS     = [];
allSig    = [];
dataE     = [];
allSta    = {};

%fnames = dir('../DeepResults/*result.mat');
fnames = dir('../../Mv03/*Measurement.mat');
%fnames = dir('./*_ShallowWithError.mat');
fnames = { fnames.name };

for k=1:length(fnames)

    fname=[ '../../Mv03/' fnames{k} ];
%    fname=[ '../DeepResults/' fnames{k} ];

    load(fname, 'ts_run', 'Traces');
    traces = Traces; clear Traces
            
    sta   = {traces.station};
    lats  = [traces.latitude];
    lons  = [traces.longitude];
    tS    = [ts_run.tStar_WF];
    
    tS=tS-mean(tS);
    
    allLats  = [ allLats lats ];
    allLons  = [ allLons lons ];
    allTS    = [ allTS tS ];
    allSta   = [ allSta sta ];
    dataE    = [ dataE lats*0+k ]; % this is just the event number
    
end

uSta=unique(allSta); %list of all station names

%%%%synthetic data
%synthetic model with two steps
%allTS(:) = normrnd(0,noise,size(allTS));

% allTS(:)                                   = 0;%normrnd(0, noise, size(allTS));
% allTS( allLons < -79 & allLons > -80 )     =  0.2;
% allTS( allLons < -81 & allLons > -82 )     = -0.2;
% 
% for k = 1:3
% 
% 	allTS(dataE == k) = allTS(dataE==k) + normrnd(0,0.075, size(allTS(dataE==k)));
% 
% end
% 
% for k = 4:6
% 
% 	allTS(dataE == k) = allTS(dataE==k) + normrnd(0, 0.15, size(allTS(dataE==k)));
% 
% end

% 
% trueModel = zeros(size(allLons));
% trueModel( allLons < -81 & allLons > -82 )  = -0.2;
% trueModel( allLons < -79 & allLons > -80 )      = 0.2;

%synthetic model with a linear increase
%allLons = linspace(min(allLons), max(allLons), length(allLons));
%allLats = linspace(min(allLats), max(allLats), length(allLats));
%allTS(:) = normrnd(0, noise, size(allTS)) + (( 0.4/(max(allLons) - min(allLons)) *(allLons - min(allLons))) - 0.2);

%make half the data random
%allTS( (length(allTS)+1):2*length(allTS) ) = 0.8*rand(size(allTS)) - 0.4;
% dataE   = [ dataE ( dataE + max(dataE)) ];
% allLats = [ allLats allLats ];
% allLons = [ allLons allLons ];
% allSta  = [ allSta allSta ];

for k = 1:length(unique(dataE))
   
    allTS(dataE==k) = allTS(dataE==k) - mean(allTS(dataE==k));
    
end

%% define a map projection and inversion grid
%origin in the center of the array footprint
parpool;

centerLat = min(allLats)+(max(allLats)-min(allLats))/2;
centerLon = min(allLons)+(max(allLons)-min(allLons))/2;

origin               = [ centerLat centerLon ];           % array center [ lat lon ]
mstruct              = defaultm('mercator');
mstruct.origin       = [ origin rotation ];%second number is rotation
mstruct              = defaultm( mstruct );
mstruct.scalefactor  = 6371;

%%%%%%%%%%
%%%consider rotating the grid for the inversion so cut down the number
%%%of unsed nodes
%%%%%%%%%%

%projected lat/lon of observations to X/Y:
%[dataX, dataY] = mfwdtran(mstruct,allLats,allLons);
dataX= allLons;
dataY= allLats;

%make the grid on which to invert

minX=min(dataX)-buffr;
maxX=max(dataX)+buffr;
minY=min(dataY)-buffr;
maxY=max(dataY)+buffr;

xVec=minX:nodeSpacing:maxX;
yVec=minY:nodeSpacing:maxY;

[xMat, yMat]=meshgrid(xVec,yVec);
    
dataStruct.allTS = allTS;
dataStruct.allLats = allLats;
dataStruct.allLons = allLons;
dataStruct.allSig = 0.4*ones(size(allTS));
dataStruct.allSta = allSta;
dataStruct.dataX = dataX;
dataStruct.dataY = dataY;
dataStruct.dataE = dataE;

%% Let's invert
parfor k = 1:n_chains

    warning off MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId
    
    %disp(['Chain #' num2str(k) ]);
    pipi=k;
    
    bestmodels(:, k) = TD_inversion_function(xVec, yVec, TD_parameters, dataStruct,pipi);
        
end

bestmodels = bestmodels(:);
[ model_stats, pdf_set_plot ] = make_models( ...
    bestmodels, xMat, yMat, TD_parameters.sig_flag, pdf_x, 0.66, TD_parameters.interp_style);

%apply event statics and plot
%es = e_static(allTS, dataE);

%for k = 1:length(allTS)
   
 %   allTS(k) = allTS(k) + es(dataE(k));
    
%end

disp(['over']);

colorvec = get(groot, 'defaultAxesColorOrder');

figure
contourf(xMat, yMat, model_stats.mean,10)
hold on
plot(dataX, dataY, 'rs')
colorbar

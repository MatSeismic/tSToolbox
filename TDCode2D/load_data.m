function [ dataStruct ] = load_data( TD_parameters )

    allLats   = [];
    allLons   = [];
    allTS     = [];
    allSig    = [];
    dataE     = [];
    allSta    = {};
    
    fnames = dir('C:\Users\Zhao\Google Drive\dropbox_bkp\Attenutation\Mv03\*Measurement.mat');
    fnames = { fnames.name };
    
    for k=1:length(fnames)

        fname=[ 'C:\Users\Zhao\Google Drive\dropbox_bkp\Attenutation\Mv03\' fnames{k} ];

        load(fname, 'ts_run', 'Traces');

        use = [Traces.QC];
        ts_run = ts_run(use);

        %loop to remove doubled lats/lons (not my preferred solution)
        for kt=1:length(ts_run)
            lt=ts_run(kt).latitude;
            ts_run(kt).latitude=lt(1);
            ln=ts_run(kt).longitude;
            ts_run(kt).longitude=ln(1);  
        end

        sta   = {ts_run.station};
        lats  = [ts_run.latitude];
        lons  = [ts_run.longitude];
        tS    = [ts_run.tStar_WF];
        sig   = 0.3*ones(size([ts_run.tStar_WF]));%[traces.sigma];

        allLats  = [ allLats lats ];
        allLons  = [ allLons lons ];
        allTS    = [ allTS tS ];
        allSig   = [ allSig sig ];
        allSta   = [ allSta sta ];
        dataE    = [ dataE lats*0+k ]; % this is just the event number

    end

    uSta=unique(allSta); %list of all station names

    for k = 1:length(unique(dataE))
        allTS(dataE==k) = allTS(dataE==k) - mean(allTS(dataE==k));
    end

    %% define a map projection and inversion grid
    %origin in the center of the array footprint

    centerLat = min(allLats)+(max(allLats)-min(allLats))/2;
    centerLon = min(allLons)+(max(allLons)-min(allLons))/2;

    origin               = [ centerLat centerLon ];           % array center [ lat lon ]
    mstruct              = defaultm('mercator');
    mstruct.origin       = [ origin TD_parameters.rotation ];%second number is rotation
    mstruct              = defaultm( mstruct );
    mstruct.scalefactor  = 6371;

    %%%%%%%%%%
    %%%consider rotating the grid for the inversion so cut down the number
    %%%of unsed nodes
    %%%%%%%%%%

    %projected lat/lon of observations to X/Y:
    [dataX, dataY] = mfwdtran(mstruct,allLats,allLons);
%     dataY = zeros(size(dataY)); %???? why

    %make the grid on which to invert

    minX = min(dataX) - TD_parameters.buffer;
    maxX = max(dataX) + TD_parameters.buffer;
    minY = min(dataY) - TD_parameters.buffer;
    maxY = max(dataY) + TD_parameters.buffer;

    xVec = minX:TD_parameters.nodeSpacing:maxX;
    yVec = minY:TD_parameters.nodeSpacing:maxY;

    dataStruct.allTS   = allTS;
    dataStruct.allLats = allLats;
    dataStruct.allLons = allLons;
    dataStruct.allSig  = allSig;
    dataStruct.allSta  = allSta;
    dataStruct.dataX   = dataX;
    dataStruct.dataY   = dataY;
    dataStruct.dataE   = dataE;
    dataStruct.xVec    = xVec;
    dataStruct.yVec    = yVec;
    dataStruct.mstruct = mstruct; %% for the purpose of plotting results on topo map
    
end


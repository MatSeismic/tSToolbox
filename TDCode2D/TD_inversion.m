%fit a smooth surface to all t* observations through an inversion scheme.
clear, close all

%parpool(12);
parpool;

name = 'CIELO';

TD_parameters = define_TDstructure( );

dataStruct  = load_data(TD_parameters);

TD_parameters.interp_style = 'nearest';
TD_parameters.sig_flag     = 3;


%% Let's invert
parfor k = 1:TD_parameters.n_chains
    
    disp(['Chain #' num2str(k) ]);
    warning off MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId
    models(:, k) = TD_inversion_function(TD_parameters, dataStruct);
        
end

if ~exist([name 'inversion1.mat'], 'file')
    
    savename = [ name 'inversion1.mat'];

else
    
    files = dir([name 'inversion*.mat']);
    n = length(files) + 1;
    
    savename = [ name 'inversion' num2str(n) '.mat'];
    
end

disp(['Saving models named ' savename]);
save(savename)

models          = models(:);
[xMat, yMat]    = meshgrid(dataStruct.xVec, dataStruct.yVec);
[ model_stats ] = make_models( models, xMat, yMat, TD_parameters.sig_flag, -0.5:0.05:0.5, 0.66, TD_parameters.interp_style);

figure
contourf(xMat, yMat, model_stats.mean, -0.2:0.02:0.2);
h = colorbar; h.Label.String = '\Deltat*';
colormap('pink')
hold on
plot(dataStruct.dataX, dataStruct.dataY, 'rs');
title('\Deltat* variation');
saveas(gcf, [name 'inversion' num2str(n) 'tstar'],'png');

figure
contourf(xMat, yMat, model_stats.std);
h = colorbar; h.Label.String = 'std';
hold on
plot(dataStruct.dataX, dataStruct.dataY, 'rs');
title('std variation');
saveas(gcf, [name 'inversion' num2str(n) 'std'],'png');

[X, P] = meshgrid(xMat(1, :), -0.5:0.05:0.5);

figure
contourf(X', P', squeeze(model_stats.pdf_set(5,:,:)), 10)
colorbar
colormap(flipud(pink))

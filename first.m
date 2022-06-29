% This script estimates quantile regressions of future GDP on current
% values of GDP and the NFCI, and matches skewed t-distributions to the 
% estimated quantiles. 
% 
% The original code/script comes from the replication files of:
% Tobias Adrian, Nina Boyarchenko, and Domenico Giannone (2019):
% "Vulnerable Growth," American Economic Review.
% 
% Any errors or omissions arising from the modification of the original script  
% of the above paper should not be attributed to the above authors

% This script requires MATLAB's Optimization Toolbox.

%% Clear workspace, set file paths and graphics/estimation settings
clear
close all
clc

addpath('azzalini')

set(0, 'defaultAxesFontName', 'Times');
set(0, 'DefaultAxesFontSize',15)
set(0, 'defaultAxesLineStyleOrder', '-|--|:', 'defaultLineLineWidth', 1.5)
setappdata(0, 'defaultAxesXTickFontSize', 1)
setappdata(0, 'defaultAxesYTickFontSize', 1)

% Folder to store figures
FigSubFolder = 'FigMainResults';
if ~exist(FigSubFolder,'dir')
    mkdir(FigSubFolder);
end

% Should saved density matching results be loaded?
% Or should density matching be reperformed?
loadSavedResults = false;

%% Load data and fix forecast settings
load Backup_Backup;
period_numeric = datenum(period_array);
% Use 1973Q1-2015Q4 subsample
jtFirst = 1;    % 2005Q1
jtLast  = 209;  % 2022Q2
period_numeric = period_numeric(jtFirst:jtLast);
main_array = main_array(jtFirst:jtLast, :);
clear('jtFirst', 'jtLast')
[T, n] = size(main_array);

% Forecast settings
H = [1, 4];          % Horizons to forecast (# of quarters ahead)
QQ = 0.05:0.05:0.95; % Quantiles to estimate in quantile regressions
% Indices for selected quantiles
[~, jq05] = min(abs(QQ - 0.05));
[~, jq25] = min(abs(QQ - 0.25));
[~, jq50] = min(abs(QQ - 0.50));
[~, jq75] = min(abs(QQ - 0.75));
[~, jq95] = min(abs(QQ - 0.95));

% Construct average growth rates
y = main_array(:, strcmp(Mnem, 'A191RL1Q225SBEA'));
for h = 1:4
    Yh(:, h) = filter(ones(1, h)/h, 1, y);
    Yh(1:(h - 1), h) = NaN;
end


%% Main results
% Estimate quantile regressions (for all quantiles and forecast horizons)
ResMain    = QRboot(main_array(:, strcmp(Mnem, 'KK')), y, H, 1, QQ);
ResGDPonly = QRboot([], y, H, 1, QQ);
ResUnc     = QRboot([], y, H, 0, QQ);



% Loop over forecast horizons
for h = H
    % Get quantiles of predictive distribution
    YQunc     = ResUnc.YQ(end, :, h);
    YQ        = ResMain.YQ(:, :, h);
    Y2        = ResMain.Y2(:, h);
    YQGDPonly = ResGDPonly.YQ(:, :, h);
    QQ        = ResMain.QQ;
    
        
    %% Fit skewed-t densities to estimated quantiles OR load saved results
    if loadSavedResults
        filename = ['ResMatch_H', num2str(h), '.mat'];
        disp(['Loading saved density matching results from file ', filename])
        load(filename)
        clear('filename')
    else
        disp('Matching densities for quantile regression using GDP and KK.')
        ResMatch = Step2match(YQ, YQunc, QQ);
        disp('Matching densities for quantile regression using GDP only.')
        ResMatchGDPonly = Step2match(YQGDPonly, YQunc, QQ);

        % Save results to .mat file
        filename = ['ResMatch_G', num2str(h), '.mat'];
        disp(['Saving results to file ', filename])
        save(filename, 'ResMatch', 'ResMatchGDPonly')
        clear('filename')
    end
       %% 3D density plot
    f = figure;
    meshc(period_numeric, ResMatch.YY', ResMatch.PST')
    datetick('x', 'yyyy')
    view(85, 50)
    set(gca, 'YLim', [-20, 20])
    set(gca, 'XLim', [period_numeric(1), period_numeric(end)])
    xlabel('Year')
    filename = fullfile(FigSubFolder, ['Dens3D_G', num2str(h), '.pdf']);
    printpdf(f, filename);
    clear('f', 'filename')
end


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
loadSavedResults = true;

%% Load data and fix forecast settings
load DataVulnerability;
% Use 1973Q1-2015Q4 subsample
jtFirst = 1;    % 1973Q1
jtLast  = 911;  % 2015Q4
Time = Time(jtFirst:jtLast);
X = X(jtFirst:jtLast, :);
clear('jtFirst', 'jtLast')
[T, n] = size(X);

% Forecast settings
H = [1, 4];          % Horizons to forecast (# of quarters ahead) CHANGE HERE
QQ = 0.05:0.05:0.95; % Quantiles to estimate in quantile regressions
% Indices for selected quantiles
[~, jq05] = min(abs(QQ - 0.05));
[~, jq25] = min(abs(QQ - 0.25));
[~, jq50] = min(abs(QQ - 0.50));
[~, jq75] = min(abs(QQ - 0.75));
[~, jq95] = min(abs(QQ - 0.95));

% Construct average growth rates
y = X(:, strcmp(Mnem, 'A191RL1Q225SBEA'));
for h = 1:4
    Yh(:, h) = filter(ones(1, h)/h, 1, y);
    Yh(1:(h - 1), h) = NaN;
end

%Dates for which you want to plot the density
JT = [134; 144; 168;];  % 2006Q2; 2008Q4; 2014Q4
    
% axis limits for various plots (these are horizon-dependent)
ylimsPredictedDistribution = [
    -15,  20;
	NaN, NaN;
	NaN, NaN;
	-15,  20;
    ];

%% Main results
% Estimate quantile regressions (for all quantiles and forecast horizons)
ResMain    = QRboot(X(:, strcmp(Mnem, 'KK')), y, H, 1, QQ);
ResGDPonly = QRboot([], y, H, 1, QQ);
ResUnc     = QRboot([], y, H, 0, QQ);


% Loop over forecast horizons
for h = 1
    % Get quantiles of predictive distribution
    YQunc     = ResUnc.YQ(end, :, h);
    YQ        = ResMain.YQ(:, :, h);
    Y2        = ResMain.Y2(:, h);
    YQGDPonly = ResGDPonly.YQ(:, :, h);
    QQ        = ResMain.QQ;
    
        
    %% Fit skewed-t densities to estimated quantiles OR load saved results
    if loadSavedResults
        filename = ['ResMatch_I', num2str(h), '.mat'];
        disp(['Loading saved density matching results from file ', filename])
        load(filename)
        clear('filename')
    else
        disp('Matching densities for quantile regression using GDP and KK.')
        ResMatch = Step2match(YQ, YQunc, QQ);
        disp('Matching densities for quantile regression using GDP only.')
        ResMatchGDPonly = Step2match(YQGDPonly, YQunc, QQ);

        % Save results to .mat file
        filename = ['ResMatch_I', num2str(h), '.mat'];
        disp(['Saving results to file ', filename])
        save(filename, 'ResMatch', 'ResMatchGDPonly')
        clear('filename')
    end
    %% 3D density plot
    f = figure;
    mesh(Time, ResMatch.YY', ResMatch.PST')
    datetick('x', 'yyyy')
    view(85, 50)
    set(gca, 'YLim', [-20, 20])
    set(gca, 'XLim', [Time(1), Time(end)])
    
    zlabel('Probability')
    ylabel("Net Portfolio Inflows (USD bn)")
    filename = fullfile(FigSubFolder, ['Dens3D_H', num2str(h), '.pdf']);
    printpdf(f, filename);
    clear('f', 'filename')
end
    
   

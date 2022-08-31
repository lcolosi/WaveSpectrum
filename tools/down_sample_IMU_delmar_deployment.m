%% Down Sampling IMU data from Deployment off Delmar
% Luke Colosi | lcolosi@ucsd.edu | February 25th, 2021

clc, clear, close all

% Set path to data
%%%% Vehicle %%%%
vehicle = 'PLANCK';
%vehicle = 'STOKES';

%%%% Root %%%%
ROOT = '../data/DELMAR2020/';

% Load data
load([ROOT vehicle '/20200908_' vehicle '_Del_Mar_IMU.mat']) 

%% Down sample Data

% Set temporal parameters:
t0 = datenum('000','FFF');
tf = datenum('050','FFF');

% Set time steps for 20 Hz time series
t_step = tf- t0;

% Create edges for 20hz time series 
time_20hz_edges = NOV.time(1):t_step:NOV.time(end);

% Obtain the indices for the bins (ranges that will be bin averaged) and
% number of bins 
[~,~,loc]=histcounts(NOV.time(1:end-5), time_20hz_edges);
nbins = numel(unique(loc));

% Create field name char variable
field = fieldnames(NOV);

% Loop through field names (except time) 
for k = 2:length(field)

    % Set field variable name
    name = char(field(k,:));

    % Compute the mean with each time step bin for each variable omitting NaNs
    eval(['nov.' name ' = accumarray(loc(:), NOV.' name '(1:end-5),[nbins,1], @(x)mean(x,''omitnan''));']);

end 

% Obtain the midpoint time step corresponding to each bin
nov.time_20hz = (0.5*(time_20hz_edges(1:end-1)+time_20hz_edges(2:end)))';

% Set path directory for saved data
path_out = '../data/DELMAR2020/';

% Save downsampled data
save([path_out vehicle '/NOVATEL_downsampled_20Hz_ALL_' vehicle '.mat'],'nov', '-v7.3')

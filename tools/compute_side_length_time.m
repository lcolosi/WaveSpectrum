%% Compute the length and time duration of square trajectory sides 
% Luke Colosi | lcolosi@ucsd.edu | August 20th, 2022 

clc, clear, close all;

% Set the current working directory
cd ~/Desktop/projects/asi_lab_summer_internship/WaveGlider/src

% Set path to data
%%%% Vehicle %%%%.
vehicle = 'STOKES';

%%%% Root %%%%
ROOT = '../data/DELMAR2020/';

%%%%%%%%%% Initial global variables %%%%%%%%%%

% Physical parameters
h = 0;                                                                      % altitude
SPHEROID = referenceEllipsoid('wgs84', 'm');                                % Reference ellipsoid
dir_con = {'CW', 'gt', 'rn'};                                               % Directional conventions

% Temporal parameters
period = 10;                                                                % Time duration for computing wave spectra
date_i = '09-Sep-2020 02:30:00';                                            % Start date 
date_f = '11-Sep-2020 16:10:00';                                            % End date; 
date_o = {date_i, date_f}; 
%------- Novatel ------- %
fe_n = 20;                                                                  % Sampling rate (Hz)
dt_n = 1/fe_n;                                                              % Time interval between measurements (s)
%------- Weather Station ------- %
fe_w = 1; 
dt_w = 1/fe_w; 

% Upload and process novatel and weather station data
[nov_s, ~, nlegs_s] = process_wg_data(vehicle, ROOT, date_o, dir_con, [dt_n, dt_w]);

%% Compute the distance in meters of the legs in the large and small squares 

%------------ Large Square ------------%  

%%%%%% Horizontal leg %%%%%%
[N, E]   = geodetic2ned(nov_s.L3.latitude(1), nov_s.L3.longitude(1), h, nov_s.L3.latitude(end), nov_s.L3.longitude(end), h, SPHEROID);
delta_x_lsh = norm([N, E]); % Units: Meters

%%%%%% Vertical leg %%%%%%
[N, E]   = geodetic2ned(nov_s.L4.latitude(1), nov_s.L4.longitude(1), h, nov_s.L4.latitude(end), nov_s.L3.longitude(end), h, SPHEROID);
delta_x_lsv = norm([N, E]); % Units: Meters

%------------ Small Square ------------%  

%%%%%% Horizontal leg %%%%%%
[N, E]   = geodetic2ned(nov_s.L70.latitude(1), nov_s.L70.longitude(1), h, nov_s.L70.latitude(end), nov_s.L70.longitude(end), h, SPHEROID);
delta_x_ssh = norm([N, E]); % Units: Meters  

%%%%%% Vertical leg %%%%%%
[N, E]   = geodetic2ned(nov_s.L71.latitude(1), nov_s.L71.longitude(1), h, nov_s.L71.latitude(end), nov_s.L71.longitude(end), h, SPHEROID);
delta_x_ssv = norm([N, E]); % Units: Meters

% Display results
disp('Large Square')
disp(['Horizontal leg: ' num2str(delta_x_lsh), ' m'])
disp(['Vertical leg: ' num2str(delta_x_lsv), ' m'])
disp(' ')
disp('Small Square')
disp(['Horizontal leg: ' num2str(delta_x_ssh), ' m'])
disp(['Vertical leg: ' num2str(delta_x_ssv), ' m'])

%% Compute the average time in minutes to complete lags in the large and small box tracks

% Loop through legs
for n = 1:nlegs_s

    % Compute time duration for ith leg
    eval(['nov_s.time_legs(n) = (nov_s.L' num2str(n) '.time_20hz(end) - nov_s.L' num2str(n) '.time_20hz(1))*(24*60);'])
     
    % Compute mean latitude
    eval(['nov_s.mlat_legs(n) = mean(nov_s.L' num2str(n) '.latitude);'])

end

% Obtain latitude indicies for large and small box cutoff
Ilat_l = find(nov_s.mlat_legs > 32.93);                                     % Indices for legs in large box
Ilat_s = find(nov_s.mlat_legs < 32.93);                                     % Indices for legs in small box

% Compute mean time duration of large and small squares

%------------ Large Square ------------%
mtime_duration_l = mean(nov_s.time_legs(Ilat_l));

%------------ Small Square ------------%
mtime_duration_s = mean(nov_s.time_legs(Ilat_s));

% Display results
disp(['Large Square: ' num2str(mtime_duration_l) ' minutes'])
disp(' ')
disp(['Small Square: ' num2str(mtime_duration_s) ' minutes'])


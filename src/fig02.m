%% Figure 2: Overview of experiment site/assets and environmental conditions for SMODE2021
% Luke Colosi | lcolosi@ucsd.edu | August 31st, 2022

%-------------------------------- Caption --------------------------------%
% (a) Trajectory of Wave Glider WHOI43 (blue) during the SMODE2021 
% experiment for the time period of 29 October at 0:00 UTC to 3 
% November at 23:00 UTC. (b) The geographic region with the location of
% the experiment site (white box). (c) Mean platform heading, (d) wind
% speed at 10 meters above the ocean surface (solid line) and wind 
% direction at 1 meter above the ocean surface (triangular markers), and
% (e) significant wave height measured by WHOI43 (blue).
%-------------------------------------------------------------------------%

clc, clearvars -except N W bathy_g bathy_sc bathy_nc, close all;

% Set text interpreter 
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex')

% Set path to data
%%%% Vehicle %%%%.
vehicle = 'WHOI43';

%%%% Root %%%%
ROOT = '../data/SMODE2021/';
ROOT_bathy = '../data/BATHY/';

% Set path for figures
fig_path = '../figs/';

%%%%%%%%%% Initial global variables %%%%%%%%%%

% Physical parameters
z_station = 1;                                                              % Height of weather station (meters)
g = 9.81;                                                                   % gravitational acceleration (m/s) 
dir_con = {'CW', 'gt', 'rn'};                                               % Directional conventions for wave glider trajectory (wind direction is CW, coming from, reference north)

% Temporal parameters
period = 20;                                                                % Time duration for computing wave spectra
date_i = '29-Oct-2021 00:00:00';                                            % Start date 
date_f = '04-Nov-2021 00:00:00';                                            % End date;
%------- Novatel ------- %
fe_n = 20;                                                                  % Sampling rate (Hz)
dt_n = 1/fe_n;                                                              % Time interval between measurements (s)
%------- Weather Station ------- %
fe_w = 1; 
dt_w = 1/fe_w; 

% Spectral parameters                                                      
fn = (1/2)*fe_n;                                                            % Nyquist frequency (Hz)
df = 0.01;                                                                  % Frequency resolution
dtheta = 5;                                                                 % Azimuthal resolution 
nfft = fe_n/df;                                                             % Window Length
f = (0:df:fn);                                                              % Observed frequency 
lambda_c = 1.5613;                                                          % Wavelength cutoff (meters): corresponding to approximately half the length of the SV3 wave glider.
f_noise = sqrt(g/(2*pi*lambda_c));                                          % Noise frequency cutoff 
freq_band = find(f >= 0.02 & f <= f_noise);                                 % Frequency band (for computing significant wave height)
toolbox = 'WAFO';                                                           % Method used to compute directional spectrum 
variables = 'heave_velocity';                                               % Heave and horizontal velocity are used to compute the direction spectrum
scaling = false;                                                            % Variance of directional spectrum is not scaled to match variance of heave spectrum 

%% Call Data

%--------- WHOI43 ---------%
if isempty(whos('N')) || isempty(whos('W'))
    N = load([ROOT vehicle '/SV3-1043_Vectornav.mat']);                     % RAW GPS/IMU Data from VectorNav (20 Hz)
    W = load([ROOT vehicle '/PLD2_1hz_ALL.mat']);                           % WXT/GILL Weather Station Data (1 hz)
    
    % Convert velocity down to velocity up via a sign change 
    N.vu = -N.vd; 
    
    % Set the time_20Hz data variable as the CR6 time
    N.time = N.time.cr6; 

end  

%--------- Bathymetry ---------%
if isempty(whos('bathy_g')) || isempty(whos('bathy_nc')) || isempty(whos('bathy_sc'))

    %--------- Global ---------%
    bathy_g.z = ncread([ROOT_bathy 'etopo1.nc'], 'z');
    bathy_g.lon = ncread([ROOT_bathy 'etopo1.nc'], 'lon');
    bathy_g.lat = ncread([ROOT_bathy 'etopo1.nc'], 'lat');

    %--------- Northern California ---------%
    bathy_nc.z = ncread([ROOT_bathy 'crm_crm_vol7.nc'], 'z');
    bathy_nc.lon = ncread([ROOT_bathy 'crm_crm_vol7.nc'], 'x');
    bathy_nc.lat = ncread([ROOT_bathy 'crm_crm_vol7.nc'], 'y');

    %--------- Southern California ---------%
    bathy_sc.z = ncread([ROOT_bathy 'crm_crm_vol6.nc'], 'z'); 
    bathy_sc.lon = ncread([ROOT_bathy 'crm_crm_vol6.nc'], 'x');
    bathy_sc.lat = ncread([ROOT_bathy 'crm_crm_vol6.nc'], 'y');
end

% Find overlapping indices for lonitude
Ilon_nc = bathy_nc.lon >= -126 & bathy_nc.lon <= -117; 
Ilon_sc = bathy_sc.lon >= -126 & bathy_sc.lon <= -117; 

% Crop local bathymetry 
bathy_nc.z_crop = bathy_nc.z(Ilon_nc,:);
bathy_sc.z_crop = bathy_sc.z(Ilon_sc,:);

% Merge bathymetry data from northern and southern california 
lon = bathy_sc.lon; 
lat = cat(1,bathy_sc.lat(1:end-1), bathy_nc.lat);
z = cat(2, bathy_sc.z_crop, bathy_nc.z_crop); 

% Change the range of longitude for global bathymetry data and shift
% bathymetry
Ilon_p = bathy_g.lon <= 180;  Ilon_n = bathy_g.lon > 180; 
bathy_g.z_shift = cat(1,bathy_g.z(Ilon_n,:), bathy_g.z(Ilon_p,:));
bathy_g.lon_shift = bathy_g.lon - 180;

% Crop global latitude and longitude 
Ilon_g =  bathy_g.lon_shift >= -126 & bathy_g.lon_shift <= -117; Ilat_g =  bathy_g.lat >= 32 & bathy_g.lat <= 44;
lon_g = bathy_g.lon_shift(Ilon_g);
lat_g = bathy_g.lat(Ilat_g);
z_g = bathy_g.z_shift(Ilon_g,Ilat_g);

%% Index Wave Glider trajectory from time period and bathymetry for experiment site 

% Find temporal indices for the experimental time period
time_idx = find(N.time >= datenum(date_i) & N.time < datenum(date_f));

% Apply to wave glider latitude and longitude trajectory
N.latitude_n = N.lat(time_idx);
N.longitude_n = N.lon(time_idx); 

% Initialize spatial parameters
lon_exp_l = -124.2; lon_exp_h = -123.4;
lat_exp_l = 36.8; lat_exp_h = 37.4;
lon_reg_l = -125.5; lon_reg_h = -121.5;
lat_reg_l = 36; lat_reg_h = 39;

% Find spatial indices for the experimental site
Ilon_exp = lon >=  lon_exp_l & lon <= lon_exp_h; Ilat_exp = lat >=  lat_exp_l & lat <= lat_exp_h; 
Ilon_reg = lon_g >=  lon_reg_l & lon_g <= lon_reg_h; Ilat_reg = lat_g >=  lat_reg_l & lat_g <= lat_reg_h; 

% Apply to bathymetry, longitude, and latitude data 
lon_exp = lon(Ilon_exp); lat_exp = lat(Ilat_exp); z_exp = z(Ilon_exp,Ilat_exp); 
lon_reg = lon_g(Ilon_reg); lat_reg = lat_g(Ilat_reg); z_reg = z_g(Ilon_reg,Ilat_reg);

% Set NaNs to max depth value 
z_exp(isnan(z_exp)) = -3800;     

%% Process data

%--- Index data for specified time period and interpolating over NaNs ---%

% Set variables
field_20hz = {'alt'; 've'; 'vn'; 'vu'; 'lon'; 'lat'; 'time'};
field_1hz = {'TWD'; 'TWS'; 'time'};

% Find indices for desired time frame
time_20hz_idx = find(N.time >= datenum(date_i) & N.time < datenum(date_f));
time_1hz_idx = find(W.PLD2_1hz.time >= datenum(date_i) & W.PLD2_1hz.time < datenum(date_f));

% Index time variables for interpolation 
n.time = N.time(time_20hz_idx);
w.time = W.PLD2_1hz.time(time_1hz_idx);

% Create the structures consisting of data from specified time frame

%%%%%%% 20 Hz data %%%%%%%

% Loop through fieldnames
for k = 1:length(field_20hz)

    % Set filenames
    name = char(field_20hz(k,:));

    % Index data
    eval(['n.' name ' = N.' name '(time_20hz_idx);'])
    
    % Compute the time elapsed from the first measurement in the 30 min
    % interval for each time recorded in seconds
    n.t =(n.time - n.time(1))*86400;
    time_int_n = 0:dt_n:round(max(n.t)); % Elasped time since first measurement. Units: seconds

    % Find indicies of non-NaN values in kth variable
    eval(['Inum = find(~isnan(n.' name '));'])

    % Linearly interpolate time steps with nans
    eval(['n.' name ' = interp1(n.t(Inum),n.' name '(Inum),time_int_n,''linear'',''extrap'');'])

end

%%%%%%% 1 Hz data %%%%%%%

% Loop through fieldnames
for k = 1:length(field_1hz)

    % Set filenames
    name = char(field_1hz(k,:));

    % Index data
    eval(['w.' name ' = W.PLD2_1hz.' name '(time_1hz_idx);'])
    
    % Compute the time elapsed from the first measurement in the 30 min
    % interval for each time recorded in seconds
    w.t =(w.time - w.time(1))*86400;
    time_int_w = 0:dt_w:round(max(w.t)); % Elasped time since first measurement. Units: seconds

    % Find indicies of non-NaN values in kth variable
    eval(['Inum = find(~isnan(w.' name '));'])

    % Linearly interpolate time steps with nans
    eval(['w.' name ' = interp1(w.t(Inum),w.' name '(Inum),time_int_w,''linear'',''extrap'');'])

end

%% Compute environmental conditions paramters for 20 minute time intervals

% Convert the 20 Hz time record into two time records with increments
% specified by period beginning at date_i and ending at date_f
dt_defined = (period + 2)*60;                                               % Units: seconds
T0_i= 0:dt_defined:max(time_int_n);                                         % Units: seconds

% Define the 1 minute time offset for the high pass filtered data
time_offset = 60;                                                           % Units: seconds
np_n = fe_n*time_offset;                                                    % Number of measurements taken in time offset duration
np_w = fe_w*time_offset;                                                    % Number of measurements taken in time offset duration

% Extend T0 to account for the offset from the truncation
i = length(T0_i) - 1 + ((5/2)*(dt_defined/time_offset) + 1);  
T0 = [T0_i(1:end-1) T0_i(end):dt_defined:T0_i(end)+(i+(i-2))*time_offset];

% Initialize counter
kr = 0;

% Loop through time increments
for is = 1:(length(T0) - 1)
    
    % Clear fast structure
    clear fast
    
    % Create index for grabbing data from time interval (data within 1 min of the edges will be removed due to effects of the filter at boundaries)
    if is == 1 
        It_n = find(time_int_n >= T0(is) & time_int_n < T0(is+1));
        It_w = find(time_int_w >= T0(is) & time_int_w < T0(is+1));
    elseif is > 1 
        It_n = find(time_int_n >= (T0(is) - (is+(is-2))*time_offset) & time_int_n < (T0(is+1) - (is+(is-2))*time_offset)); 
        It_w = find(time_int_w >= (T0(is) - (is+(is-2))*time_offset) & time_int_w < (T0(is+1) - (is+(is-2))*time_offset)); 
    end
    
    % Create the fast structure consisting of data from the time increment
    
    %%%%%%% 20 Hz data %%%%%%%
    % Loop through fieldnames
    for k = 1:length(field_20hz)

        % Set filenames
        name = char(field_20hz(k,:));

        % Index data
        %------- Planck -------%
        eval(['fast.' name ' = n.' name '(:,It_n);'])

    end

    %%%%%%% 1 Hz data %%%%%%%
    % Loop through fieldnames
    for k = 1:length(field_1hz)

        % Set filenames
        name = char(field_1hz(k,:));

        % Index data
        eval(['fast.' name ' = w.' name '(:,It_w);'])

    end
    
    % Computing platform speed and direction 
    [fast.ground_speed_ms, fast.true_course] = platform_speed_dir(fast.ve, fast.vn, dir_con); 

    % Index non-datenum time series 
    ttime_n = time_int_n(It_n);
    
    % Set transfer function coefficients for high pass butterworth filter 
    [Bcut,Acut] = butter(4,df/fn,'high');
    
    %%%%%%% 20 Hz data %%%%%%%
    
    % Loop through field names (except lon, lat, and time) 
    for k = 1:(length(field_20hz) - 3)
        
        % Set field variable name
        name = char(field_20hz(k,:));
        
        % Detrend and high pass filter data
        eval(['hp = filtfilt(Bcut,Acut,detrend(double(fast.' name ')));'])
        
        % Truncate high pass time series by one minute on either side to avoid
        % boundary effects from filter
        eval(['fast.high_pass.' name ' = hp(np_n:end-np_n);'])
        
    end 
    
    % Truncate time series ground_speed_ms, and true_course
    fast.high_pass.ground_speed_ms = fast.ground_speed_ms(np_n:end-np_n);
    fast.high_pass.true_course = fast.true_course(np_n:end-np_n);
    fast.high_pass.time = fast.time(np_n:end-np_n);
    fast.high_pass.ttime_n = ttime_n(np_n:end-np_n);
    
    %%%%%%% 1 Hz data %%%%%%%
    
    % Loop through field names  
    for k = 1:length(field_1hz)
        
        % Set field variable name
        name = char(field_1hz(k,:));
        
        % Truncate time series by one minute on either side to be
        % consistent with high pass data
        eval(['fast.high_pass.' name ' = fast. ' name '(np_w:end-np_w);'])
        
    end
    
    % Set counter
    kr = kr + 1;
    
    %Compute Directional Spectrum 
    [Sd, f, ~] = compute_directional_spectrum(fast.high_pass.alt, fast.high_pass.ve, fast.high_pass.vn, fast.high_pass.vu, fast.high_pass.ttime_n, f, df, dtheta, nfft, fe_n, toolbox, variables, scaling, dir_con); 

    % Compute the Omni-directional Spectra from directional wave spectra 
    omni_spectrum = sum(Sd * dtheta, 1);
    
    % Compute the half way point of the time series
    T_leg = n.time(It_n);
    T_half = T_leg(round(length(It_n) - (1/2)*length(It_n)));
    
    % Floor time step 
    T_half_n = floor(T_half*1440)/1440;
    
    % Create EC structure:  
    EC.time(1,kr) = T_half_n; 
    if kr == 1
        EC.f = f;  % frequency (same for all time intervals)
    end
    
    % Compute mean platform speed and direction
    [EC.mD(:,kr), EC.std_D(:,kr)] = direction_stats(fast.high_pass.true_course, dt_n, 0);                      
    EC.mspeed(:,kr) = mean(fast.high_pass.ground_speed_ms); EC.std_speed(:,kr) = std(fast.high_pass.ground_speed_ms); 

    % Compute mean wind speed and direction
    [EC.mTWD(:,kr), EC.std_TWD(:,kr)] = direction_stats(fast.high_pass.TWD, dt_w, false);
    fast.high_pass.mTWS = mean(fast.high_pass.TWS); EC.std_TWS(:,kr) = std(fast.high_pass.TWS);

    %Compute the wind speed at 10 meters above the ocean surface
    [~,~,EC.mTWS_10(:,kr)] = FUNC_ustar_Z(z_station,fast.high_pass.mTWS);
    
    % Compute Significant Wave height: 
    EC.Hs_spec(:,kr) = 4 * sqrt(trapz(EC.f(freq_band), omni_spectrum(freq_band))); % Hs (Spectral method)

    %Display 20 min interval 
    disp(datestr(T_half_n))
    
end

%% Plot platform trajectory and environmental conditions
clc, close all

% Set plotting parameters
cb_l = -4000; cb_h = -1; cb_land = 0;
t_ticks = datetime('29-Oct-2021 00:00:00'):days(1):datetime('04-Nov-2021 00:00:00');
red = [0.6350 0.0780 0.1840]; 
blue = [0 0.4470 0.7410]; 
fontsize = 11;

% Create Figure and axes
fig = figure('units','normalized','outerposition',[0 0 0.45 1], 'Name', 'Experiment site/assets and environmental conditions for SMODE2021');

%--------------- Subplot 1 ---------------%
subplot(5,1,[1 2])

% Plot bathymetry 
pcolor(lon_exp,lat_exp,z_exp');

hold on 

% Plot WHOI43 Wave Glider trajectory
pc2 = plot(N.longitude_n, N.latitude_n, '.', 'color', blue, 'MarkerSize',3);

% Plot markers for beginning and end of trajectory
pc3 = scatter(N.longitude_n(1),  N.latitude_n(1), 30, 'w','o', 'filled', 'MarkerEdgeColor', 'k');
pc4 = scatter(N.longitude_n(end),  N.latitude_n(end), 30, 'w','d', 'filled', 'MarkerEdgeColor', 'k');

hold off

% Set figure attributes
title('(a)')
axis equal
xlabel('Longitude')
ylabel('Latitude')
xlim([lon_exp_l, lon_exp_h])
ylim([lat_exp_l, lat_exp_h])
xticks(-124.2:0.2:-123.4)
yticks(36.8:0.2:37.4)
xl = xticks; yl = yticks;
xticklabels({[num2str(xl(1)) '$^\circ$']; [num2str(xl(2)) '$^\circ$']; [num2str(xl(3)) '$^\circ$']; [num2str(xl(4)) '$^\circ$']; [num2str(xl(5)) '$^\circ$']})
yticklabels({[num2str(yl(1)) '$^\circ$']; [num2str(yl(2)) '$^\circ$']; [num2str(yl(3)) '$^\circ$']; [num2str(yl(4)) '$^\circ$']})
legend([pc2, pc3, pc4], 'WHOI43', 'Initial Position', 'Final Position', 'Location', 'northeast', 'Fontsize', fontsize-3)
grid on
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out','TickLength', [0.015,0.75]);
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Color', 'w')
shading flat

%------- set the colormap of the zoomed in region -------%
% Set the bathymetric colormap
cmap = colormap([flipud(cmocean('deep')); [0.5, 0.5, 0.5]]);                % flipud(cmocean('deep'))

% Create a depth vector corresponding to the levels of the colormap
cmap_z = [linspace(cb_l,cb_h, length(cmap)-1), cb_land]; 

% Find the max and min values of the bathymetry in the experiment site
z_max = max(z_exp,[],'All'); z_min = min(z_exp,[],'All'); 

% Index the colormap depth vector 
idx_bathy = find(cmap_z >= z_min & cmap_z <= z_max);
cmap_z_exp = cmap_z(idx_bathy);

% Initialize new colormap 
cmap_exp = cmap(idx_bathy,:);

% Set colormap
colormap(cmap_exp)

%Display Colorbar
cb = colorbar; 
caxis([cb_l, -1000]);
cb.Label.Interpreter = 'Latex';
cb.FontSize = fontsize;
cb.Label.String = 'Depth (m)';
cb.TickDirection = 'out'; 
cb.TickLabelInterpreter = 'latex';

% Find the position of the current axes
p_ax1 = get(gca, 'Position');                                               % Note: position vector: [left bottom width height]

% Shift axis up and get new position
y_shift = 0.02;
set(gca, 'Position', [p_ax1(1) p_ax1(2)+y_shift p_ax1(3) p_ax1(4)])         % x-offset to align axes labels and titles
p_ax1 = get(gca, 'Position');

% Create a new position vector for the inset axis using that position data 
% from the larger axis.                                             % Note: position vector: [left bottom width height]
p_ax2 = [p_ax1(1)+0.15 p_ax1(2)+0.03 p_ax1(3)-.5 p_ax1(4)-.2];

% Apply empirical corrections to x and y end points for annotate function in
% order to account for the modulation due to the axis square command
x_i = p_ax2(1) + 0.008;
x_f = p_ax2(1) + p_ax2(3) - 0.064;
y_i = p_ax2(2);
y_f = p_ax2(2) + p_ax2(4);

% Set the position vector of the h axes 
x = linspace(x_i, x_f, length(lon_reg));
y = linspace(y_i, y_f, length(lat_reg));

% Apply the experimental site longitude and latitude indices  
Ilon_exp_g = lon_reg >=  lon_exp_l & lon_reg <= lon_exp_h; Ilat_exp_g  = lat_reg >=  lat_exp_l & lat_reg <= lat_exp_h;
x_exp = x(Ilon_exp_g); 
y_exp = y(Ilat_exp_g);

% Create the new position vector for the annotation box
p_box = [x_exp(1),y_exp(1),x_exp(end)- x_exp(1),y_exp(end)- y_exp(1)]; 

% Create a box for the experimental area
annotation('rectangle',p_box,'Color','w', 'LineWidth', 1.5)

% Create new axis 
h = axes('Parent', gcf, 'Position', p_ax2);

% Plot bathymetry 
pcolor(lon_reg,lat_reg,z_reg')

hold on 

% Plot coastline
contour(lon_reg,lat_reg,z_reg',[0,0], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '-')

% Set figure attributes
title('(b)', 'color', 'w')
axis square
xl = xticks; yl = yticks;
xticklabels({[num2str(xl(1)) '$^\circ$']; [num2str(xl(2)) '$^\circ$']})
yticklabels({[num2str(yl(1)) '$^\circ$']; [num2str(yl(2)) '$^\circ$']; [num2str(yl(3)) '$^\circ$']; [num2str(yl(4)) '$^\circ$']})
set(gca, 'FontSize', 8)
set(gca,'TickDir','out','TickLength', [0.025,0.75]);
set(gca,'XColor','w')
set(gca,'YColor','w')
set(gcf, 'InvertHardcopy', 'off')
set(gca, 'color', 'w')
set(gca,'TickLabelInterpreter','latex')
shading flat
cmap = colormap(h, cmap);

%Display Colorbar
cb = colorbar; 
caxis([cb_l, cb_land]);
cb.Label.Interpreter = 'Latex';
cb.FontSize = 8;
cb.Label.String = 'Depth (m)';
cb.TickDirection = 'out';
cb.Color = 'w';
cb.Ticks = linspace(cb_l,cb_land,5);
cb.TickLabels = num2cell(linspace(cb_l,cb_land,5)); 
cb.TickLabelInterpreter = 'latex';

% Set the distance for 
L = 5;                                                                      % Edge length of small box track (units: kilometers)
phi = 33;                                                                   % Rounded Latitude of experiment site 
lon_L = (1/(111.32*cosd(phi)))*L;                                           % Degrees of longitude for 500 meters at latitude phi

% Reset the position vector to account for the modulation due to the axis square command and
% colorbar
p_ax3 = [p_ax1(1)+0.115 p_ax1(2) p_ax1(3)-0.23 p_ax1(4)];                 % Note: position vector: [left bottom width height]

% Set the x and y for ax3
x_ax1 = linspace(p_ax3(1), p_ax3(1) + p_ax3(3), length(lon_exp));
y_ax1 = linspace(p_ax3(2), p_ax3(2) + p_ax3(4), length(lat_exp));

% Find the end point in degree longitude of the 500 m distance from the
% right left corner of the map
lon_f = lon_exp(1) + lon_L; 

% Find the index of the closest longitude grid point to the lon_f value
[~,idx_lon]=min(abs(lon_exp - lon_f));

% Find the positions on the ax1 figure axis for the length scale 
px_i = x_ax1(1); 
px_f = x_ax1(idx_lon); 
py_i = y_ax1(1)+0.281; 
py_f = y_ax1(1)+0.281; 

% Create scale bar
% annotation('textbox',...
%     [px_i py_i px_f - px_i  0.0161290322580648],...
%     'String',{'5$\:$km'},'Interpreter','latex', ...
%     'HorizontalAlignment','center','FontSize',7,'FitBoxToText','off',...
%     'FaceAlpha',0.7,'EdgeColor','none','BackgroundColor',[0.8 0.8 0.8]);

%----- Left Vertical Bar -----%
annotation('line',[0.6 0.6],...
    [0.828+0.02 0.838+0.02]+y_shift,'Color',[1 1 1],'LineWidth',1);

%----- Right Vertical Bar -----%
annotation('line',[0.6+(px_f-px_i) 0.6+(px_f-px_i)],...
    [0.828+0.02 0.839+0.02]+y_shift,'Color',[1 1 1],'LineWidth',1);

%----- Horizontal Bar -----%
annotation('line',[0.6  0.6+(px_f-px_i)],...
    [0.833+0.02 0.833+0.02]+y_shift,'Color',[1 1 1],'LineWidth',1);

%----- Text Box -----%
annotation('textbox',...
    [0.627543209876543 0.829301075268817+0.02+y_shift 0.0591851851851851 0.013440860215054],...
    'Color',[1 1 1],'String',{'5$\:$km'},'Interpreter','latex',...
    'FontSize',8,'FitBoxToText','off','EdgeColor','none');

%--------------- Subplot 2 ---------------%
subplot(5,1,3)

% Plot the Wave Glider headin2
plot(EC.time,EC.mD, '^-', 'Color', blue, 'MarkerSize', 4, 'MarkerFaceColor',blue, 'LineWidth', 0.5);

% Set axis attributes
yticks([0,90,180,270,360])
ylim([0,360])
ylabel({'Wave Glider'; 'Heading ($^{\circ}$)'}, 'Interpreter', 'latex')
set(gca,'YColor','k')

% Set figure attributes
title('(c)')
xlim([EC.time(1), EC.time(end)])
xticks(datenum(t_ticks))
datetick('x', 'mmm dd', 'keepticks')
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
set(gca, 'YGrid', 'off', 'XGrid', 'on') 

%--------------- Subplot 3 ---------------%
subplot(5,1,4)

% Plot the Wind direction
yyaxis left
pc1 = plot(EC.time(1:1:end),EC.mTWD(1:1:end), '^', 'Color', blue, 'MarkerSize', 4, 'MarkerFaceColor',blue);

% Set axis attributes
yticks([0,90,180,270,360])
ylim([0,360])
ylabel('Wind Direction ($^{\circ}$)', 'Interpreter', 'latex')
set(gca,'YColor','k')
set(gca,'TickLabelInterpreter','latex')

% Plot the Wind speeds
yyaxis right
pc2 = plot(EC.time,EC.mTWS_10, '-', 'Color', 'k', 'LineWidth', 1.5);

% Set axis attributes
yticks(0:2:12)
ylim([0 12])
ylabel('U$_{10}$ (ms$^{-1}$)', 'Interpreter', 'latex')
set(gca,'YColor','k')

% Set figure attributes
title('(d)')
xlim([EC.time(1), EC.time(end)])
xticks(datenum(t_ticks))
datetick('x', 'mmm dd', 'keepticks')
legend([pc1, pc2], 'Direction', 'Speed', 'Interpreter', 'Latex', 'Orientation', 'vertical', 'Position',[0.140500922813826 0.291315668179691 0.113203113461718 0.0318548390942235], 'Fontsize', fontsize-3);
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
set(gca, 'YGrid', 'off', 'XGrid', 'on') 

%------------- Subplot 4 -------------%
subplot(5,1,5);

% plot significant wave height 
pc1 = plot(EC.time, EC.Hs_spec, '-', 'Color', blue, 'LineWidth', 1.5);
 
% Set figure attributes
title('(e)')
xlabel('UTC time from 29 Oct 2021', 'Interpreter', 'latex')
ylabel('H$_s$ (m)')
xlim([EC.time(1), EC.time(end)])
ylim([1, 5])
xticks(datenum(t_ticks))
datetick('x', 'mmm dd', 'keepticks')
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
set(gca, 'YGrid', 'off', 'XGrid', 'on') 

% Save Figure
saveas(fig, [fig_path 'fig02.png'])

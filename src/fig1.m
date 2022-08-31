%% Figure 1: Overview of experiment site/assets and environmental conditions for DELMAR2020
% Luke Colosi | lcolosi@ucsd.edu | October 5th, 2021

%--------------- Notes ---------------%
% 1. Cutoff wavelength corresponding to approximately half the length of 
%    the SV3 wave glider.   
%-------------------------------------%

% Caption: (a) Trajectories of Planck (blue) and Stokes (red) Wave 
% Gliders during the Delmar2020 experiment for the time period of 
% September 9th at 2:30:00 UTC to September 11th at 16:10:00 UTC. (b)
% Geographic region with location of experiment site (white box). Mean 
% platform direction (c), Wind speed and direction (d), and significant
% wave height (e) display the environmental conditions and additional
% pertinent trajectory information.  

clc, clearvars -except Np Ns Wp Ws bathy, close all;

% Set the current working directory
cd ~/Desktop/projects/asi_lab_summer_internship/WaveSpectrum/src

% Set text interpreter 
set(0,'defaultTextInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')

% Set path to data
%%%% Vehicle %%%%.
vehicle1 = 'PLANCK';
vehicle2 = 'STOKES';

%%%% Root %%%%
ROOT = '../data/DELMAR2020/';
ROOT_bathy = '../data/sccoos/';

% Set path for figures
fig_path = '../figs/';

%%%%%%%%%% Initial global variables %%%%%%%%%%

% Physical parameters
z_station = 1;                                                              % Height of weather station (meters)
g = 9.81;                                                                   % gravitational acceleration (m/s) 
dir_con = {'CW', 'cf', 'rn'};                                               % Directional conventions for wave glider trajectory (wind direction is CW, coming from, reference north)

% Temporal parameters
period = 20;                                                                % Time duration for computing wave spectra
date_i = '09-Sep-2020 02:30:00';                                            % Start date 
date_f = '11-Sep-2020 16:10:00';                                            % End date;
%------- Novatel ------- %
fe_n = 20;                                                                  % Sampling rate (Hz)
dt_n = 1/fe_n;                                                              % Time interval between measurements (s)
%------- Weather Station ------- %
fe_w = 1; 
dt_w = 1/fe_w; 

% Spectral parameters                                                      
fn = (1/2)*fe_n;                                                            % Nyquist Frequency (Hz)
df = 0.01;                                                                  % Frequency resolution
nfft = fe_n/df;                                                             % Window Length
f = (0:df:fn);                                                              % Observed frequency 
lambda_c = 1.5613;                                                          % Wavelength cutoff (meters)
f_noise = sqrt(g/(2*pi*lambda_c));                                          % noise frequency cutoff 
freq_band = find(f >= 0.02 & f <= f_noise);                                 %Frequency band (for computing significant wave height)

%% Call Data

%--------- PLANCK ---------% 
if isempty(whos('Np')) || isempty(whos('Wp'))
    Np = load([ROOT vehicle1 '/NOVATEL_downsampled_20Hz_ALL_' vehicle1 '.mat']); % Downsampled GPS/IMU Data (20 Hz)
    Wp = load([ROOT vehicle1 '/PLD2_1hz_ALL.mat']); % WXT/GILL Weather Station Data (1 hz)
    
    % Transpose time_20hz field: 
    Np.nov.time_20hz = Np.nov.time_20hz';

end 

%--------- STOKES ---------% 
if isempty(whos('Ns')) || isempty(whos('Ws'))
    Ns = load([ROOT vehicle2 '/NOVATEL_downsampled_20Hz_ALL_' vehicle2 '.mat']); % RAW downsampled GPS/IMU Data (20 Hz)
    Ws = load([ROOT vehicle2 '/PLD2_1hz_ALL.mat']); % WXT/GILL Weather Station Data (1 hz)
    
    % Transpose time_20hz field: 
    Ns.nov.time_20hz = Ns.nov.time_20hz';

end

%--------- Bathymetry ---------%
if isempty(whos('bathy'))
    bathy = importdata([ROOT_bathy '/sdNorth.xyz']);
end

% Get coordinate variables
lon = unique(bathy(:,1)); 
lat = unique(bathy(:,2)); 

% Obtain dimensions of the data
nx = length(lon) ; 
ny = length(lat) ;

% Frame matrix in grid 
z = flip(reshape(bathy(:,3),[nx,ny]), 2);

%% Index Wave Glider trajectory from time period and bathymetry for experiment site   

% Find temporal indices for the experimental time period
time_idxp = find(Np.nov.time_20hz >= datenum(date_i) & Np.nov.time_20hz < datenum(date_f));
time_idxs = find(Ns.nov.time_20hz >= datenum(date_i) & Ns.nov.time_20hz < datenum(date_f));

% Apply to Wave glider latitude and longitude trajectory
Np.nov.latitude_n = Np.nov.latitude(time_idxp);
Np.nov.longitude_n = Np.nov.longitude(time_idxp);
Ns.nov.latitude_n = Ns.nov.latitude(time_idxs);
Ns.nov.longitude_n = Ns.nov.longitude(time_idxs);

% Initialize spatial parameters
lon_exp_l = -117.70 - 0.0167; lon_exp_h = -117.56 + 0.0167;
lat_exp_l = 32.91 - 0.0167; lat_exp_h = 33.01 + 0.0167;

% Find spatial indices for the experimental site
lon_idx_exp = find(lon >= lon_exp_l & lon <= lon_exp_h);
lat_idx_exp = find(lat >= lat_exp_l & lat <= lat_exp_h);

% Apply to bathymetry and lonitude data 
lon_exp = lon(lon_idx_exp); lat_exp = lat(lat_idx_exp); z_exp = z(lon_idx_exp,lat_idx_exp); 

%% Compute environmental conditions paramters for 20 minute time intervals

% Set variables
field_20hz = {'heave';'VEL_east'; 'VEL_north'; 'time_20hz'};
field_1hz = {'TWD'; 'TWS'; 'time'};

% Convert the 20 Hz time record into two time records with increments
% specified by period
T0 = unique( floor(Np.nov.time_20hz*24*60/period)/(24*60/period) );
T0(T0 < datenum(date_i)) = []; 
T0(T0 > datenum(date_f)) = [];
T1 = T0 + 1/(24*60/period);

% Initialize counter
kr = 0;

% Loop through time increments
for is = 1:length(T0)
    
    % Clear nov and w structures
    clear nov_p nov_s w_p w_s
    
    % Create index for grabbing data from ith time interval
    It_p = find(Np.nov.time_20hz >= T0(is) & Np.nov.time_20hz < T1(is));
    It_s = find(Ns.nov.time_20hz >= T0(is) & Ns.nov.time_20hz < T1(is));
    It_wp = find(Wp.PLD2_1hz.time >= T0(is) & Wp.PLD2_1hz.time < T1(is));
    It_ws = find(Ws.PLD2_1hz.time >= T0(is) & Ws.PLD2_1hz.time < T1(is));
    
    %Display ith time interval 
    disp(datestr(T0(is)))
    
    % Create the nov and w structures consisting of data from the ith time
    % interval
    
    %%%%%%% 20 Hz data %%%%%%%
    % Loop through fieldnames
    for k = 1:length(field_20hz)

        % Set filenames
        name = char(field_20hz(k,:));

        % Index data
        %------- Planck -------%
        eval(['nov_p.' name ' = Np.nov.' name '(It_p,:);'])
        %------- Stokes -------%
        eval(['nov_s.' name ' = Ns.nov.' name '(It_s,:);'])

    end

    %%%%%%% 1 Hz data %%%%%%%
    % Loop through fieldnames
    for k = 1:length(field_1hz)

        % Set filenames
        name = char(field_1hz(k,:));

        % Index data
        %------- Planck -------%
        eval(['w_p.' name ' = Wp.PLD2_1hz.' name '(:,It_wp);'])
        %------- Stokes -------%
        eval(['w_s.' name ' = Ws.PLD2_1hz.' name '(:,It_ws);'])

    end
    
    % Compute the time elapsed from the first measurement in the time
    % interval for each time recorded in seconds
    nov_p.t = (nov_p.time_20hz - nov_p.time_20hz(1))*86400; time_int_np = 0:dt_n:round(max(nov_p.t,[],'omitnan'));
    nov_s.t = (nov_s.time_20hz - nov_s.time_20hz(1))*86400; time_int_ns = 0:dt_n:round(max(nov_s.t,[],'omitnan'));
    w_p.t = (w_p.time - w_p.time(1))*86400; time_int_wp = 0:dt_w:round(max(w_p.t,[],'omitnan'));
    w_s.t = (w_s.time - w_s.time(1))*86400; time_int_ws = 0:dt_w:round(max(w_s.t,[],'omitnan'));

    % Check if more than 6000 ground speed measurements (10 mins) exists to
    % continue
    Inum_p = find(~isnan(nov_p.VEL_east)); Inum_s = find(~isnan(nov_s.VEL_east)); 
    if length(Inum_p) < 6000 && length(Inum_s) < 6000
        continue
    end
    
    %%%%%%% 20 Hz data %%%%%%%
    
    % Loop through field names
    for k = 1:length(field_20hz)
        
        % Set field variable name
        name = char(field_20hz(k,:));
        
        %------- Planck -------%
        % Find indices of non-nan data points
        eval(['Inum = find(~isnan(nov_p.' name '));'])
        
        % Interpolate time steps with nans
        eval(['nov_p.' name ' = interp1(nov_p.t(Inum),nov_p.' name '(Inum),time_int_np,''linear'',''extrap'');'])
        
        %------- Stokes -------%
        % Find indices of non-nan data points
        eval(['Inum = find(~isnan(nov_s.' name '));'])
        
        % Interpolate time steps with nans
        eval(['nov_s.' name ' = interp1(nov_s.t(Inum),nov_s.' name '(Inum),time_int_ns,''linear'',''extrap'');'])
        
    end 
    
    %%%%%%% 1 Hz data %%%%%%%
    
    % Loop through field names (except time) 
    for k = 1:length(field_1hz)
        
        % Set field variable name
        name = char(field_1hz(k,:));
        
        %------- Planck -------%
        % Find indices of non-nan data points
        eval(['Inum = find(~isnan(w_p.' name '));'])
        
        % Interpolate time steps with nans
        eval(['w_p.' name ' = interp1(w_p.t(Inum),w_p.' name '(Inum),time_int_wp,''linear'',''extrap'');'])
        
        %------- Stokes -------%
        % Find indices of non-nan data points
        eval(['Inum = find(~isnan(w_s.' name '));'])
        
        % Interpolate time steps with nans
        eval(['w_s.' name ' = interp1(w_s.t(Inum),w_s.' name '(Inum),time_int_ws,''linear'',''extrap'');'])
        
    end
    
    % Set counter
    kr = kr + 1;

    %--- Computing platform speed and direction ---% 
    [nov_p.ground_speed_ms, nov_p.true_course] = platform_speed_dir(nov_p.VEL_east, nov_p.VEL_north, dir_con);
    [nov_s.ground_speed_ms, nov_s.true_course] = platform_speed_dir(nov_s.VEL_east, nov_s.VEL_north, dir_con);

    % Compute the power density spectrum for the vertical position with 50% overlap:
    [Szz_p,~]  = pwelch(nov_p.heave',nfft,nfft/2,nfft,fe_n);
    [Szz_s,~]  = pwelch(nov_s.heave',nfft,nfft/2,nfft,fe_n);
    
    % Create EC structure:  
    EC.time(1,kr) = T0(is);
    if kr == 1
        EC.f = f;
    end

    % Compute mean platform speed and direction
    [EC.mD_p(:,kr), EC.std_D_p(:,kr)] = direction_stats(nov_p.true_course, dt_n, 0);                      
    EC.mspeed_p(:,kr) = mean(nov_p.ground_speed_ms); EC.std_speed_p(:,kr) = std(nov_p.ground_speed_ms); 
    [EC.mD_s(:,kr), EC.std_D_s(:,kr)] = direction_stats(nov_s.true_course, dt_n, 0);                      
    EC.mspeed_s(:,kr) = mean(nov_s.ground_speed_ms); EC.std_speed_s(:,kr) = std(nov_s.ground_speed_ms);  
    
    % Compute mean wind speed and direction
    [EC.mTWD_p(:,kr), EC.std_TWD_p(:,kr), ~] = direction_stats(w_p.TWD, dt_w, 0); % Convention: CW, coming From, ref north
    [EC.mTWD_s(:,kr), EC.std_TWD_s(:,kr), ~] = direction_stats(w_s.TWD, dt_w, 0); % Convention: CW, coming From, ref north
    w_p.mTWS = mean(w_p.TWS); EC.std_TWS_p(:,kr) = std(w_p.TWS);
    w_s.mTWS = mean(w_s.TWS); EC.std_TWS_p(:,kr) = std(w_s.TWS);

    %Compute the wind speed at 10 meters above the ocean surface
    [~,~,EC.mTWS_10_p(:,kr)] = FUNC_ustar_Z(z_station,w_p.mTWS);
    [~,~,EC.mTWS_10_s(:,kr)] = FUNC_ustar_Z(z_station,w_s.mTWS);
    
    % Compute significant Wave height: 
    EC.Hs_spec_p(:,kr) = 4 * sqrt(trapz(EC.f(freq_band),Szz_p(freq_band))); EC.Hs_temp1_p(:,kr) = 4 * std(nov_p.heave,'omitnan'); EC.Hs_temp2_p(:,kr) = 4 * mean(nov_p.heave.^2,'omitnan');
    EC.Hs_spec_s(:,kr) = 4 * sqrt(trapz(EC.f(freq_band),Szz_s(freq_band))); EC.Hs_temp1_s(:,kr) = 4 * std(nov_p.heave,'omitnan'); EC.Hs_temp2_s(:,kr) = 4 * mean(nov_p.heave.^2,'omitnan');
    
end 

%% Plot platform trajectory and environmental conditions
clc, close all

% Set plotting parameters
cb_l = -1000; cb_h = 1000;
t_ticks = datetime('09-Sep-2020 00:00:00'):hours(12):datetime('12-Sep-2020 00:00:00');
red = [0.6350 0.0780 0.1840]; 
blue = [0 0.4470 0.7410]; 
fontsize = 11;

% Create Figure and axes
fig = figure('units','normalized','outerposition',[0 0 0.45 1]); % [0 0 1 1]

%--------------- Subplot 1 ---------------%
subplot(5,1,[1 2]) % subplot(3,2,[1 3 5])

% Plot bathymetry 
pcolor(lon_exp,lat_exp,z_exp');

hold on 

% Plot PLANCK Wave Glider trajectory
pc2 = plot(Np.nov.longitude_n, Np.nov.latitude_n, '.', 'color', blue, 'MarkerSize',3);

% Plot STOKES Wave Glider trajectory
pc3 = plot(Ns.nov.longitude_n, Ns.nov.latitude_n, '.', 'color', red, 'MarkerSize',3);

% Plot markers for beginning and end of trajectory
pc4 = scatter(Np.nov.longitude_n(1),  Np.nov.latitude_n(1), 30, 'w','o', 'filled', 'MarkerEdgeColor', 'k');
pc5 = scatter(Np.nov.longitude_n(end),  Np.nov.latitude_n(end), 30, 'w','d', 'filled', 'MarkerEdgeColor', 'k');
scatter(Ns.nov.longitude_n(1), Ns.nov.latitude_n(1), 30, 'w', 'o', 'filled', 'MarkerEdgeColor', 'k');
scatter(Ns.nov.longitude_n(end),  Ns.nov.latitude_n(end), 30, 'w','d', 'filled', 'MarkerEdgeColor', 'k');

hold off

% Set figure attributes
title('(a)')
axis equal
xlim([lon_exp_l, lon_exp_h])
ylim([lat_exp_l, lat_exp_h])
%xlabel('Longitude', 'Interpreter', 'latex')
%ylabel('Latitude', 'Interpreter', 'latex')
legend([pc2, pc3, pc4, pc5], 'Planck', 'Stokes', 'Initial Position', 'Final Position', 'Location', 'northeast', 'Fontsize', fontsize)
grid on
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out','TickLength', [0.015,0.75]);
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Color', 'w')
shading flat

%------- set the colormap of the zoomed in region -------%
% Set the topographic colormap
cmap = colormap(cmocean('topo'));

% Create a depth vector corresponding to the levels of the colormap
cmap_z = linspace(cb_l,cb_h, 256); 

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
caxis([cb_l, -700]);
cb.Label.Interpreter = 'Latex';
cb.FontSize = fontsize;
cb.Label.String = 'Elevation (m)';
cb.TickDirection = 'out';
cb.Ticks = linspace(cb_l,-700,6);
cb.TickLabels = num2cell(linspace(cb_l,-700,6)); 
cb.TickLabelInterpreter = 'latex';

% Find the position of the current axes and create a new position vector 
% for the new axis using that position data.
p_ax1 = get(gca, 'Position');                                               % Note: position vector: [left bottom width height]
p_ax2 = [p_ax1(1)+0.165 p_ax1(2)+0.035 p_ax1(3)-.5 p_ax1(4)-.2];           


% Apply empirical corrections to x and y end points for annotate function in
% order to account for the modulation due to the axis square command and
% colorbar
x_i = p_ax2(1) +0.01;                                                      
x_f = p_ax2(1) + p_ax2(3) - 0.078;                                         
y_i = p_ax2(2) + 0.001;                                                     
y_f = p_ax2(2) + p_ax2(4);                                                  

% Set the position vector of the h axes 
x = linspace(x_i, x_f, length(lon));
y = linspace(y_i, y_f, length(lat));

% Apply the experimental site longitude and latitude indices  
x_exp = x(lon_idx_exp); 
y_exp = y(lat_idx_exp);

% Create the new position vector for the annotation box
p_box = [x_exp(1),y_exp(1),x_exp(end)- x_exp(1),y_exp(end)- y_exp(1)]; 

% Create a box for the experimental area
annotation('rectangle',p_box,'Color','w', 'LineWidth', 1.5)

% Create new axis 
h = axes('Parent', gcf, 'Position', p_ax2);

% Plot bathymetry 
pcolor(lon,lat,z')

hold on 

% Plot coastline
contour(lon,lat,z',[0,0], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '-')

% Set figure attributes
title('(b)', 'color', 'w')
axis square
xticks(-117.7:0.2:-117.1)
set(gca, 'FontSize', 8)
set(gca,'TickDir','out','TickLength', [0.025,0.75]);
set(gca,'XColor','w')
set(gca,'YColor','w')
set(gcf, 'InvertHardcopy', 'off')
set(gca, 'color', 'w')
set(gca,'TickLabelInterpreter','latex')
shading flat
cmap = colormap(h, cmocean('topo'));

%Display Colorbar
cb = colorbar; 
caxis([cb_l, cb_h]);
cb.Label.Interpreter = 'Latex';
cb.FontSize = 8;
cb.Label.String = 'Elevation (m)';
cb.TickDirection = 'out';
cb.Color = 'w';
cb.Ticks = linspace(cb_l,cb_h,9);
cb.TickLabels = num2cell(linspace(cb_l,cb_h,9)); 
cb.TickLabelInterpreter = 'latex';

% Create winding number arrows annotations
%----- Planck Large Box -----%
xa = [0.373456790123456 0.337962962962963];
ya = [0.874163319946452 0.874163319946452];
annotation('line',xa,ya,'Color','w', 'LineWidth', 1.5)
xa = [0.33891975308642 0.33891975308642];
ya = [0.875163319946452 0.845373493975903];
annotation('textarrow',xa,ya,'String','10 Rev ','Color','w', 'LineWidth', ...
           1.5, 'Interpreter', 'latex', 'FontSize', 8,'HeadWidth',9,...
           'HeadLength',6)

%----- Planck Small Box -----%
xa = [0.632716049382716 0.612654320987653];
ya = [0.661311914323962 0.661311914323961];
annotation('line',xa-0.001,ya,'Color','w', 'LineWidth', 1.5)
xa = [0.611111111111109 0.611111111111109];
ya = [0.660311914323962 0.681392235609102];
annotation('textarrow',xa,ya,'String','9 Rev ','Color','w', 'LineWidth', ...
           1.5, 'Interpreter', 'latex', 'FontSize', 8,'HeadWidth',9,...
           'HeadLength',6)

%----- Stokes Large Box -----%
xa = [0.449074074074072 0.413580246913579];
ya = [0.875502008032128 0.875502008032128];
annotation('line',xa,ya,'Color','w', 'LineWidth', 1.5)
xa = [0.412993827160493 0.412993827160493];
ya = [0.876502008032128 0.846712182061579];
annotation('textarrow',xa,ya,'String','11 Rev ','Color','w', 'LineWidth', ...
            1.5, 'Interpreter', 'latex', 'FontSize', 8,'HeadWidth',9,...
            'HeadLength',6)

%----- Stokes Small Box -----%
xa = [0.666666666666666 0.666666666666667];
ya = [0.659973226238286 0.682730923694779];
annotation('line',xa,ya+0.001,'Color','w', 'LineWidth', 1.5)
xa = [0.666666666666666 0.645061728395059];
ya = [0.661311914323962 0.661311914323961];
annotation('textarrow',xa,ya,'String','9 Rev ', 'TextRotation', 0, 'Color', ...
           'w', 'LineWidth', 1.5, 'Interpreter', 'latex', 'FontSize', 8,...
           'HorizontalAlignment','left', 'VerticalAlignment','top',...
           'HeadWidth',9, 'HeadLength',6)

%--------------- Subplot 2 ---------------%
subplot(5,1,3)

% Plot the Wave Glider heading
plot(EC.time,EC.mD_p, '^-', 'Color', blue, 'MarkerSize', 4, 'MarkerFaceColor',blue, 'LineWidth', 0.5);

hold on 
plot(EC.time,EC.mD_s, '^-', 'Color', red, 'MarkerSize', 4, 'MarkerFaceColor',red, 'LineWidth', 0.5);
hold off 

% Set axis attributes
yticks([0,90,180,270,360])
ylim([0,360])
ylabel('Wave Glider Heading ($^{\circ}$)', 'Interpreter', 'latex')
set(gca,'YColor','k')

% Set figure attributes
title('(c)')
xticks(datenum(t_ticks))
datetick('x', 'mmm dd', 'keepticks')
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlim([EC.time(1), EC.time(end)])
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
set(gca, 'YGrid', 'off', 'XGrid', 'on') 

%--------------- Subplot 3 ---------------%
subplot(5,1,4)

% Plot the Wind direction
yyaxis left
plot(EC.time(1:1:end),EC.mTWD_p(1:1:end), '^', 'Color', blue, 'MarkerSize', 4, 'MarkerFaceColor',blue)

hold on 
plot(EC.time(1:1:end),EC.mTWD_s(1:1:end), '^', 'Color', red, 'MarkerSize', 4, 'MarkerFaceColor',red)
hold off 

% Set axis attributes
yticks([0,90,180,270,360])
ylim([0,360])
ylabel('Wind Direction ($^{\circ}$)', 'Interpreter', 'latex')
set(gca,'YColor','k')
set(gca,'TickLabelInterpreter','latex')

% Plot the Wind speeds
yyaxis right
plot(EC.time,EC.mTWS_10_p, '-', 'Color', blue, 'LineWidth', 1.5);

hold on 
plot(EC.time,EC.mTWS_10_s, '-', 'Color', red, 'LineWidth', 1.5);
hold off 

% Set axis attributes
yticks(0:2:10)
ylim([0,10])
ylabel('U$_{10}$ (ms$^{-1}$)', 'Interpreter', 'latex')
set(gca,'YColor','k')

% Set figure attributes
title('(d)')
xticks(datenum(t_ticks))
datetick('x', 'mmm dd', 'keepticks')
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlim([EC.time(1), EC.time(end)])
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
set(gca, 'YGrid', 'off', 'XGrid', 'on') 

%------------- Subplot 4 -------------%
subplot(5,1,5);

% plot significant wave height 
pc1 = plot(EC.time, EC.Hs_spec_p, '-', 'Color', blue, 'LineWidth', 1.5);

hold on 
pc2 = plot(EC.time, EC.Hs_spec_s, '-', 'Color', red, 'LineWidth', 1.5);
hold off
 
% Set figure attributes
title('(e)')
ylabel('H$_s$ (m)')
ylim([0.8, 1.4])
xticks(datenum(t_ticks))
datetick('x', 'mmm dd', 'keepticks')
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlim([EC.time(1), EC.time(end)])
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
legend([pc1, pc2], 'Planck', 'Stokes', 'Interpreter', 'Latex', 'Location', 'northeast');
set(gca, 'YGrid', 'off', 'XGrid', 'on') 

% Save Figure
saveas(fig, [fig_path 'figure_1.png'])

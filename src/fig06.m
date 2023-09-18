%% Figure 6: Example of a directional spectrum 
% Luke Colosi | lcolosi@ucsd.edu | August 20th, 2022

%-------------------------------- Caption --------------------------------%
% An observed directional wave spectrum computed from measurements 
% collected onboard Wave Glider Stokes during a single fixed heading leg 
% of the large 1000 m box trajectory in the DELMAR2020 experiment 
% (9 September 2020 at 23:57:06 to 10 September 2020 at 00:25:18).
%-------------------------------------------------------------------------%

clc, clear, close all;

% Set default interpreter to latex for all text in figures
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

% Set path to data
%%%% Vehicle %%%%.
vehicle = 'STOKES';

%%%% Root %%%%
ROOT = '../data/DELMAR2020/';

% Set path for figures
fig_path = '../figs/';

%%%%%%%%%% Initial global variables %%%%%%%%%%

% Physical parameters
g = 9.81;                                                                   % gravitational acceleration (m/s) 
dir_con = {'CW', 'cf', 'rn'};                                               % Directional conventions

% Temporal parameters
period = 10;                                                                % Time duration for computing wave spectra
date_i = '09-Sep-2020 02:30:00';                                            % Start date 
date_f = '11-Sep-2020 16:10:00';                                            % End date
date_o = {date_i, date_f}; 
%------- Novatel ------- %
fe_n = 20;                                                                  % Sampling rate (Hz)
dt_n = 1/fe_n;                                                              % Time interval between measurements (s)
%------- Weather Station ------- %
fe_w = 1; 
dt_w = 1/fe_w; 

% Spectral parameters                                                      
fn = (1/2)*fe_n;                                                            % Nyquist Frequency (Hz)
df = 0.01;                                                                  % Frequency resolution
nfft = fe_n/df;                                                             % Window length 
f = (0:df:fn);                                                              % Observed frequency (Hz)
dtheta = 5;                                                                 % Angular resolution (degrees)
ntheta = 360/dtheta+1;                                                      % Number of angles
lambda_c = 1.5613;                                                          % Wavelength cutoff (meters)
f_noise = sqrt(g/(2*pi*lambda_c));                                          % Noise frequency cutoff 
toolbox = 'WAFO';                                                           % Method used to compute directional spectrum 
variables = 'heave_velocity';                                               % Heave and horizontal velocity are used to compute the direction spectrum
scaling = false;                                                            % Variance of directional spectrum is not scaled to match variance of heave spectrum 
method = 'MEM';                                                             % Method used to compute directional spectrum 

% Upload and process novatel and weather station data
[nov_s, ~, nlegs_s] = process_wg_data(vehicle, ROOT, date_o, dir_con, [dt_n, dt_w]);

%% Compute Wave Spectra

% Loop through legs
for n = 1:nlegs_s
    
    % Compute Directional Spectrum  
    eval(['[nov_s.Sd(:,:,n), nov_s.f, nov_s.theta] = compute_directional_spectrum(nov_s.L' num2str(n) '.heave, nov_s.L' num2str(n) '.VEL_east, nov_s.L' num2str(n) '.VEL_north, nov_s.L' num2str(n) '.VEL_up, nov_s.L' num2str(n) '.time_20hz, f, df, dtheta, nfft, fe_n, toolbox, variables, scaling, dir_con, method);']) 
    
end

% Obtain indicies frequency below high frequency cut off
Inoise_s = find(nov_s.f < f_noise);

% Remove Frequencies beyond 1 Hz from correction
nov_s.f_ob = nov_s.f(Inoise_s);

% Remove spectral densities beyond 1 Hz 
nov_s.Sd_f_ob = nov_s.Sd(:,Inoise_s,:);

%% Plot directional wave spectrum 

% Set plotting parameters
Nspokes = 13;                                                               % Number of radial lines
Ncircles = 6;                                                               % Number of azimuthal circles
pos = linspace(nov_s.f_ob(2), nov_s.f_ob(end), Ncircles);                   % Position of frequency labels
theta_s = linspace(0,360,length(nov_s.theta));                              % Set theta for plotting
Rticks = {'0.01 Hz','0.03 Hz','0.06 Hz','0.16 Hz', '0.40 Hz', '1 Hz'};      % Frequency labels (obtained by logspace(log10(nov_p.f_in(2)), log10(nov_p.f_in(end)), Ncircles))
Contours = logspace(-9,-6.5,30);                                                   
cmap = colormap(flipud(cbrewer2('RdYlBu', numel(Contours)))); 
Contours_n = [10^-16, Contours];
cmap_n = cat(1,cmap(1,:), cmap);
fontsize = 28;
LineWidth = 0.5;
LineColor = 'k';
LineStyle = '-';
TextColor = 'w';
itime = 18; 

% Smooth out directional spectra
Sd_smooth = movmean(nov_s.Sd_f_ob(:,:,itime),3,2);

% Create figure
fig = figure('units','normalized','outerposition',[0 0 1 1]);

% Plot Directional spectra  
[~,cb] = polarcontourf(nov_s.f_ob(2:end), theta_s, Sd_smooth(:,2:end), ...
                       'circlesPos', pos, 'Rscale', 'log',...
                       'Ncircles', Ncircles, 'Nspokes', Nspokes,...
                       'typeRose' ,'meteo', 'Contours', Contours_n,...
                       'colBar', 3, 'fontsize', fontsize, 'colormap', cmap_n,...
                       'LineWidth', LineWidth, 'LineColor', LineColor,...
                       'LineStyle', LineStyle, 'RtickLabel', Rticks, ...
                       'TextColor', TextColor);

% Set logarithmic colorbar
cb.Ticks = log10([10^-9, 10^-8, 10^-7]);
cb.TickLabels = {'10$^{-9}$'; '10$^{-8}$'; '10$^{-7}$'}; 
cb.Label.String = 'S$_{ob}$($f_{ob}$, $\theta$) (m$^2$ Hz$^{-1}$ deg$^{-1}$)';
cb.TickDirection = 'out';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = fontsize;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = fontsize;

% Set other figure attributes 
set(gcf,'color','w')
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')

% Display time frame directional spectrum is computed over
t_initial = eval(['datestr(nov_s.L' num2str(itime) '.time_20hz(1))']); 
t_final = eval(['datestr(nov_s.L' num2str(itime) '.time_20hz(end))']);
disp(['Time Frame: ' t_initial ' to ' t_final])

% Save Figure
print(fig,'-depsc', [fig_path 'fig06.eps'], '-r300');                        % Resolution units: dpi (dots per inch); default resolution: 150 dpi 

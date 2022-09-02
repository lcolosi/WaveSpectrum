%% Figure 7: Omni-directional and saturation spectrograms for Stokes DELMAR2020 Experiment 
% Luke Colosi | lcolosi@ucsd.edu | August 20th, 2022

%-------------------------------- Caption --------------------------------%
% Observed omni-directional (a) and saturation (b) wave spectrograms 
% computed from measurement collected aboard the Stokes wave glider's 
% during the large box trajectory (07:57:22 to 08:42:02 September 10, 2020). 
%-------------------------------------------------------------------------%

clc, clearvars -except Ns, close all;

% Set text interpreter 
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
period = 5;                                                                 % Time duration for computing wave spectra
date_i = '09-Sep-2020 02:30:00';                                            % Start date 
date_f = '11-Sep-2020 16:10:00';                                            % End date
date_o = {date_i, date_f}; 

% Spectral parameters
fe = 20;                                                                    % Sampling rate of novatel (Hz)
dt = 1/fe;                                                                  % Time interval between measurements (s)
fn = (1/2)*fe;                                                              % Nyquist Frequency (Hz)
df = 0.01;                                                                  % Frequency resolution
nfft = fe/df;                                                               % Window length 
f = (0:df:fn);                                                              % Observed frequency (Hz)
dtheta = 5;                                                                 % Angular resolution (degrees)
ntheta = 360/dtheta+1;                                                      % Number of angles
lambda_c = 1.5613;                                                          % % Wavelength cutoff (meters)
f_noise = sqrt(g/(2*pi*lambda_c));                                          % Noise frequency cutoff 
toolbox = 'WAFO';                                                           % Method used to compute directional spectrum 
variables = 'heave_velocity';                                               % Heave and horizontal velocity are used to compute the direction spectrum
scaling = false;                                                            % Variance of directional spectrum is not scaled to match variance of heave spectrum 
freq_band = find(f > 0.02 & f < f_noise);                                   %Frequency band (for computing significant wave height)

%% Call Data

%--------- STOKES ---------% 
if isempty(whos('Ns'))
    Ns = load([ROOT vehicle '/NOVATEL_downsampled_20Hz_ALL_' vehicle '.mat']); % RAW downsampled GPS/IMU Data (20 Hz)
    
    % Transpose time_20hz field: 
    Ns.nov.time_20hz = Ns.nov.time_20hz';

end

%% Compute spectrogram 

% Set variables
field_20hz = {'heave'; 'VEL_east'; 'VEL_north'; 'VEL_up'; 'time_20hz'};

% Convert the 20 Hz time record into two 20 min increment time records 
T0 = unique( floor(Ns.nov.time_20hz*24*60/period)/(24*60/period) );
T0(T0 < datenum(date_i)) = []; 
T0(T0 > datenum(date_f)) = [];
T1 = T0 + 1/(24*60/period);

% Initialize counter
kr = 0;

% Loop through time increments
for is = 1:length(T0)
    
    % Create index for grabbing data from ith time interval
    It_s = find(Ns.nov.time_20hz >= T0(is) & Ns.nov.time_20hz < T1(is));
    
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
        eval(['nov_s.' name ' = Ns.nov.' name '(It_s,:);'])

    end
    
    % Compute the time elapsed from the first measurement in the time
    % interval for each time recorded in seconds
    nov_s.t = (nov_s.time_20hz - nov_s.time_20hz(1))*86400; 
    time_int_ns = 0:dt:round(max(nov_s.t,[],'omitnan'));
    
    % Check if more than 6000 ground speed measurements (10 mins) exists to
    % continue
    Inum_s = find(~isnan(nov_s.VEL_east)); 
    if length(Inum_s) < 6000
        continue
    end
    
    %%%%%%% 20 Hz data %%%%%%%
    
    % Loop through field names
    for k = 1:length(field_20hz)
        
        % Set field variable name
        name = char(field_20hz(k,:));
        
        % Find indices of non-nan data points
        eval(['Inum = find(~isnan(nov_s.' name '));'])
        
        % Interpolate time steps with nans
        eval(['nov_s.' name ' = interp1(nov_s.t(Inum),nov_s.' name '(Inum),time_int_ns,''linear'',''extrap'');'])
        
    end 
    
    % Set counter
    kr = kr + 1;
    
    % Compute directional spectrum  
    [nov_s.Sd(:,:,kr), nov_s.f, nov_s.theta] = compute_directional_spectrum(nov_s.heave, nov_s.VEL_east, nov_s.VEL_north, nov_s.VEL_up, nov_s.time_20hz, f, df, dtheta, nfft, fe, toolbox, variables, scaling, dir_con);

    % Compute mean wave direction 
    [nov_s.mwd(:,kr), ~, ~, ~] = dir_spec_parameters(nov_s.VEL_east, nov_s.VEL_north, nov_s.VEL_up, nfft, fe, dir_con); 
    
    % Compute the Omni-directional Spectra from directional wave spectra 
    nov_s.spectrogram_omni(:,kr) = sum(nov_s.Sd(:,:,kr) * dtheta, 1);

    % Compute the half way point of the time series
    nov_s.time_legs(kr) = T1(kr) - minutes(period/2);
    
end 

% Obtain indicies frequency below high frequency cut off
Inoise_s = find(nov_s.f <= f_noise & nov_s.f >= 0.02);

% Remove frequencies and beyond 1 Hz 
nov_s.f_ob = nov_s.f(Inoise_s);
nov_s.spectrogram_omni_f_ob = nov_s.spectrogram_omni(Inoise_s,:);
nov_s.mwd_f_ob = nov_s.mwd(Inoise_s,:);

% Compute equivalent saturation spectrum 
nov_s.sat_spectrogram_omni_f_ob = nov_s.spectrogram_omni_f_ob .* (nov_s.f_ob').^(5);

%% Plot omni-directional and saturation wave spectrogram
clc, close all; 

% Set power spectral densities equal to zero to NaN
idx_zero = nov_s.spectrogram_omni_f_ob == 0; 
nov_s.spectrogram_omni_f_ob(idx_zero) = NaN; 
nov_s.sat_spectrogram_omni_f_ob(idx_zero) = NaN; 

% Find the time steps closest to the intital and final time of the leg
% trajectory from figures 5 and 6. 
time_i = datenum('10-Sep-2020 07:57:22'); 
time_f = datenum('10-Sep-2020 08:42:02'); 
idx_int = find(T0 >= time_i); 
idx_fin = find(T1 <= time_f); 

% Set plotting parameters
t_ticks = datetime('09-Sep-2020 00:00:00'):hours(12):datetime('12-Sep-2020 00:00:00'); 
red = [0.6350 0.0780 0.1840]; 
blue = [0 0.4470 0.7410];
fontsize = 18;

% Create figure
figure('units','normalized','outerposition',[0 0 1 0.8])

%------------ Subplot 1 ------------% 
subplot(2,1,1)

% Plot omni-directional spectrogram
pc = pcolor(T0, nov_s.f_ob, nov_s.spectrogram_omni_f_ob);

hold on 
xline(T0(idx_int(1)),'-', 'LineWidth', 1.5, 'color', [0 0 0])
xline(T1(idx_fin(end)),'-', 'LineWidth', 1.5, 'color', [0 0 0])
hold on 

% Set figure attributes
title('(a)')
pc.EdgeColor = 'none';
ylabel('$f_{ob}$ (Hz)', 'Interpreter', 'latex')
xticks(datenum(t_ticks))
datetick('x', 'mmm dd', 'keepticks')
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlim([T0(1), T1(end)])
grid on
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca, 'TickLength', [0.007, 0.007]) 
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
annotation('rectangle',[0.678472222222223 0.581903276131046 0.195833333333333 0.341653666146646],...
           'LineWidth',1.5,'LineStyle','--');

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu')))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = 'S($t,f_{ob}$) (m$^2$ Hz$^{-1}$)';
caxis([10^-9, 10^-4]);
cb.Ticks = [10^-9; 10^-8; 10^-7; 10^-6; 10^-5; 10^-4] ;
cb.TickLabels = {'$10^{-9}$'; '$10^{-8}$'; '$10^{-7}$'; '$10^{-6}$'; '$10^{-5}$'; '$10^{-4}$'} ; 
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.03;
cb.FontSize = fontsize;

%------------ Subplot 2 ------------% 
subplot(2,1,2)

% Plot saturation spectrogram
pc = pcolor(T0, nov_s.f_ob, nov_s.sat_spectrogram_omni_f_ob);

hold on 
xline(T0(idx_int(1)),'-', 'LineWidth', 1.5, 'color', [0 0 0])
xline(T1(idx_fin(end)),'-', 'LineWidth', 1.5, 'color', [0 0 0])
hold on

% Set figure attributes
title('(b)')
pc.EdgeColor = 'none';
xlabel('UTC time from Sep 9$^{\textrm{th}}$, 2020', 'Interpreter', 'latex')
ylabel('$f_{ob}$ (Hz)', 'Interpreter', 'latex')
xticks(datenum(t_ticks))
datetick('x', 'mmm dd', 'keepticks')
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlim([T0(1), T1(end)])
grid on
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca, 'TickLength', [0.007, 0.007]) 
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')
annotation('rectangle',[0.678472222222222 0.110764430577223 0.195138888888889 0.34009360374415],...
           'LineWidth',1.5,'LineStyle','--');

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu')))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = '$f_{ob}^5 \cdot$ S($t,f_{ob}$) (m$^2$ Hz$^{4}$)';
caxis([10^-12, 10^-7.5]);
cb.Ticks = [ 10^-12; 10^-11; 10^-10; 10^-9; 10^-8] ;
cb.TickLabels = { '$10^{-12}$'; '$10^{-11}$'; '$10^{-10}$'; '$10^{-9}$'; '$10^{-8}$'} ; 
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.03;
cb.FontSize = fontsize;

% Save Figure
saveas(gcf, [fig_path 'figure_7.png'])

%% Plot the Mean Wave Direction 

% Create Figure and axes
figure('units','normalized','outerposition',[0 0 1 0.6])

% Plot MWD
pc = pcolor(T0, nov_s.f_ob, nov_s.mwd_f_ob);

% Set figure attributes
title('Mean Wave Direction', 'Interpreter', 'latex', 'FontSize', fontsize)
pc.EdgeColor = 'none';
xlabel('UTC time from Sep 9$^{\textrm{th}}$, 2020', 'Interpreter', 'latex')
ylabel('$f_{ob}$ (Hz)', 'Interpreter', 'latex')
xticks(datenum(t_ticks))
datetick('x', 'mmm dd', 'keepticks')
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlim([T0(1), T1(end)])
grid on
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')

% Set colorbar attributes
cb = colorbar;
colormap(cmocean('phase'))
cb.Label.Interpreter = 'Latex';
cb.Label.String = '$\overline{\theta}(t,f_{ob})$ (deg)';
caxis([0 360]);
cb.Ticks = linspace(0,360,5);
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.02;
cb.FontSize = fontsize;

% Save Figure
saveas(gcf, [fig_path 'mean_wave_direction_stokes_delmar2020.png'])
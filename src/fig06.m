%% Figure 6: Examples of omni-directional and saturation spectra 
% Luke Colosi | lcolosi@ucsd.edu | August 20th, 2022

%-------------------------------- Caption --------------------------------%
% Observed omni-directional (a) and saturation (b) wave spectra computed
% from the directional wave spectrum in Figure~\ref{f5}. The dashed and
% dot-dashed lines a have $f_{ob}^{-4}$ and $f_{ob}^{-5}$ spectral slopes
% respectively.
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
z_station = 1;                                                              % Height of weather station (meters)
g = 9.81;                                                                   % gravitational acceleration (m/s) 
r = 9.7 * 10^(-3);                                                          % Phillip's constant (dimensionless)
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

% Upload and process novatel and weather station data
[nov_s, w_s, nlegs_s] = process_wg_data(vehicle, ROOT, date_o, dir_con, [dt_n, dt_w]);

%% Compute Wave Spectra

% Set mean wave direction of high frequency waves (0.2 Hz <= f_ob <= 0.5 Hz)
mwd = mod(180 + 303, 360);

% Loop through legs
for n = 1:nlegs_s
    
    % Compute Directional Spectrum  
    eval(['[nov_s.Sd(:,:,n), nov_s.f, nov_s.theta] = compute_directional_spectrum(nov_s.L' num2str(n) '.heave, nov_s.L' num2str(n) '.VEL_east, nov_s.L' num2str(n) '.VEL_north, nov_s.L' num2str(n) '.VEL_up, nov_s.L' num2str(n) '.time_20hz, f, df, dtheta, nfft, fe_n, toolbox, variables, scaling, dir_con);']) 
    
    % Compute the Omni-directional Spectra from directional wave spectra 
    nov_s.spectrogram_omni(:,n) = sum(nov_s.Sd(:,:,n) * dtheta, 1);

    % Compute mean wind speed over time interval and friction velocity 
    eval(['w_s.mTWS(n) = mean(w_s.L' num2str(n) '.TWS);'])
    [w_s.mfv(n),~,~] = FUNC_ustar_Z(z_station,w_s.mTWS(n));

    % Compute transition frequencies
    nov_s.fst(n) = ((g*sqrt(r))/(2*pi*w_s.mfv(n)));                         % Equilibrium to saturation range frequency transition

    % Compute mean direction of platform 
    eval(['[nov_s.mD_legs(n), ~] = direction_stats(nov_s.L' num2str(n) '.true_course, dt_n, 0);']) % Convention: CW, going towards, ref north
    
    % Compute relative angle between platform and high wave direction
    nov_s.rel_theta_wave_legs(n) = mod(mwd - nov_s.mD_legs(n), 360);  % 1D-method (Convention: CW, going towards, ref north)

end

% Obtain indicies frequency below high frequency cut off
Inoise_s = find(nov_s.f <= f_noise);

% Remove frequencies and beyond 1 Hz 
nov_s.f_ob = nov_s.f(Inoise_s);
nov_s.spectrogram_omni_f_ob = nov_s.spectrogram_omni(Inoise_s,:);
nov_s.Sd_f_ob = nov_s.Sd(:,Inoise_s,:);

% Compute equivalent saturation spectrum 
nov_s.sat_spectrogram_omni_f_ob = nov_s.spectrogram_omni_f_ob .* (nov_s.f_ob').^(5);

%% Plot omni-directional and saturation wave spectra
clc, close all; 

% Set plotting parameters
red = [0.6350 0.0780 0.1840]; 
blue = [0 0.4470 0.7410];
fontsize = 24;
itime = 18; 

% Set variables 
fi = nov_s.fst(itime);                                                     
idx_f = find(nov_s.f_ob >= fi); 
Si = nov_s.spectrogram_omni_f_ob(idx_f(1),itime);                          

% Compute f^-5 slope lines
slope4 = Si*((nov_s.f_ob/fi).^(-4));
slope5 = Si*((nov_s.f_ob/fi).^(-5));

% Create figure
figure('units','normalized','outerposition',[0 0 1 1])

%------------- Subplot 1 -------------%
subplot(1,2,1);

% Plot Omni-directional Wave Spectrum
loglog(nov_s.f_ob, nov_s.spectrogram_omni_f_ob(:,itime), '-', 'LineWidth', 2, 'Color', blue);

hold on 
    pc1 = loglog(nov_s.f_ob,slope4, '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    pc2 = loglog(nov_s.f_ob,slope5, '-.', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
hold off

% Set figure Attributes
title('(a)')
ylabel('S$_{ob}(f_{ob})$  (m$^2$ Hz$^{-1}$)','Interpreter','Latex')
xlabel('$f_{ob}$ (Hz)','Interpreter','Latex')
xlim([10^-2 10^0])
ylim([10^-10 10^-4])
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
set(gca,'TickLabelInterpreter','latex')
axis square
grid on 
set(gca,'Box','on') 
set(gca, 'color', 'w')
legend([pc1, pc2], {'$f_{ob}^{-4}$', '$f_{ob}^{-5}$'}, 'Location', 'NorthEast', 'Fontsize', fontsize)

%------------- Subplot 2 -------------%
subplot(1,2,2);

% Plot Saturation Wave Spectrum
loglog(nov_s.f_ob, nov_s.sat_spectrogram_omni_f_ob(:,itime), '-', 'LineWidth', 2, 'Color', red)

% Set figure Attributes
title('(b)')
ylabel('B$_{ob}(f_{ob})$ (m$^2$ Hz$^{4}$)','Interpreter','Latex')
xlabel('$f_{ob}$ (Hz)','Interpreter','Latex')
xlim([10^-2 10^0])
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
set(gca,'TickLabelInterpreter','latex')
axis square
grid on 
set(gca,'Box','on') 
set(gca, 'color', 'w')

% Display time frame omni-directional-Spectrum is computed over
t_initial = eval(['datestr(nov_s.L' num2str(itime) '.time_20hz(1))']); 
t_final = eval(['datestr(nov_s.L' num2str(itime) '.time_20hz(end))']);
disp(['Time Frame: ' t_initial ' to ' t_final])

% Save Figure
saveas(gcf, [fig_path 'fig06.png'])

% Display pertinent variables
disp(['Friction velocity: ' num2str(w_s.mfv(itime))])
disp(['Equilibrium to saturation range frequency transition ' num2str(nov_s.fst(itime))])
disp(['Relative angle between mean wave direction for high frequency waves and platform heading: ' num2str(nov_s.rel_theta_wave_legs(itime))])


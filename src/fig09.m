%% Figure 9: Map Observed frequency wave spectrogram to Intrinsic frequency space
% Luke Colosi | lcolosi@ucsd.edu | August 20th, 2022

%-------------------------------- Caption --------------------------------%
% Saturation wave spectrogram for the Stokes Wave Glider small box
% trajectory in (a) observed frequency space and intrinsic frequency space
% using the (b) 1D- and (c) 2D-methods. The black curve is the
% $6\times10^{-10}$ m$^2$Hz$^4$ saturation spectral level contour where
% the correlation coefficient $r$ is computed. 
%-------------------------------------------------------------------------% 

clc, clear, close all;

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
z_station = 1;                                                              % Height of weather station (meters)
g = 9.81;                                                                   % gravitational acceleration (m/s) 
r = 9.7 * 10^(-3);                                                          % Phillip's constant (dimensionless)
dir_con = {'CW', 'gt', 'rn'};                                               % Directional conventions

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
[nov_s, w_s, nlegs_s] = process_wg_data(vehicle, ROOT, date_o, dir_con,...
                                        [dt_n, dt_w]);

%% Compute Wave Spectra
clc, close all 

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
    nov_s.fp(n) = nov_s.f(nov_s.spectrogram_omni(:,n) == max(nov_s.spectrogram_omni(:,n), [], 'omitnan'));  % Peak frequency
    nov_s.feq(n) = sqrt(2.25)*nov_s.fp(n);                                                                  % Spectral peak to equilibrium range frequency transition
    nov_s.fst(n) = ((g*sqrt(r))/(2*pi*w_s.mfv(n)));                                                         % Equilibrium to saturation range frequency transition

    % Compute mean direction of platform and wind
    eval(['[nov_s.mD_legs(n), ~] = direction_stats(nov_s.L' num2str(n) '.true_course, dt_n, 0);']) % Convention: CW, going towards, ref north
    eval(['[w_s.mTWD_legs(n), ~] = direction_stats(w_s.L' num2str(n) '.TWD, dt_w, 0);'])           % Convention: CW, going towards, ref north

    % Compute mean speed using instantaneous speed
    eval(['nov_s.mspeed_legs(n) = mean(nov_s.L' num2str(n) '.ground_speed_ms);'])

    % Compute mean speed using projected speed 
    eval(['nov_s.ground_speed_ms_proj = nov_s.L' num2str(n) '.ground_speed_ms .* cosd(nov_s.L' num2str(n) '.true_course - nov_s.mD_legs(n));'])
    nov_s.mspeed_proj(n) = mean(nov_s.ground_speed_ms_proj);
    
    % Compute relative angle between platform and wave direction
    nov_s.rel_theta_wind_legs(n) = mod(w_s.mTWD_legs(n) - nov_s.mD_legs(n), 360);  % 1D-method (Convention: CW, going towards, ref north)
    nov_s.rel_theta_legs(:,n) = mod(nov_s.theta - nov_s.mD_legs(n), 360);          % 2D-method (Convention: CW, going towards, ref north)

    % Compute the half way point of the time series
    eval(['nov_s.time_legs(n) = nov_s.L' num2str(n) '.time_20hz(round(length(nov_s.L' num2str(n) '.time_20hz) - (1/2)*length(nov_s.L' num2str(n) '.time_20hz)));']) 
    
    % Compute mean latitude
    eval(['nov_s.mlat_legs(n) = mean(nov_s.L' num2str(n) '.latitude);'])

end

% Obtain frequency indicies below high frequency cut off
Inoise_s = find(nov_s.f <= f_noise);

% Remove frequencies and beyond 1 Hz 
nov_s.f_ob = nov_s.f(Inoise_s);
nov_s.spectrogram_omni_f_ob = nov_s.spectrogram_omni(Inoise_s,:);
nov_s.Sd_f_ob = nov_s.Sd(:,Inoise_s,:);

% Compute equivalent saturation spectrum 
nov_s.sat_spectrogram_omni_f_ob = nov_s.spectrogram_omni_f_ob .* (nov_s.f_ob').^(5);

%% Account for Doppler shifts by mapping observed frequency to intrinsic frequency
clc, close all; 

% Set variables for mapping
tail = [zeros(size(nov_s.fst)); nov_s.feq; nov_s.fst]; 

%------ Map wave spectra ------%

% Loop through legs
for ileg = 1:nlegs_s
    
    % Map Omni-Directional Spectra
    [nov_s.spectrogram_omni_f_in_1d(:,ileg), nov_s.f_in_1d(:,ileg), nov_s.spectrogram_omni_f_ob(:,ileg), ~, nov_s.fb_1d(:,ileg), ~, ~, ~] = map_omni_dir_spectrum(nov_s.spectrogram_omni_f_ob(:,ileg), nov_s.f_ob, f_noise, df, nov_s.mspeed_proj(ileg), nov_s.rel_theta_wind_legs(ileg), tail(:,ileg));
    
    % Map Directional Spectra 
    [nov_s.dir_spectrogram_f_in(:,:,ileg), nov_s.f_in_2d(:,:,ileg), nov_s.dir_spectrogram_f_ob(:,:,ileg), ~, nov_s.fb_2d(:,:,ileg), ~, ~, ~]  = map_dir_spectrum(nov_s.Sd_f_ob(:,:,ileg)', nov_s.f_ob, f_noise, df, dtheta, nov_s.mspeed_proj(ileg), nov_s.rel_theta_legs(:,ileg));
    
    % Compute Omni-directional Spectra from directional spectra 
    spectrogram_omni_f_in_2d_nt = sum(nov_s.dir_spectrogram_f_in(:,:,ileg) * dtheta, 2);

    % Find non-nan values in S_fin after mapping
    idx_nan = ~isnan(spectrogram_omni_f_in_2d_nt); 

    % Attach a spectral tail if specified
    if tail(1,ileg) == true
        S_truc = spectrogram_omni_f_in_2d_nt(idx_nan);
        f_truc = nov_s.f_ob(idx_nan);
        f_tail = [f_truc(end), nov_s.f_ob(~idx_nan)];
        [nov_s.spectrogram_omni_f_in_2d(:,ileg), ~, fit, f_fit] = omnidir_spectral_tail(S_truc, f_truc, f_tail, tail(2,ileg), tail(3,ileg), f_noise);
    else
        nov_s.spectrogram_omni_f_in_2d(:,ileg) = spectrogram_omni_f_in_2d_nt;
    end

end

% Compute equivalent saturation spectrum for intrinsic frequency spectrum
nov_s.sat_spectrogram_omni_f_in_1d = nov_s.spectrogram_omni_f_in_1d .* (nov_s.f_ob').^(5); 
nov_s.sat_spectrogram_omni_f_in_2d = nov_s.spectrogram_omni_f_in_2d .* (nov_s.f_ob').^(5); 

%% Plot Observed and Intrinsic Saturation Spectrograms

% Set frequency low cutoff and latitude for small box cutoff
Ifreq_l = find(nov_s.f_ob > 0.03);                                          % Frequency low cutoff 
Ifreq_h = find(nov_s.f_ob > 0.03 & nov_s.f_ob < 0.4);                       % Frequency low and high cutoff
Ilat_s = find(nov_s.mlat_legs < 32.93);                                     % Indices for legs in small box

% Set plotting variables
[T,F] = meshgrid(nov_s.time_legs(Ilat_s), nov_s.f_ob(Ifreq_h));
contour_mod = 6*[1*10^(-10), 1*10^(-10)];
first = datetime(2020, 09, 10,23,0,0);
last = datetime(2020, 09, 11,13,0,0);
t_ticks = datenum(first:hours(1):last);
fontsize = 18;

% Create Figure and axes
figure('units','normalized','outerposition',[0 0 0.45 1])

%------------------- Observed Frequency Saturation Spectrogram -------------------%
ax1 = subplot(3,1,1);

% Plot Spectrogram
pc = pcolor(nov_s.time_legs(Ilat_s), nov_s.f_ob(Ifreq_l), nov_s.sat_spectrogram_omni_f_ob(Ifreq_l,Ilat_s));

hold on 
[ct1,~] = contour(T,F,nov_s.sat_spectrogram_omni_f_ob(Ifreq_h,Ilat_s), contour_mod, 'LineWidth', 2, 'LineColor', 'k');
hold off

% Set figure attributes
title('(a)')
pc.EdgeColor = 'none';
ylabel('$f_{ob}$ (Hz)', 'Interpreter', 'latex')
xlim([nov_s.time_legs(Ilat_s(1)), nov_s.time_legs(end)])
xticks(datenum(t_ticks(2:end)))
datetick('x', 'HH', 'keepticks')
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu', 100)))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = 'B($t,f$) (m$^2$ Hz$^{4}$)';
caxis([10^-12, 10^-7.5]);
cb.Ticks = [10^-12; 10^-11; 10^-10; 10^-9; 10^-8] ;
cb.TickLabels = { '$10^{-12}$'; '$10^{-11}$'; '$10^{-10}$'; '$10^{-9}$'; '$10^{-8}$'} ; 
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.02;
cb.FontSize = fontsize-2;
cb.Position = [0.860243103645498 0.270161290322581 0.0286457852433912 0.478648761961135]; 

%------------------- 1D-method Intrinsic Frequency Saturation Spectrogram -------------------%
ax2 = subplot(3,1,2);

% Plot Spectrogram
pc = pcolor(nov_s.time_legs(Ilat_s), nov_s.f_ob(Ifreq_l), nov_s.sat_spectrogram_omni_f_in_1d(Ifreq_l,Ilat_s));

hold on 
[ct2,~] = contour(T,F,nov_s.sat_spectrogram_omni_f_in_1d(Ifreq_h,Ilat_s), contour_mod, 'LineWidth', 2, 'LineColor', 'k');
hold off

% Set figure attributes
title('(b)')
pc.EdgeColor = 'none';
ylabel('$f_{in}$ (Hz)', 'Interpreter', 'latex')
xlim([nov_s.time_legs(Ilat_s(1)), nov_s.time_legs(end)])
xticks(datenum(t_ticks(2:end)))
datetick('x', 'HH', 'keepticks')
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu')))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = 'B$_{in}$($t,f_{in}$) (m$^2$ Hz$^{4}$)';
caxis([10^-12, 10^-7.5]);
cb.Ticks = [10^-12; 10^-11; 10^-10; 10^-9; 10^-8] ;
cb.TickLabels = { '$10^{-12}$'; '$10^{-11}$'; '$10^{-10}$'; '$10^{-9}$'; '$10^{-8}$'} ;  
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.02;
cb.FontSize = fontsize;

% Turn off colorbar
colorbar('off')

%------------------- 2D-method Intrinsic Frequency Saturation Spectrogram -------------------%
ax3 = subplot(3,1,3);

% Plot Spectrogram
pc = pcolor(nov_s.time_legs(Ilat_s), nov_s.f_ob(Ifreq_l), nov_s.sat_spectrogram_omni_f_in_2d(Ifreq_l,Ilat_s));

hold on 
[ct3,~] = contour(T,F,nov_s.sat_spectrogram_omni_f_in_2d(Ifreq_h,Ilat_s), contour_mod, 'LineWidth', 2, 'LineColor', 'k');
hold off

% Set figure attributes
title('(c)')
pc.EdgeColor = 'none';
xlabel('UTC time from Sep 11$^{\textrm{th}}$, 2020 (hrs)')
ylabel('$f_{in}$ (Hz)', 'Interpreter', 'latex')
xlim([nov_s.time_legs(Ilat_s(1)), nov_s.time_legs(end)])
xticks(datenum(t_ticks(2:end)))
datetick('x', 'HH', 'keepticks')
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca,'FontSize',fontsize)
set(gca,'TickLabelInterpreter','latex')

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu')))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = 'B$_{in}$($t,f_{in}$) (m$^2$ Hz$^{4}$)';
caxis([10^-12, 10^-7.5]);
cb.Ticks = [10^-12; 10^-11; 10^-10; 10^-9; 10^-8] ;
cb.TickLabels = { '$10^{-12}$'; '$10^{-11}$'; '$10^{-10}$'; '$10^{-9}$'; '$10^{-8}$'} ;  
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.02;
cb.FontSize = fontsize;

% Turn off colorbar
colorbar('off')

%-------- Set the position of the subplots --------%
% find current position [x,y,width,height]
pos1 = get(ax1,'Position');
pos2 = get(ax2,'Position');
pos3 = get(ax3,'Position');

% set width of the axes equal to a set width
width = 0.6986; 
pos1(3) = width; 
pos2(3) = width;
pos3(3) = width;
set(ax1,'Position',pos1)
set(ax2,'Position',pos2)
set(ax3,'Position',pos3)

% Save Figure
saveas(gcf, [fig_path 'figure_9.png'])

%% Compute correlation coefficient to quantify skill of the methods

% Set significance level and mean wave direction of high frequency waves (0.2 <= f_ob <= 0.5 Hz)
alpha = 0.05;
mwd = mod(180 + 303, 360);                                                  % Convert to directional convention CW, going towards, and reference north (see figure 7 code for calculation of mean wave direction 303). 

% Compute projected velocity for high frequency waves (0.2 <= f_ob <= 0.5 Hz)
nov_s.proj_speed_legs = nov_s.mspeed_legs.*cosd(mod(mwd - nov_s.mD_legs, 360));
 
% Set the frequency band to integrate/average over
idx_low = nov_s.f_ob >= 1*10^(-1);
idx_high = nov_s.f_ob <= 10*10^(-1);
Ifreq_band = logical(idx_low.*idx_high);

% Average over a frequency band in the saturation spectrum
av_sat_spec_fob = mean(nov_s.sat_spectrogram_omni_f_ob(Ifreq_band,Ilat_s),1,'omitnan');
av_sat_spec_fin_1d = mean(nov_s.sat_spectrogram_omni_f_in_1d(Ifreq_band,Ilat_s),1,'omitnan');
av_sat_spec_fin_2d = mean(nov_s.sat_spectrogram_omni_f_in_2d(Ifreq_band,Ilat_s),1,'omitnan');

% Compute the correlation between the projected velocity for high frequency waves and averaged unmapped and mapped saturation wave spectrogram
[r_av_sat_fob, pval_av_sat_fob] = corr(nov_s.proj_speed_legs(Ilat_s)',av_sat_spec_fob'); 
[r_av_sat_fin_1d, pval_av_sat_fin_1d] = corr(nov_s.proj_speed_legs(Ilat_s)',av_sat_spec_fin_1d');  
[r_av_sat_fin_2d, pval_av_sat_fin_2d] = corr(nov_s.proj_speed_legs(Ilat_s)',av_sat_spec_fin_2d');

% Display results 
disp('Mean Saturation Spectrum correlation with projected speed')
disp(['Correlation: Observed ' num2str(r_av_sat_fob) ', 1d method ' num2str(r_av_sat_fin_1d) ', and 2d method ' num2str(r_av_sat_fin_2d)])
disp(['Significance of correlation at 95% CL: Observed ' num2str(pval_av_sat_fob < alpha) ', 1d method ' num2str(pval_av_sat_fin_1d < alpha) ', and 2d method ' num2str(pval_av_sat_fin_2d < alpha)])

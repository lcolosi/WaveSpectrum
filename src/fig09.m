%% Figure 9: Map Observed frequency wave spectrogram to Intrinsic frequency space
% Luke Colosi | lcolosi@ucsd.edu | August 20th, 2022

%-------------------------------- Caption --------------------------------%
% Saturation wave spectrogram for the Stokes Wave Glider small box
% trajectory in (a) observed frequency space and intrinsic frequency space
% using the (b) 1D- and (c) 2D-methods. The black curve is the
% $6\times10^{-10}$ m$^2$Hz$^4$ saturation spectral level contour where
% the correlation coefficient $r$ is computed. 
%-------------------------------------------------------------------------% 

%--------------------------------- Notes ---------------------------------%
% (1) Explanation for what's going wrong with attaching a spectral tail:
% The spectral tail fitting works descently well when the blocking
% frequency is within the saturation range and I fit a line in log10 space
% from the equilibrium to saturation range transition frequency to
% the highest accurately resolved frequency, the spectral tail is pretty
% good. But when I try to fit the spectral tail in the equilibrium range,
% the spectral tail is really bad! There are two reasons for this: 
%
%   (a) The transition frequency from the spectral peak to the equilibrium
%   range is pretty low and close to the spectral peak. This transition
%   frequency's function form was used in Lenain et al. 2017 (references
%   Phillips 1985). I need to understand exactly where this comes from. 
% 
%   (b) As a result of (a), the linear fit includes a large band of the
%   wave spectrum. For the delmar experiment, there exist multiple wave
%   systems in the wave field (swell and developing wind seas). This biases
%   the linear fit and makes it least squares approach unaffective. 
% 
% To resolve this problem, reduce the range of frequencies that we preform 
% the least squares fit for the equilibrium range case. 
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
    [nov_s.spectrogram_omni_f_in_1d(:,ileg), nov_s.f_in_1d(:,ileg), nov_s.spectrogram_omni_f_ob(:,ileg), ~, nov_s.fb_1d(:,ileg), ~, ~, ~] = map_omni_dir_spectrum(nov_s.spectrogram_omni_f_ob(:,ileg), nov_s.f_ob, f_noise, df, nov_s.mspeed_legs(ileg), nov_s.rel_theta_wind_legs(ileg), tail(:,ileg));
    
    % Map Directional Spectra 
    [nov_s.dir_spectrogram_f_in(:,:,ileg), nov_s.f_in_2d(:,:,ileg), nov_s.dir_spectrogram_f_ob(:,:,ileg), ~, nov_s.fb_2d(:,:,ileg), ~, ~, ~]  = map_dir_spectrum(nov_s.Sd_f_ob(:,:,ileg)', nov_s.f_ob, f_noise, df, dtheta, nov_s.mspeed_legs(ileg), nov_s.rel_theta_legs(:,ileg), tail(:,ileg));
    
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

%     % Plot omni-directional spectra with and without tails
%     figure('units','normalized','outerposition',[0.6 0.5 0.4 0.6]);
%     loglog(nov_s.f_ob, spectrogram_omni_f_in_2d_nt, '-b', 'LineWidth', 2)
% 
%     % Set figure attributes
%     xlabel('f (Hz)')
%     ylabel('S(f) ($m^2 Hz^{-1}$)')
%     xlim([10^-2 10^0])
%     ylim([10^-10 10^-4])
%     set(gca, 'FontSize', 15)
% 
%     % Save Figure
%     saveas(gcf, '../../WaveGlider/figs/SMODE_P2021/omnidir_spec_tail_1.png')
% 
%     % Plot the least-square fit and the omni-directional spectrum with
%     % spectral tail 
%     hold on 
%     loglog(nov_s.f_ob, nov_s.spectrogram_omni_f_in_2d(:,ileg), ':r', 'LineWidth', 2)
%     loglog(f_fit, fit, '.-k', 'LineWidth', 2)
%     hold off

end

% Compute equivalent saturation spectrum for intrinsic frequency spectrum
nov_s.sat_spectrogram_omni_f_in_1d = nov_s.spectrogram_omni_f_in_1d .* (nov_s.f_ob').^(5); 
nov_s.sat_spectrogram_omni_f_in_2d = nov_s.spectrogram_omni_f_in_2d .* (nov_s.f_ob').^(5); 

%% Plot Observed and Intrinsic Saturation Spectrograms
clc

% Set frequency low cutoff and latitude for small box cutoff
Ifreq_l = find(nov_s.f_ob > 0.03);                                          % Frequency low cutoff 
Ifreq_h = find(nov_s.f_ob > 0.03 & nov_s.f_ob < 0.4);                       % Frequency low and high cutoff
Ilat_s = find(nov_s.mlat_legs < 32.93);                                     % Indices for legs in small box

% Set plotting variables
[T,F] = meshgrid(nov_s.time_legs(Ilat_s), nov_s.f_ob(Ifreq_h));
contour_mod = [6*10^(-10), 6*10^(-10)];
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

%% Compute correlation coefficient to quantify skill of the technique
clc

% Set significance level
alpha = 0.05;

% Clean contour data

%----------- Observed Saturation Spectrogram -----------% 

% Remove outlier contour values 
idx_time_ct1 = ct1(1,:) >= nov_s.time_legs(Ilat_s(1)); 
idx_freq_ct1 = ct1(2,:) <= 0.22; 
idx_ct1 = logical(idx_time_ct1 .* idx_freq_ct1); 
ct1_clean = ct1(:,idx_ct1); 

% Interpolate onto time interval of time_legs
ct1_interp = interp1(ct1_clean(1,:),ct1_clean(2,:),nov_s.time_legs(Ilat_s),'linear','extrap');

%----------- 1D-Method Saturation Spectrogram -----------% 

% Remove outlier contour values 
idx_time_ct2 = ct2(1,:) >= nov_s.time_legs(Ilat_s(1)); 
idx_freq_ct2 = ct2(2,:) <= 0.22; 
idx_ct2 = logical(idx_time_ct2 .* idx_freq_ct2); 
ct2_clean = ct2(:,idx_ct2); 

% Interpolate onto time interval of time_legs
ct2_interp = interp1(ct2_clean(1,:),ct2_clean(2,:),nov_s.time_legs(Ilat_s),'linear','extrap');

%----------- 2D-Method Saturation Spectrogram -----------% 

% Remove outlier contour values 
idx_time_ct3 = ct3(1,:) >= nov_s.time_legs(Ilat_s(1)); 
idx_freq_ct3 = ct3(2,:) <= 0.22; 
idx_ct3 = logical(idx_time_ct3 .* idx_freq_ct3); 
ct3_clean = ct3(:,idx_ct3); 

% Interpolate onto time interval of time_legs
ct3_interp = interp1(ct3_clean(1,:),ct3_clean(2,:),nov_s.time_legs(Ilat_s),'linear','extrap');

% Compute the correlation between the heading and integrated/averaged unmapped and mapped wave spectrogram
[r_ct_fob, pval_ct_fob] = corr(nov_s.mD_legs(Ilat_s)',ct1_interp'); 
[r_ct_fin_1d, pval_ct_fin_1d] = corr(nov_s.mD_legs(Ilat_s)',ct2_interp');  
[r_ct_fin_2d, pval_ct_fin_2d] = corr(nov_s.mD_legs(Ilat_s)',ct3_interp');

% Display results 
disp('Power Spectral density Contour Saturation Spectrum')
disp(['Correlation: Observed ' num2str(r_ct_fob) ', 1d method ' num2str(r_ct_fin_1d) ', and 2d method ' num2str(r_ct_fin_2d)])
disp(['Significance of correlation at 95% CL: Observed ' num2str(pval_ct_fob < alpha) ', 1d method ' num2str(pval_ct_fin_1d < alpha) ', and 2d method ' num2str(pval_ct_fin_2d < alpha)])

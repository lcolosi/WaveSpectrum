%% Compute Omni-directional wave spectra in terms of Intrinsic frequency
% Luke Colosi | lcolosi@ucsd.edu | June 5th, 2022

clc, clearvars -except N W, close all;

% Set the current working directory
cd ~/Desktop/projects/asi_lab_summer_internship/WaveGlider/src/

% Set text interpreter 
set(0,'defaultTextInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')

% Set path to data
%%%% Vehicle %%%%
vehicle = 'WHOI43';

%%%% Root %%%%
ROOT = '../data/SMODE_P2021/';

% Set path for figures
fig_path = '../figs/SMODE_P2021/';

%%%%%%%%%% Initial global variables %%%%%%%%%%

% Physical parameters
z_station = 1;                                                              % Height of weather station (meters)
g = 9.81;                                                                   % gravitational acceleration (m/s) 
r = 9.7 * 10^(-3);                                                          % Phillip's constant (dimensionless)
dir_con = {'CW', 'gt', 'rn'};                                               % Directional conventions

% Temporal parameters
period = 10 + 2; % Length of record to high pass filter Units: minutes
date_i = '29-Oct-2021 00:00:00';  % Start date 
date_f = '04-Nov-2021 00:00:00';  % End date;
date_o = {date_i, date_f}; 
%------- Vectornav ------- %
fe_n = 20;                                                                  % Sampling rate (Hz)
dt_n = 1/fe_n;                                                              % Time interval between measurements (s)
%------- Weather Station ------- %
fe_w = 1; 
dt_w = 1/fe_w;

% Spectral parameters                                                      
df = 0.01;                                                                  % Frequency resolution (Hz)
dtheta = 5;                                                                 % Angular resolution (degrees)
fn = (1/2)*fe_n;                                                            % Nyquist Frequency (Hz)
nfft = fe_n/df;                                                             % Window length 
f = (0:df:fn);                                                              % Observed frequency (Hz)
ntheta = 360/dtheta+1;                                                      % Number of angles
lambda_c = 1.5613;                                                          % Wavelength cutoff (meters)
f_noise = sqrt(g/(2*pi*lambda_c));                                          % Noise frequency cutoff 
toolbox = 'WAFO';                                                           % Method used to compute directional spectrum 
variables = 'heave_velocity';                                               % Heave and horizontal velocity are used to compute the direction spectrum
scaling = false;                                                            % Variance of directional spectrum is not scaled to match variance of heave spectrum  

%% Call Data

%--------- WHOI43 ---------%
if isempty(whos('N')) && isempty(whos('W'))
    N = load([ROOT vehicle '/SV3-1043_Vectornav.mat']);                     % RAW GPS/IMU Data from VectorNav (20 Hz)
    W = load([ROOT vehicle '/PLD2_1hz_ALL.mat']);                           % WXT/GILL Weather Station Data (1 hz)
    
    % Convert velocity down to velocity up via a sign change 
    N.vu = -N.vd; 
    
    % Set the time_20Hz data variable as the CR6 time
    N.time = N.time.cr6; 

end 

%% Process data

%--- Index data for specified time period and interpolating over NaNs ---%

% Set variables
field_20hz = {'alt'; 've'; 'vn'; 'vu'; 'lon'; 'lat'; 'time'; 'true_course'; 'ground_speed_ms'};
field_1hz = {'TWD'; 'TWS'; 'time'};

% Find indices for desired time frame
time_20hz_idx = find(N.time >= datenum(date_o(1)) & N.time < datenum(date_o(2)));
time_1hz_idx = find(W.PLD2_1hz.time >= datenum(date_o(1)) & W.PLD2_1hz.time < datenum(date_o(2)));

% Index time variables for interpolation 
n.time = N.time(time_20hz_idx);
w.time = W.PLD2_1hz.time(time_1hz_idx);

% Create the structures consisting of data from specified time frame

%%%%%%% 20 Hz data %%%%%%%

% Loop through fieldnames (except true_course and ground_speed_ms )
for k = 1:length(field_20hz) - 2

    % Set filenames
    name = char(field_20hz(k,:));

    % Index data
    eval(['n.' name ' = N.' name '(time_20hz_idx);'])
    
    % Compute the time elapsed from the first measurement in the 30 min
    % interval for each time recorded (units: seconds) 
    n.t = (n.time - n.time(1))*86400;
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

%--- Computing platform speed and direction ---% 
[n.ground_speed_ms, n.true_course] = platform_speed_dir(n.ve, n.vn, dir_con); 

%% Compute spectrogram and platform propagation speed and direction over specified time intervals
clc, close all;

% Convert the 20 Hz time record into two time records with increments
% specified by period beginning at date_i and ending at date_f
dt_defined = period*60; % Units: seconds
T0_i= 0:dt_defined:max(time_int_n); % Units: seconds

% Define the 1 minute time offset for the high pass filtered data
time_offset = 60; % Units: seconds
np_n = fe_n*time_offset; % Number of non-weather measurements taken in time offset duration
np_w = fe_w*time_offset; % Number of weather measurements taken in time offset duration

% Extend T0 to account for the offset from the truncation
i = length(T0_i) - 1 + ((5/2)*(dt_defined/time_offset) + 1);  
T0 = [T0_i(1:end-1) T0_i(end):dt_defined:T0_i(end)+(i+(i-2))*time_offset];

% Loop through legs
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
    
    % Index non-datenum time series
    fast.ttime_n = time_int_n(It_n);
    
    % Set transfer function coefficients for high pass butterworth filter 
    [Bcut,Acut] = butter(4,df/fn,'high');
    
    %%%%%%% 20 Hz data %%%%%%%
    
    % Loop through field names (except lon, lat, time) 
    for k = 1:(length(field_20hz) - 5)
        
        % Set field variable name
        name = char(field_20hz(k,:));
        
        % Detrend and high pass filter data
        eval(['hp = filtfilt(Bcut,Acut,detrend(double(fast.' name ')));'])
        
        % Truncate high pass time series by one minute on either side to avoid
        % boundary effects from filter
        eval(['fast.high_pass.' name ' = hp(np_n:end-np_n);'])
        
    end 
    
    % Truncate time series, non-datenum time series and platform ground speed and direction
    fast.high_pass.time = fast.time(np_n:end-np_n);
    fast.high_pass.ttime_n = fast.ttime_n(np_n:end-np_n);
    fast.high_pass.ground_speed_ms = fast.ground_speed_ms(np_n:end-np_n);
    fast.high_pass.true_course = fast.true_course(np_n:end-np_n);
    
    %%%%%%% 1 Hz data %%%%%%%
    
    % Loop through field names  
    for k = 1:length(field_1hz)
        
        % Set field variable name
        name = char(field_1hz(k,:));
        
        % Truncate time series by one minute on either side to be
        % consistent with high pass data
        eval(['fast.high_pass.' name ' = fast. ' name '(np_w:end-np_w);'])
        
    end
    
    % Compute mean wind speed over time interval and friction velocity 
    fast.high_pass.mTWS = mean(fast.high_pass.TWS);
    [fast.high_pass.mfv(is),~,~] = FUNC_ustar_Z(z_station,fast.high_pass.mTWS);
    
    %Compute Directional Spectrum 
    [n.Sd(:,:,is), n.f, n.theta] = compute_directional_spectrum(fast.high_pass.alt, fast.high_pass.ve, fast.high_pass.vn, fast.high_pass.vu, fast.high_pass.ttime_n, f, df, dtheta, nfft, fe_n, toolbox, variables, scaling, dir_con);

    % Compute the Omni-directional Spectra from directional wave spectra 
    n.spectrogram_omni(:,is) = sum(n.Sd(:,:,is) * dtheta, 1);
    
    % Compute heave spectra
    [n.spectrogram_heave(:,is), ~, ~] = pwelch(fast.high_pass.alt,nfft,nfft/2,nfft,fe_n);
    
    % Compute transition frequencies
    n.fp(is) = n.f(n.spectrogram_omni(:,is) == max(n.spectrogram_omni(:,is), [], 'omitnan')); % Peak frequency
    n.feq(is) = sqrt(2.25)*n.fp(is);                                                          % Spectral peak to equilibrium range frequency transition
    n.fst(is) = ((g*sqrt(r))/(2*pi*fast.high_pass.mfv(is)));                                  % Equilibrium to saturation range frequency transition

    % Compute Mean and standard deviation for Direction of platform
    [n.mD(is), n.std_D(is)] = direction_stats(fast.high_pass.true_course, dt_n, false); % Convention: CW, going towards, ref north
    
    % Compute mean speed using instantaneous velocity and projected
    % velocity
    n.mspeed(is) = mean(fast.high_pass.ground_speed_ms);
    n.stdspeed(is) = std(fast.high_pass.ground_speed_ms);
    n.ground_speed_ms_proj = fast.high_pass.ground_speed_ms .* cosd(fast.high_pass.true_course - n.mD(is));
    n.mspeed_proj(is) = mean(n.ground_speed_ms_proj);
    
    % Compute relative angle between platform and wave direction
    n.rel_theta(:,is) = mod(n.theta - n.mD(is), 360);                    % Convention: CW, going towards, ref wave glider mean direction
    n.rel_theta_std(is) = n.std_D(is)^2;
    
    % Compute the half way point of the time series
    T_leg = n.time(It_n);
    T_half = T_leg(round(length(It_n) - (1/2)*length(It_n)));
    
    % Floor time step 
    n.time_legs(is) = floor(T_half*1440)/1440;
    
    % Display center interval 
    disp(datestr(n.time_legs(is)))
    
end

% Obtain indicies frequency below high frequency cut off
Inoise = find(n.f <= f_noise);

% Remove frequencies beyond 1 Hz from correction
n.f_ob = n.f(Inoise);

% Remove spectral densities beyond 1 Hz from Omni-directional and directional spectra  
n.spectrogram_heave_f_ob = n.spectrogram_heave(Inoise,:);
n.spectrogram_omni_f_ob = n.spectrogram_omni(Inoise,:);
n.Sd_f_ob = n.Sd(:,Inoise,:);

% Compute saturation spectrogram
n.sat_spectrogram_heave_f_ob = n.spectrogram_heave_f_ob .* (n.f_ob').^(5);
n.sat_spectrogram_omni_f_ob = n.spectrogram_omni_f_ob .* (n.f_ob').^(5);

%% Correct wave spectra by mapping observed frequency to intrinsic frequency
clc, close all; 

% Set variables for mapping 
tail = [zeros(size(n.fst)); n.feq; n.fst];
   
% Doppler correct wave spectra

% ----------------- Instantaneous Speed ----------------- % 

% Loop through legs 
for i = 1:(length(T0)-1 ) 
    
    % Display center interval 
    disp(datestr(n.time_legs(i)))
      
    % Correct Directional Spectra 
    [n.inst_speed.dir_spectrogram_f_in(:,:,i), n.inst_speed.f_in_2d(:,:,i), n.inst_speed.dir_spectrogram_f_ob(:,:,i), ~, n.inst_speed.fb(:,:,i), ~, n.inst_speed.J(:,:,i), n.inst_speed.var_2d(:,i)]  = doppler_correct_dir_spectrum(n.Sd_f_ob(:,:,i)', n.f_ob, f_noise, df, dtheta, n.mspeed(i), n.rel_theta(:,i), tail(:,i)); 
    
    % Compute Omni-directional Spectra from directional spectra 
    n.inst_speed.spectrogram_omni_f_in(:,i) = sum(n.inst_speed.dir_spectrogram_f_in(:,:,i) * dtheta, 2);
    
end

% Compute equivalent saturation spectrum for corrected spectra

% --------- Instant Speed --------- % 
n.inst_speed.sat_spectrogram_omni_f_in = n.inst_speed.spectrogram_omni_f_in .* (n.f_ob').^(5); 

% Save spectrograms to mat file
time = n.time_legs; f = n.f_ob; spec_ob = n.sat_spectrogram_omni_f_ob; spec_in = n.inst_speed.sat_spectrogram_omni_f_in;
save([ROOT vehicle '/SMODE2021_spec_' vehicle '.mat'], 'time', 'f', 'spec_ob', 'spec_in')

%% Compute frequency band integrated/averaged unmapped and mapped wave spectrogram 
clc

% Set significance level
alpha = 0.05; 

% Set the frequency band to integrate/average over
idx_low = n.f_ob >= 2*10^(-1);
idx_high = n.f_ob <= 5*10^(-1);
Ifreq_band = logical(idx_low.*idx_high);

% Average over a frequency band 
av_sat_spec_fob = mean(n.sat_spectrogram_omni_f_ob(Ifreq_band,:),1,'omitnan');
av_sat_spec_fin = mean(n.inst_speed.sat_spectrogram_omni_f_in(Ifreq_band,:),1,'omitnan');

% Compute the integral over frequency band for observed spectrogram: 
int_spec_fob = trapz(n.f_ob(Ifreq_band),n.spectrogram_omni_f_ob(Ifreq_band,:),1); 

% Loop through legs 
for ileg = 1:(length(T0)-1)
    
    % Find the indices of NaN values for 1
    idx_nan_2d = ~isnan(n.inst_speed.spectrogram_omni_f_in(:,ileg)); 
    
    % Find the indices where frequency band and non-nans overlap
    idx_2d = logical(idx_nan_2d .* Ifreq_band'); 
    
    % Integrate over a frequency band for intrinsic frequency spectrograms:
    int_spec_fin(:,ileg) = trapz(n.f_ob(idx_2d),n.inst_speed.spectrogram_omni_f_in(idx_2d,ileg));

end

% Compute the correlation between the heading and integrated/averaged unmapped and mapped wave spectrogram
[r_av_fob, pval_av_fob] = corr(n.mD',av_sat_spec_fob'); 
[r_av_fin, pval_av_fin] = corr(n.mD',av_sat_spec_fin'); 
[r_int_fob, pval_int_fob] = corr(n.mD',int_spec_fob'); 
[r_int_fin, pval_int_fin] = corr(n.mD',int_spec_fin'); 

% Compute correlation with my own function 
r_mf = corr_coef(n.mD', av_sat_spec_fob');

% Compute the coefficient of determination 
[r2_av_fob] = coef_det(n.mD,av_sat_spec_fob);
[r2_av_fin] = coef_det(n.mD,av_sat_spec_fin); 
[r2_int_fob] = coef_det(n.mD,int_spec_fob); 
[r2_int_fin] = coef_det(n.mD,int_spec_fin); 

% Display results 
disp('Averaged Saturation Spectrum')
disp(['Correlation: before ' num2str(r_av_fob) ' and after ' num2str(r_av_fin)])
disp(['Significance of correlation at 95% CL: before ' num2str(pval_av_fob < alpha) ' and after ' num2str(pval_av_fin < alpha) ' mapping'])
disp(['R^2: before ' num2str(r2_av_fob) ' and after ' num2str(r2_av_fin)])
disp(' ')
disp('Integrated Saturation Spectrum')
disp(['Correlation: before ' num2str(r_int_fob) ' and after ' num2str(r_int_fin)])
disp(['Significance of correlation at 95% CL: before ' num2str(pval_int_fob < alpha) ' and after ' num2str(pval_int_fin < alpha) ' mapping'])
disp(['R^2: before ' num2str(r2_int_fob) ' and after ' num2str(r2_int_fin)])

% Initialize variable
freqs = [0.1, 0.2, 0.3, 0.4, 0.5];
r_fob = zeros(size(freqs)); r_fin = zeros(size(freqs)); 
pval_fob = zeros(size(freqs)); pval_fin = zeros(size(freqs)); 
r_squared_fob = zeros(size(freqs)); r_squared_fin = zeros(size(freqs));

% Loop through frequencies
for ifreq = 1:length(freqs)
    
    %Find index of ith frequency 
    idx_f = n.f_ob == freqs(ifreq);

    % Compute correlation and coefficient of determination before mapping
    [r_fob(ifreq),pval_fob(ifreq)]  = corr(n.mD',n.sat_spectrogram_omni_f_ob(idx_f,:)');
    r_squared_fob(ifreq) = coef_det(n.mD,n.sat_spectrogram_omni_f_ob(idx_f,:));

    % Compute correlation and coefficient of determination after mapping
    [r_fin(ifreq),pval_fin(ifreq)]  = corr(n.mD',n.inst_speed.spectrogram_omni_f_in(idx_f,:)');
    r_squared_fin(ifreq) = coef_det(n.mD,n.inst_speed.spectrogram_omni_f_in(idx_f,:));

    % Display results 
    disp(' ')
    disp(['Frequency: ' num2str(freqs(ifreq)) ])
    disp(['Correlation: before ' num2str(r_fob(ifreq)) ' and after ' num2str(r_fin(ifreq))])
    disp(['Significance of correlation at 95% CL: before ' num2str(pval_fob(ifreq) < alpha) ' and after ' num2str(pval_fin(ifreq) < alpha) ' mapping'])
    disp(['Coefficient of determination: before ' num2str(r_squared_fob(ifreq)) ' and after '  num2str(r_squared_fin(ifreq))])
     
end

%% Plot Observed and Intrinsic Frequency Saturation Spectra 
clc, close all; 

% Set plotting variables
first = datetime(2021, 10, 29,0,0,0);
last = datetime(2021, 11, 04,0,0,0);
t_ticks = datenum(first:days(1):last);

% Set frequency low cutoff and latitude for small box cutoff
Ifreq = find(n.f_ob > 0.02);  %Frequency low cutoff 

% Create Figure and axes
figure('units','normalized','outerposition',[0 0 1 1])

%------------------- Observed Frequency Saturation Spectrogram -------------------%
ax1 = subplot(2,1,1);

% Plot Spectrogram
pc = pcolor(n.time_legs, n.f_ob(Ifreq), n.sat_spectrogram_omni_f_ob(Ifreq,:));

% Set figure attributes
pc.EdgeColor = 'none';
ylabel('$f_{ob} \;(Hz)$', 'Interpreter', 'latex')
xlim([datenum(first), n.time_legs(end)])
datetick('x', 'mm/dd', 'keeplimits')
%xlim([datenum(2021,11,1,0,0,0), datenum(2021,11,2,0,0,0)])
%datetick('x', 'HH:MM', 'keeplimits')
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca,'FontSize',fontsize-2)
set(gca,'TickLabelInterpreter','latex')

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu')))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = '$f_{ob}^{5} \cdot S(f_{ob}) \;(m^2 Hz^{4})$';
caxis([10^-7, 2*10^-3]);
cb.Ticks = [10^-7; 10^-6; 10^-5; 10^-4; 10^-3];
cb.TickLabels = {'$10^{-7}$'; '$10^{-6}$'; '$10^{-5}$'; '$10^{-4}$'; '$10^{-3}$'}; 
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.02;
cb.FontSize = 18;

%------------------- Intrinsic Frequency Saturation Spectrogram -------------------%
ax2 = subplot(2,1,2);

% Plot Spectrogram
pc1 = pcolor(n.time_legs, n.f_ob(Ifreq), n.inst_speed.sat_spectrogram_omni_f_in(Ifreq,:));

hold on 
pc2 = plot(n.time_legs,nanmin(squeeze(n.inst_speed.fb(:,2,:)),[],1), '-k', 'LineWidth', 1);
% pc3 = plot(n.time_legs,n.fst, '-', 'color', 'w', 'LineWidth', 1);
hold off

% Set figure attributes
pc1.EdgeColor = 'none';
xlabel('UTC time from Oct. $29^{\textrm{th}}$, $2021$', 'Interpreter', 'latex')
ylabel('$f_{in} \;(Hz)$', 'Interpreter', 'latex')
xlim([datenum(first), n.time_legs(end)])
datetick('x', 'mm/dd', 'keeplimits')
%xlim([datenum(2021,11,1,0,0,0), datenum(2021,11,2,0,0,0)])
%datetick('x', 'HH:MM', 'keeplimits')
set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca,'FontSize',fontsize-2)
set(gca,'TickLabelInterpreter','latex')

% Set colorbar attributes
cb = colorbar;
colormap(flipud(cbrewer2('RdYlBu')))
set(gca,'ColorScale','log')
cb.Label.Interpreter = 'Latex';
cb.Label.String = '$f_{in}^{5} \cdot S(f_{in})\cdot \frac{df_{ob}}{df_{in}} \;(m^2 Hz^{4})$';
caxis([10^-7, 2*10^-3]);
cb.Ticks = [10^-7; 10^-6; 10^-5; 10^-4; 10^-3];
cb.TickLabels = {'$10^{-7}$'; '$10^{-6}$'; '$10^{-5}$'; '$10^{-4}$'; '$10^{-3}$'}; 
cb.TickLabelInterpreter = 'latex';
cb.TickDirection = 'out';
cb.TickLength = 0.02;
cb.FontSize = 18;

% Save Figure
saveas(gcf, [fig_path 'whoi_wg_smode_pilot_map_sat_specgram_heave_velocity_filter.png'])

%% Plot heading and wave glider speed
close all

% Set plotting variables
t_ticks = datetime(date_i):days(1):datetime(date_f);

% Create Figure and axes
fig = figure('units','normalized','outerposition',[0 0 0.7 1]);

%------------- Subplot 1 -------------%
ax1 = subplot(2,1,1);

% Plot the Wind direction
pc1 = plot(n.time_legs,n.mD, '.', 'Color', [0 0.4470 0.7410], 'MarkerSize', 5); 

% Set figure attributes
ylabel('$\textrm{Direction} (\textrm{Going Towards}\; ^{\circ})$', 'Interpreter', 'latex')
xlim([n.time_legs(1), n.time_legs(end)])
ylim([0,360])
xticks(datenum(t_ticks))
yticks([0,90,180,270,360])
datetick('x', 'mm/dd', 'keepticks')
set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
title('$\bf{(a)}$', 'FontSize', 15)

%------------- Subplot 2 -------------%
ax2 = subplot(2,1,2);

% Plot Temporally computed Significant Wave Height 
plot(n.time_legs, n.mspeed, '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.5)
 
% Set figure attributes
xlabel('UTC time from Oct. $29^{\textrm{th}}$, $2021$', 'Interpreter', 'latex')
ylabel('$| \vec{u} | \;(ms^{-1})$', 'Interpreter', 'latex')
xlim([n.time_legs(1), n.time_legs(end)])
ylim([0, 1.5])
xticks(datenum(t_ticks))
datetick('x', 'mm/dd', 'keepticks')
set(gca,'FontSize',12)
set(gca,'TickLabelInterpreter','latex')
set(ax2, 'box', 'on', 'Visible', 'on')
title('$\bf{(b)}$', 'FontSize', 15)

% Save Figure
saveas(fig, [fig_path 'smode_pilot_platform_motion.png'])

%% 
close all 
figure('units','normalized','outerposition',[0 0 1 1])

for iplot = 500:520
hold on 
loglog(n.f_ob(Ifreq), n.inst_speed.sat_spectrogram_omni_f_in(Ifreq,iplot), '-b', 'LineWidth', 2)
hold off
xlabel('$f_{obs} \;\textrm{or}\; f_{int} $')
ylabel('$S(f_{int})$')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
title()

% Display center interval 
disp(datestr(n.time_legs(iplot)))

pause 

end

%% Plot frequency band integrated and averaged power spectral density as a scatter plot before and after mapping 
close all

idx_sbox = find(n.time_legs > datenum(2021,11,01) & n.time_legs < datenum(2021,11,02)); % Indices for small box

%------------------- Average -------------------%

% Create Figure and axes
figure('units','normalized','outerposition',[0 0 1 0.6])

%------------------- f_ob Saturation Spectrogram -------------------%
plot(n.time_legs, av_sat_spec_fob, 'LineWidth', 1, 'Color', [0.6350 0.0780 0.1840]);

%------------------- f_in_1d Saturation Spectrogram -------------------%
hold on 
plot(n.time_legs, av_sat_spec_fin, 'LineWidth', 1, 'Color', [0 0.4470 0.7410]);
hold off

% Set figure attributes
ylabel('$I(t) \;(m^2 Hz^{4})$', 'Interpreter', 'latex')
xlim([n.time_legs(1), n.time_legs(end)])
ylim([0, 3*10^-3])  %ylim([3*10^-5, 5*10^-3])
xticks(datenum(t_ticks(2:end)))
datetick('x', 'mm/dd', 'keepticks')
%set(gca,'Yscale','log')
set(gca,'TickDir','out');
set(gca,'FontSize',fontsize-2)
set(gca,'TickLabelInterpreter','latex')
legend({'Observed Frequency', 'Intrinsic frequency'})

% Save Figure
saveas(gcf, [fig_path 'freq_average_saturation_spectrogram_heave_velocity_filter.png'])

%------------------- Integral -------------------%

% Create Figure and axes
figure('units','normalized','outerposition',[0 0 1 0.6])

%------------------- f_ob Spectrogram -------------------%
plot(n.time_legs, int_spec_fob, 'LineWidth', 1, 'Color', [0.6350 0.0780 0.1840]);

%------------------- f_in_2d Spectrogram -------------------%
hold on 
plot(n.time_legs, int_spec_fin, 'LineWidth', 1, 'Color', [0 0.4470 0.7410]);
hold off

% Set figure attributes
ylabel('$I(t) \;(m^2)$', 'Interpreter', 'latex')
xlim([n.time_legs(1), n.time_legs(end)])
%ylim([0.01, 0.07])
xticks(datenum(t_ticks(2:end)))
datetick('x', 'mm/dd', 'keepticks')
set(gca,'TickDir','out');
set(gca,'FontSize',fontsize-2)
set(gca,'TickLabelInterpreter','latex')
legend({'Observed Frequency', 'Intrinsic frequency'})

% Save Figure
saveas(gcf, [fig_path 'freq_integral_spectrogram_heave_velocity_filter.png'])

%--- Scatter plot ---%
% Create figure axes 
fig = figure('units','normalized','outerposition',[0 0 1 1]);

% Plot a scatter plot of the U and V current at mooring PS30M
scatter(av_sat_spec_fob(idx_sbox),n.mD(idx_sbox),15,'b', 'filled')

% Set figure attributes
xlabel('$I(t) \;(m^2 Hz^{4})$', 'Interpreter', 'latex')
ylabel('$\overline{\theta} \;(degrees)$', 'Interpreter', 'latex')
%xlim([-1 1])
%ylim([-1 1])
grid on
set(gca,'TickDir','out');
set(gca,'FontSize',14)
set(gca,'TickLabelInterpreter','latex')

% Save Figure
%saveas(fig, [fig_path 'figure_5.png'])

figure(5);
subplot(2,1,1)
plot(n.time_legs, av_sat_spec_fob, '-b')
subplot(2,1,2)
plot(n.time_legs, n.mD, 'r')

%% Developmental Code 
% %% Plot transition frequencies 
% close all 
% 
% % Set variables for plotting 
% POS = [100 100 900 900];
% 
% %  Create Figure 
% figure('Name', 'Bin Averaging check')
% set(gcf,'color',[1 1 1])
% set(gcf,'position',POS)
% 
% % Plot
% semilogy(n.time_legs, n.feq, '.-r', 'Markersize', 8)
% hold on 
% semilogy(n.time_legs, n.fst, '.-b', 'Markersize', 8)
% semilogy(n.time_legs, squeeze(n.inst_speed.fb(1,2,:)), '.-g', 'Markersize', 8)
% yline(f_noise, '--k', 'LineWidth', 1.5, 'Label', '$f_{noise}$', 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'FontSize',18)
% hold off
% 
% % Set figure attributes 
% xlabel('UTC time')
% ylabel('$f_{eq} \; \textrm{or} \; f_{sat}$ (Hz)', 'Interpreter','latex')
% datetick('x','mm/dd','keeplimits')
% set(gca,'FontSize',18)

% fe = 20;  % Sampling rate (Hz) of the Wave Glider (non-weather)
% fe_w = 1;  % Sampling rate (Hz) of the Wave Glider (weather)
% dt = 1/fe;  % Time interval between measurements (s) (non-weather)
% dt_w = 1/fe_w;  % Time interval between measurements (s) (weather)
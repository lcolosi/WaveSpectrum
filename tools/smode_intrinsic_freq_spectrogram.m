%% Compute Omni-directional wave spectra in terms of Intrinsic frequency
% Luke Colosi | lcolosi@ucsd.edu | June 5th, 2022

clc, clearvars -except N W, close all;

% Set text interpreter 
set(0,'defaultTextInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')

% Set path to data
%%%% Vehicle %%%%
vehicle = 'WHOI43';

%%%% Root %%%%
ROOT = '../data/SMODE2021/';

% Set path for figures
fig_path = '../figs/SMODE2021/';

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
dt_defined = period*60;                                                     % Units: seconds
T0_i= 0:dt_defined:max(time_int_n);                                         % Units: seconds

% Define the 1 minute time offset for the high pass filtered data
time_offset = 60;                                                           % Units: seconds
np_n = fe_n*time_offset;                                                    % Number of non-weather measurements taken in time offset duration
np_w = fe_w*time_offset;                                                    % Number of weather measurements taken in time offset duration

% Extend T0 to account for the offset from the truncation
i = 865;                                                                    % Old code used to compute time offset (hard coded this instead): length(T0_i) - 1 + ((5/2)*(dt_defined/time_offset) + 1);  
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

    % Compute the Omnidirectional Spectra from directional wave spectra 
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
   
%------ Map wave spectra ------%

% Loop through legs 
for i = 1:(length(T0)-1 ) 
      
    % Correct Directional Spectra 
    [n.dir_spectrogram_f_in(:,:,i), n.f_in_2d(:,:,i), n.dir_spectrogram_f_ob(:,:,i), ~, n.fb(:,:,i), ~, n.J(:,:,i), n.var_2d(:,i)]  = map_dir_spectrum(n.Sd_f_ob(:,:,i)', n.f_ob, f_noise, df, dtheta, n.mspeed_proj(i), n.rel_theta(:,i)); 
    
    % Compute Omni-directional Spectra from directional spectra 
    spectrogram_omni_f_in_nt = sum(n.dir_spectrogram_f_in(:,:,i) * dtheta, 2);

    % Find non-nan values in S_fin after mapping
    idx_nan = ~isnan(spectrogram_omni_f_in_nt);

    % Attach a spectral tail if specified
    if tail(1,i) == true
        S_truc = spectrogram_omni_f_in_nt(idx_nan);
        f_truc = n.f_ob(idx_nan);
        f_tail = [f_truc(end), n.f_ob(~idx_nan)];
        [n.spectrogram_omni_f_in(:,i), ~, fit, f_fit] = omnidir_spectral_tail(S_truc, f_truc, f_tail, tail(2,i), tail(3,i), f_noise);
    else
        n.spectrogram_omni_f_in(:,i) = spectrogram_omni_f_in_nt;
    end
    
end

% Compute equivalent saturation spectrum for corrected spectra
n.sat_spectrogram_omni_f_in = n.spectrogram_omni_f_in .* (n.f_ob').^(5); 

% Save spectrograms to mat file
time = n.time_legs; f = n.f_ob; spec_ob = n.sat_spectrogram_omni_f_ob; spec_in = n.sat_spectrogram_omni_f_in;
save([ROOT vehicle '/SMODE2021_spec_' vehicle '.mat'], 'time', 'f', 'spec_ob', 'spec_in')

%% Compute the mean platform speed and minimum blocking frequency 

% Set the time period to average over
idx_i = n.time_legs >= datenum(2021, 11, 1);
idx_f = n.time_legs <= datenum(2021, 11, 2);
Itime = logical(idx_i.*idx_f);

% Compute mean platform speed
mean_plat_speed = mean(n.mspeed_proj(Itime));
std_plat_speed = std(n.mspeed_proj(Itime)); 
stde_plat_speed = std_plat_speed/sqrt(length(n.mspeed_proj(Itime)));
max_plat_speed = max(n.mspeed_proj(Itime));

% Compute minimum blocking frequency
min_fb_ts = nanmin(squeeze(n.fb(:,2,:)),[],1);
min_fb = min(min_fb_ts(Itime));

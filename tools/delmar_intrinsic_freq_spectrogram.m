%% Compute Intrinsic frequency spectrum for DELMAR2020 experiment 
% Luke Colosi | lcolosi@ucsd.edu | August 22th, 2022  

clc, clearvars -except Ns Ws, close all;

% Set text interpreter 
set(0,'defaultTextInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')

% Set path to data
%%%% Vehicle %%%%.
vehicle = 'STOKES';

%%%% Root %%%%
ROOT = '../data/DELMAR2020/';

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
method = 'MEM';                                                             % Method used to compute directional spectrum 

%% Call Data

%--------- STOKES ---------% 
if isempty(whos('Ns')) || isempty(whos('Ws'))
    Ns = load([ROOT vehicle '/NOVATEL_downsampled_20Hz_ALL_' vehicle '.mat']); % RAW downsampled GPS/IMU Data (20 Hz)
    Ws = load([ROOT vehicle '/PLD2_1hz_ALL.mat']);                             % WXT/GILL Weather Station Data (1 hz)
    
    % Transpose time_20hz field: 
    Ns.nov.time_20hz = Ns.nov.time_20hz';

end

%% Compute spectrogram 

% Set variables
field_20hz = {'heave'; 'VEL_east'; 'VEL_north'; 'VEL_up'; 'time_20hz'};
field_1hz = {'TWS'; 'time'};

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
    It_ws = find(Ws.PLD2_1hz.time >= T0(is) & Ws.PLD2_1hz.time < T1(is));
    
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

    %%%%%%% 1 Hz data %%%%%%%
    % Loop through fieldnames
    for k = 1:length(field_1hz)

        % Set filenames
        name = char(field_1hz(k,:));

        % Index data
        eval(['w_s.' name ' = Ws.PLD2_1hz.' name '(:,It_ws);'])

    end
    
    % Compute the time elapsed from the first measurement in the time
    % interval for each time recorded in seconds
    nov_s.t = (nov_s.time_20hz - nov_s.time_20hz(1))*86400; time_int_ns = 0:dt_n:round(max(nov_s.t,[],'omitnan'));
    w_s.t = (w_s.time - w_s.time(1))*86400; time_int_ws = 0:dt_w:round(max(w_s.t,[],'omitnan'));
    
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

    %%%%%%% 1 Hz data %%%%%%%
    
    % Loop through field names (except time) 
    for k = 1:length(field_1hz)
        
        % Set field variable name
        name = char(field_1hz(k,:));
        
        % Find indices of non-nan data points
        eval(['Inum = find(~isnan(w_s.' name '));'])
        
        % Interpolate time steps with nans
        eval(['w_s.' name ' = interp1(w_s.t(Inum),w_s.' name '(Inum),time_int_ws,''linear'',''extrap'');'])
        
    end

    % Computing platform speed and direction 
    [nov_s.ground_speed_ms, nov_s.true_course] = platform_speed_dir(nov_s.VEL_east, nov_s.VEL_north, dir_con);
    
    % Set counter
    kr = kr + 1;
    
    % Compute directional spectrum  
    [nov_s.Sd(:,:,kr), nov_s.f, nov_s.theta] = compute_directional_spectrum(nov_s.heave, nov_s.VEL_east, nov_s.VEL_north, nov_s.VEL_up, nov_s.time_20hz, f, df, dtheta, nfft, fe_n, toolbox, variables, scaling, dir_con, method);

    % Compute the Omni-directional Spectra from directional wave spectra 
    nov_s.spectrogram_omni(:,kr) = sum(nov_s.Sd(:,:,kr) * dtheta, 1);

    % Compute mean wind speed over time interval and friction velocity 
    w_s.mTWS(kr) = mean(w_s.TWS);
    [w_s.mfv(kr),~,~] = FUNC_ustar_Z(z_station,w_s.mTWS(is));

    % Compute transition frequencies
    nov_s.fp(kr) = nov_s.f(nov_s.spectrogram_omni(:,is) == max(nov_s.spectrogram_omni(:,is), [], 'omitnan'));  % Peak frequency
    nov_s.feq(kr) = sqrt(2.25)*nov_s.fp(is);                                                                   % Spectral peak to equilibrium range frequency transition
    nov_s.fst(kr) = ((g*sqrt(r))/(2*pi*w_s.mfv(is)));                                                          % Equilibrium to saturation range frequency transition

    % Compute mean direction of platform
    [nov_s.mD(kr), ~] = direction_stats(nov_s.true_course, dt_n, 0);        % Convention: CW, going towards, ref north
    
    % Compute mean speed using instantaneous speed and projected
    % velocity
    nov_s.mspeed(kr) = mean(nov_s.ground_speed_ms);
    nov_s.ground_speed_ms_proj = nov_s.ground_speed_ms .* cosd(nov_s.true_course - nov_s.mD(kr));
    nov_s.mspeed_proj(kr) = mean(nov_s.ground_speed_ms_proj);
    
    % Compute relative angle between platform and wave direction
    nov_s.rel_theta(:,kr) = mod(nov_s.theta - nov_s.mD(kr), 360);           % Convention: CW, going towards, ref north

    % Compute the half way point of the time series
    T_leg = Ns.nov.time_20hz(It_s);
    T_half = T_leg(round(length(It_s) - (1/2)*length(It_s)));
    
    % Floor time step 
    nov_s.time_legs(is) = round(T_half*1440)/1440;
    
    % Display center interval 
    disp(datestr(nov_s.time_legs(is)))
    
end 

% Obtain indicies frequency below high frequency cut off
Inoise_s = find(nov_s.f <= f_noise);

% Remove frequencies and beyond 1 Hz 
nov_s.f_ob = nov_s.f(Inoise_s);
nov_s.spectrogram_omni_f_ob = nov_s.spectrogram_omni(Inoise_s,:);
nov_s.Sd_f_ob = nov_s.Sd(:,Inoise_s,:);

% Compute equivalent saturation spectrum 
nov_s.sat_spectrogram_omni_f_ob = nov_s.spectrogram_omni_f_ob .* (nov_s.f_ob').^(5);

%% Account for Doppler shifts by mapping observed frequency to intrinsic frequency
clc

% Set variables for mapping
tail = [zeros(size(nov_s.fst)); nov_s.feq; nov_s.fst]; 

%------ Map wave spectra ------%

% Loop through time increments
for it = 1:length(T0)
    
    % Map Directional Spectra 
    [nov_s.dir_spectrogram_f_in(:,:,it), nov_s.f_in(:,:,it), nov_s.fb(:,:,it)]  = map_dir_spectrum(nov_s.Sd_f_ob(:,:,it)', nov_s.f_ob, f_noise, nov_s.mspeed_proj(it), nov_s.rel_theta(:,it));

    % Compute Omni-directional Spectra from directional spectra 
    spectrogram_omni_f_in_nt = sum(nov_s.dir_spectrogram_f_in(:,:,it) * dtheta, 2);

    % Find non-nan values in S_fin after mapping
    idx_nan = ~isnan(spectrogram_omni_f_in_nt);

    % Attach a spectral tail if specified
    if tail(1,it) == true
        S_truc = spectrogram_omni_f_in_nt(idx_nan);
        f_truc = nov_s.f_ob(idx_nan);
        f_tail = [f_truc(end), nov_s.f_ob(~idx_nan)];
        [nov_s.spectrogram_omni_f_in(:,it), ~, fit, f_fit] = omnidir_spectral_tail(S_truc, f_truc, f_tail, tail(2,it), tail(3,it), f_noise);
    else
        nov_s.spectrogram_omni_f_in(:,it) = spectrogram_omni_f_in_nt;
    end

end

% Compute equivalent saturation spectrum for intrinsic frequency spectrum
nov_s.sat_spectrogram_omni_f_in = nov_s.spectrogram_omni_f_in .* (nov_s.f_ob').^(5);

% Save spectrograms to mat file
time = nov_s.time_legs; f = nov_s.f_ob; spec_ob = nov_s.sat_spectrogram_omni_f_ob; spec_in = nov_s.sat_spectrogram_omni_f_in;
save([ROOT vehicle '/DELMAR2020_spec_' vehicle '.mat'], 'time', 'f', 'spec_ob', 'spec_in')

% Save data for 2D method tutorial 
f_ob = nov_s.f_ob;  S_fob = permute(nov_s.Sd_f_ob,[2,1,3]); f_cut = f_noise; U_bar = nov_s.mspeed_proj; theta_r = nov_s.rel_theta;
save([ROOT vehicle '/2D_method_tutorial_data.mat'], 'time', 'f_ob', 'S_fob', 'f_cut', 'U_bar', 'theta_r')

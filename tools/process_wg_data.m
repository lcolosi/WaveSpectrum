function [nov, w, nlegs] = process_wg_data(vehicle, ROOT, date_o, dir_con, dt)
    
    %%%%
    % [nov, w, nlegs] = process_wg_data(vehicle, ROOT, date_o, dir_con, dt)
    %
    % Function processing wave glider novatel and meterological data for
    % the Del Mar experiment. Uploads novatel down sampled 20 Hz data,
    % payload 2 (1 Hz) meterological data, and big and small square
    % trajectory data and preforms the following operations: 
    %
    %   (1) Grabs the following data: 
    %           (1) Heave
    %           (2) VEL_east 
    %           (3) VEL_north
    %           (4) VEL_up
    %           (5) latitude
    %           (6) longitude
    %           (7) time_20hz
    %           (8) TWD
    %           (9) TWS
    %       from a time interval specified by date_o. 
    %   (2) Interpolates over gaps in the meterological data (Novatel data 
    %       does not have gaps as a result of the downsampling).
    %   (3) Computes ground speed and true course of wave glider from
    %       velocity components.
    %   (4) Segments data per leg that will be used for the computing wave 
    %       spectra and mapping it into intrinsic frequency space. 
    %
    %   Parameters
    %   ----------
    %   Vehicle : Specifies the vehicle to process data for. Keyword 
    %             argument is a stings. Options include: 'PLANCK', and
    %             'STOKES'.
    %   ROOT : Specifies the relative or absolute path to data. Example: 
    %          ROOT = '../data/DELMAR2020/';
    %   date_o : Specifies the date to index data. Format is a matrix with
    %            the elements  date_o = [date_i, date_f].
    %   dir_con : Specifies the directional convention for the platform's 
    %             progation direction. Format of keyword argument is: 
    %             dir_con = {'rotation direction', 
    %                        'coming from/going towards',
    %                        'zero angle reference direction'}
    %             Options include: 
    %                 (1) rotation direction = 'CW' (Clockwise), 'CCW' (counter clockwise) 
    %                 (2) coming from/going towards = 'cf', 'gt'
    %                 (3) zero angle reference direction = 'rn' (refernece
    %                     north), 're' (reference east) 
    %             Note, dir_con must be a cell array.
    %   dt : Specifies the time elapsed between measurements. Keyword
    %        argument is a two element matrix of the form: dt = [dt_n, dt_w]
    %        where dt_n is for novatel data ad dt_w is for weather station
    %        data. Recall dt = 1/fe where fe is the sampling frequency. 
    % 
    %   Returns
    %   -------
    %   nov : Novatel data structure. 
    %   w : Meterological data structure
    %   nlegs : Number of legs
    %
    %   Notes
    %   -----
    %   (1) Original directional conventions: 
    %        (i) nov.true_course : CCW, going towards, reference east
    %        (ii) w.TWD : CW, coming from, reference north
    %
    %%%%
    
    %% Call Data 
    
    %---------- Novatel downsampled Data (20 Hz) ----------%
    NOV = load([ROOT vehicle '/NOVATEL_downsampled_20Hz_ALL_' vehicle '.mat']); 
    
    %---------- Weather Station Data (1 hz) ----------%
    W = load([ROOT vehicle '/PLD2_1hz_ALL.mat']); 
    
    %---------- Big Square Leg Coordinates ----------%
    LEGB = load([ROOT vehicle '/LEGS_' vehicle '_LARGE_SQUARE.mat']); 
    
    %---------- Small Square Leg Coordinates ----------%
    LEGS = load([ROOT vehicle '/LEGS_' vehicle '_SMALL_SQUARE.mat']); 

    % Transpose time_20hz field: 
    NOV.nov.time_20hz = NOV.nov.time_20hz';
    
    %% Index data for specified time period and interpolating over NaNs for Weather data
    
    % Set variables
    field_20hz = {'heave'; 'VEL_east'; 'VEL_north'; 'VEL_up'; 'latitude'; 'longitude';'time_20hz'; 'true_course'; 'ground_speed_ms'};
    field_1hz = {'TWD'; 'TWS'; 'time'};

    % Find indices for desired time frame
    time_20hz_idx = find(NOV.nov.time_20hz >= datenum(date_o(1)) & NOV.nov.time_20hz < datenum(date_o(2)));
    time_1hz_idx = find(W.PLD2_1hz.time >= datenum(date_o(1)) & W.PLD2_1hz.time < datenum(date_o(2)));
    
    % Create the structures consisting of data from specified time frame

    %%%%%%% 20 Hz data %%%%%%%

    % Loop through fieldnames
    for k = 1:length(field_20hz) - 2

        % Set filenames
        name = char(field_20hz(k,:));

        % Index data
        eval(['nov.' name ' = NOV.nov.' name '(time_20hz_idx);'])

    end
    
    %%%%%%% 1 Hz data %%%%%%%

    % Loop through fieldnames
    for k = 1:length(field_1hz)

        % Set filenames
        name = char(field_1hz(k,:));

        % Index data
        eval(['w.' name ' = W.PLD2_1hz.' name '(time_1hz_idx);'])

    end

    % Compute the time elapsed from the first measurement in the time
    % interval for each time recorded in seconds
    nov.t = (nov.time_20hz - nov.time_20hz(1))*86400; time_int_n = 0:dt(1):round(max(nov.t,[],'omitnan'));
    w.t = (w.time - w.time(1))*86400; time_int_w = 0:dt(2):round(max(w.t,[],'omitnan'));
    
    %%%%%%% 20 Hz data %%%%%%%
    
    % Loop through field names
    for k = 1:length(field_20hz)- 2
        
        % Set field variable name
        name = char(field_20hz(k,:));
        
        % Find indices of non-nan data points
        eval(['Inum = find(~isnan(nov.' name '));'])
        
        % Interpolate time steps with nans
        eval(['nov.' name ' = interp1(nov.t(Inum),nov.' name '(Inum),time_int_n,''linear'',''extrap'');'])
        
    end 
    
    %%%%%%% 1 Hz data %%%%%%%
    
    % Loop through field names (except time) 
    for k = 1:length(field_1hz)
        
        % Set field variable name
        name = char(field_1hz(k,:));
        
        % Find indices of non-nan data points
        eval(['Inum = find(~isnan(w.' name '));'])
        
        % Interpolate time steps with nans
        eval(['w.' name ' = interp1(w.t(Inum),w.' name '(Inum),time_int_w,''linear'',''extrap'');'])
        
    end


    %% Computing platform speed and direction of progagtion 

    % Compute speeed
    nov.ground_speed_ms = sqrt(nov.VEL_east.^2 + nov.VEL_north.^2); 

    % Compute direction and align directional convention for wind and
    % platform

    %---------- CW, going towards, ref north ----------%
    if strcmp(dir_con(1), 'CW') && strcmp(dir_con(2), 'gt') && strcmp(dir_con(3), 'rn')

        nov.true_course = mod(90 - atan2d(nov.VEL_north,nov.VEL_east), 360);  
        w.TWD = mod(180 + w.TWD, 360);
        
    %---------- CW, coming from, ref north ----------%
    elseif strcmp(dir_con(1), 'CW') && strcmp(dir_con(2), 'cf') && strcmp(dir_con(3), 'rn')

        nov.true_course = mod(270 - atan2d(nov.VEL_north,nov.VEL_east), 360);
        
    end

    
    %% Segment data into trajectory legs  
    
    % Concatinate start and end times for all legs
    LEGs = [LEGB.LEG_000.time_start_leg, LEGB.LEG_090.time_start_leg,...    % Big square start times
            LEGB.LEG_180.time_start_leg, LEGB.LEG_270.time_start_leg,...
            LEGS.LEG_000.time_start_leg, LEGS.LEG_090.time_start_leg,...    % Small square start times end times
            LEGS.LEG_180.time_start_leg, LEGS.LEG_270.time_start_leg];
    LEGe = [LEGB.LEG_000.time_end_leg, LEGB.LEG_090.time_end_leg,...        % Big square end times
            LEGB.LEG_180.time_end_leg, LEGB.LEG_270.time_end_leg,... 
            LEGS.LEG_000.time_end_leg, LEGS.LEG_090.time_end_leg,...        % Small square end times
            LEGS.LEG_180.time_end_leg, LEGS.LEG_270.time_end_leg];

    % Set number of legs
    nlegs = length(LEGs);
    
    % Sort start times
    [LEGs_s, Idx_legs] = sort(LEGs);

    %Sort end times based on start times
    LEGe_s = LEGe(Idx_legs);
    
    % Loop through time start times 
    for ileg = 1:length(LEGs_s)

        % Find indices for leg
        Idx_20hz = find(nov.time_20hz >= LEGs_s(ileg) & nov.time_20hz <= LEGe_s(ileg));
        Idx_1hz = find(w.time >= LEGs_s(ileg) & w.time <= LEGe_s(ileg));

        % Loop through 20hz field variables
        for k = 1:length(field_20hz)

            % Set filenames
            name = char(field_20hz(k,:));

            % Index data
            eval(['nov.L' num2str(ileg) '.' name ' = nov.' name '(Idx_20hz);'])

        end

        % Loop through 1hz field variables
        for k = 1:length(field_1hz)

            % Set filenames
            name = char(field_1hz(k,:));

            % Index data 
            eval(['w.L' num2str(ileg) '.' name ' = w.' name '(Idx_1hz);'])

        end
    end
    
end

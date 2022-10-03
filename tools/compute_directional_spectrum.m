function [S, freq, theta] = compute_directional_spectrum(heave, vel_east, vel_north, vel_up, time, f, df, dtheta, nfft, fe, toolbox, variables, scaling, dir_con)

    %%%%
    % [S, freq, theta] = compute_directional_spectrum(heave, vel_east, vel_north, vel_up, time, f, df, dtheta, nfft, fe, toolbox, variables, scaling, dir_con)
    %
    % Function for computing the directional wave spectrum from the motion 
    % of an autonomous platform (vertical dispacement (heave), and
    % vertical and horizontal velocity ). Directional Spectrum is 
    % calculated using one of two toolboxes: CDIP toolbox (maximum entropy
    % method only) or WAFO toolbox. The maximum entropy method is used for
    % spectrum estimate. Heave and velocity data must be properly
    % interoplated (NO NAN VALUES!) and filtered (high pass or low pass 
    % filter depending your data; heave and velocity components time
    % series are detrended within function) for spectral analysis. 
    %
    %   Parameters
    %   ----------
    %   heave : Verticle displacement of the instrument. Quantity must only
    %           have real values. Time series must be detrended (mean and 
    %           linear trend removed from time series). However, if not the
    %           heave is detrended within function. Units: m.
    %   vel_east : Eastward component of velocity. Quantity must only
    %              have real values. Units: ms^-1.
    %   vel_north : Northward component of velocity. Quantity must only
    %               have real value. Units: ms^-1.
    %   vel_up : Vertical component of velocity. Quantity must only
    %            have real values. Units: ms^-1.
    %   time : Time vector. Units : s.
    %   f : Cyclical Frequency where f = 0:df:fn where df is the frequency
    %       resolution and fn is the nyquist frequency where fn = 1/2*fe. 
    %       Units: s^-1.
    %   df : Frequency resolution. Units: s^-1.
    %   nfft : Window Length for Welch's power spectrum estimate.
    %   fe : Sampling rate of sensor. Units: s^-1.
    %   toolbox : Specifies the toolbox used to compute directional spectra. 
    %             Options include: 'CDIP' or 'WAFO'.
    %   variables : Specifies the variables used to compute the
    %               direction spectrum for the WAFO toolbox. Options
    %               include: 'velocity' and 'heave_velocity'. 'velocity'
    %               argument means the three components of velocity are
    %               used. 'heave_velocity' means heave and the
    %               horizontal components of velocity are used. If the CDIP
    %               toolbox is used set variables = []; 
    %   scaling : Specifies whether the variance of the directional wave
    %             spectrum will be scaled to match the variance of the
    %             heave spectrum. If true, then directional wave spectrum
    %             will be scaled; if false, it will not. 
    %   dir_con : Specifies the directional convention for wave directions.
    %             Format of keyword argument is:
    %             dir_con = {'rotation direction', 'coming from/going towards',
    %                        'zero angle reference direction'}
    %             Options include: 
    %                 (1) rotation direction = 'CW' (Clockwise), 'CCW' (counter clockwise) 
    %                 (2) coming from/going towards = 'cf', 'gt'
    %                 (3) zero angle reference direction = 'rn' (reference
    %                     north), 're' (reference east)
    %             Note, dir_con must be a cell array.
    %   
    %   Returns
    %   -------
    %   S : Directional Spectrum. Units: m^2 (Hz deg)^(-1).
    %   freq : Cyclical frequency. Units: Hz.
    %   theta : Azimutal direction for wave spectrum. Units: degrees 
    %           (directional convention set by dir_con).
    % 
    %   Notes
    %   -----
    %   (1) Original directional convention of theta before converting
    %       convention to dir_con: 
    %
    %           CDIP -> [CCW, going towards, reference East]
    %           WAFO -> [CCW, going towards, reference East] when bet = 1 (default). 
    % 
    %   (2) The WAFO toolbox can be found for download at the wafo website:
    %        
    %           https://www.maths.lth.se/matstat/wafo/download/index.html
    % 
    %       Further inforation may be found on thier github page: 
    % 
    %           https://github.com/wafo-project/wafo
    %       
    %
    %%%%
    
    % Remove mean and linear trend from heave time series 
    heave = detrend(heave); 
    vel_east = detrend(vel_east);
    vel_north = detrend(vel_north);
    vel_up = detrend(vel_up);
    
    %%%%%%%%% Toolbox 1: CDIP %%%%%%%%%
    if strcmp(toolbox, 'CDIP')
    
        % Compute the power density spectrum for the vertical position,
        % horizontal velocity, and vertical velocity with 50% overlap:
        [Pz,~]  = pwelch(heave,nfft,nfft/2,nfft,fe);
        [Pe,~]  = pwelch(vel_east,nfft,nfft/2,nfft,fe);
        [Pn,~]  = pwelch(vel_north,nfft,nfft/2,nfft,fe);
        [Pu,~]  = pwelch(vel_up,nfft,nfft/2,nfft,fe);

        % Compute the cross spectrum for the vertical position,
        % horizontal velocity, and vertical velocity with 50% overlap:
        [Peu,~] = cpsd(vel_up,vel_east,nfft,nfft/2,nfft,fe);
        [Pnu,~] = cpsd(vel_up,vel_north,nfft,nfft/2,nfft,fe);
        [Pen,~] = cpsd(vel_north,vel_up,nfft,nfft/2,nfft,fe);

        % Rename the Power densities and set co- and quadrature- spectra
        Szz = Pz;
        See = Pe;
        Snn = Pn;
        Suu = Pu;
        Ceu = real(Peu);
        Qeu = imag(Peu);
        Cnu = real(Pnu);
        Qnu = imag(Pnu);
        Cen = real(Pen);
        Qen = imag(Pen);

        % Compute first- and second-order Fourier coefficients of the directional 
        % distribution of wave energy
        a1 = Qeu ./ sqrt( (See + Snn).*Suu );
        b1 = Qnu ./ sqrt( (See + Snn).*Suu );
        a2 = (See - Snn) ./ (See + Snn) ;
        b2 =  2*Cen      ./ (See + Snn);

        %%%%%% Compute Directional Spectra %%%%%%

        % Create in structure with field variables for mem_est function
        in.freq = f;  % Frequency 
        in.dir = ones(size(in.freq));                                       % Directions 
        in.ener_dens = Szz*df;                                              % Energy  
        in.a1 = a1;                                                         % First-order fourier coefficient
        in.b1 = b1;                                                         % First-order fourier coefficient
        in.a2 = a2;                                                         % Second-order fourier coefficient
        in.b2 = b2;                                                         % Second-order fourier coefficient

        % Compute directional spectrum 
        [cdip,~] = mem_est(in);

        % Rename fields in cdip structure
        freq = cdip.freq;                                                   % Frequency
        theta = cdip.dir';                                                  % Directions
        S = cdip.ds;                                                        % Directional Wave Spectrum 

        % Set directional convention for 
        %---- Clockwise/Counter-Clockwise ----%
        if strcmp(dir_con(1), 'CW')
            theta = -theta;
        end
    
        %---- Going-to/Coming-from ----%
        if strcmp(dir_con(2), 'cf')
            theta = theta + 180; 
        end

        %---- Reference North/East ----%
        if strcmp(dir_con(3), 'rn')
            theta = theta + 90;
        end
        
        % Make theta range from 0 to 360 
        theta = mod(theta,360);

        % Sort directions and the directional index of the Directional Wave Spectrum
        % in ascending order 0 to 360 degrees
        [theta,Is] = sort(theta);
        S = S(Is,:); 

        % Convert the Directional Wave Spectrum to units of m^2/Hz/deg (scale with E)
        if scaling == true
            S = S * sum(Szz*df) / sum(S(:)*freq(2)*theta(2)); 
        end
    
    %%%%%%%%% Toolbox 2: Wave Analysis for Fatigue an Oceanography (WAFO) %%%%%%%%%
    elseif strcmp(toolbox, 'WAFO')
    
        % Initialize variables
        method_comp = 'MEM';
        bet = 1;                                                            % Theta given in terms of directions toward which waves travel
        ntheta = round(360/dtheta)+1;
        xyz = zeros(3,3);
        bfs = [ 0; 0; 1];
        % Case 1: heave and horizontal velocity
        if strcmp(variables, 'heave_velocity') 
            def = [10;11;1];
            W = [time(:) vel_east(:) vel_north(:) heave(:)]; 
        % Case 2: Horizontal and vertical velocity
        elseif strcmp(variables, 'velocity') 
            def = [10;11;12];
            W = [time(:) vel_east(:) vel_north(:) vel_up(:)];
        end
        POS = [xyz def bfs];

        % Compute Directional Wave Spectrum
        [wafo,~,~,~]  = dat2dspec(W,POS,6000,nfft,ntheta,method_comp,'bet',bet);
        
        % Compute the power density spectrum for the vertical position with
        % 50% overlap:
        [E,~]  = pwelch(heave,nfft,nfft/2,nfft,fe);

        % Rename frequency, directions, and directional spectrum variables
        freq = f;                                                           % Units: Hz
        theta = wafo.theta * 180 / pi;                                      % Units: degrees
        S = wafo.S;                                                         % Units: m^2/(Hz deg)

        % Set directional convention for 
        %---- Clockwise/Counter-Clockwise ----%
        if strcmp(dir_con(1), 'CW')
            theta = -theta;
        end
    
        %---- Going-to/Coming-from ----%
        if strcmp(dir_con(2), 'cf')
            theta = theta + 180; 
        end

        %---- Reference North/East ----%
        if strcmp(dir_con(3), 'rn')
            theta = theta + 90;
        end
        
        % Make theta range from 0 to 360 
        theta = mod(theta,360);

        % Sort directions in ascending order 0 to 360 degrees and grabs unique
        % directions. 
        [theta,Isort] = sortrows(theta);
        [theta,Iunique] = unique(theta);

        % Convert the Directional Wave Spectrum to units of m^2/Hz/deg
        S = S(Isort(Iunique),:) / (2*pi);                                   % PSD in m^2/Hz/deg.
        
        % Scale ditrectional spectrum with variance from the heave power spectrum.
        if scaling == true
            S = S * sum(E*df) / sum(sum(S*df,2)*dtheta);                    % scale with rectangular integration of the heave power spectrum 
        end
        
    end
end
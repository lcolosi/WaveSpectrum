function [S_fin, f_int, f_b, variance] = map_omni_dir_spectrum(S_fob, f_obs, f_cut, U, theta_r, tail)
    
    %%%%
    % [S_fin, f_int, f_b, variance] = map_omni_dir_spectrum(S_fob, f_obs, f_cut, U, theta_r, tail)
    %
    % Function for mapping observed frequency (frequency measured in the 
    % moving reference frame of the platform) to intrinsic frequency 
    % (frequency measured in the absence of platform motion) for deep water
    % surface gravity waves. Function also maps the observed 1D spectrum 
    % S_ob(f_ob) (i.e., the omni-direction frequency spectrum) into 
    % intrinsic frequency space as the intrinsic 1D spectrum S_in(f_in). 
    % The observed frequency spectrum may be computed by either integrating
    % the observed 2D spectrum over all directions or computing the 
    % frequency spectrum of the platform s vertical displacement using the
    % Welch method   
    %
    %   Parameters
    %   ----------
    %   S_fob : Observed 1D spectrum S(f_obs). Units: m^2 Hz^(-1).
    %   f_obs : Observed Frequency. Units: Hz. 
    %   f_cut : Noise cutoff frequency for the platform. For the SV3 wave
    %           gliders, f_cut is nominally ~1 Hz (corresponds to waves 
    %           with a wavelength half the length of the waveglider).
    %           Units Hz. 
    %   U : Propagation speed of the observer.  Units: ms^(-1). 
    %   theta_r : Angle between the propagation velocity of the platform
    %             and the direction of propagation of waves. Here, this is
    %             scalar value such that the wave field is assumed to be 
    %             uni-directional. The wave direction may assumed to be 
    %             aligned with the mean wind direction. However, other wave
    %             direction may be used. Units: Degrees.
    %   tail : Specifies whether a spectral tail with a given slope is 
    %          added to the intrinsic frequency spectrum after the 
    %          blocking frequency (if the f_b < f_cut). Options include:
    %
    %                   tail = [true, f_eq, f_sat] or [false, f_eq, f_sat] 
    %
    %          where true or false determines if the tail will be added or
    %          not. f_eq is the transition from the spectral peak to the 
    %          equilibrium range f_sat os the transition between the 
    %          equilibrium and saturation ranges. Nomially, the equilibrium
    %          range has a spectral slope of f^{-4} and the saturation 
    %          range has a spectral slope of f^{-5}. Frequency units: Hz. 
    % 
    %   Returns
    %   -------
    %   S_fin : Intrinsic 1D spectrum. Units: m^2 Hz^(-1).
    %   f_int : Intrinsic Frequency. Units: Hz.
    %   f_b : Bifurcation frequency and last frequency not an NaN (blocking 
    %         frequency). The bifircation frequency exists for branches 3-5
    %         (platform moving in direction of wave propagation) and is NaN
    %         for branches 1-2. The bocking frequency is set by the
    %         jacobian, the linear interpolation and the bifurcation
    %         frequency. The Bifurcation frequency is saved in first column
    %         f_bifurcation = f_b(:,1) and blocking frequency is saved in
    %         second column f_block = f_b(:,2). Units: Hz.
    %   variance : Matrix containing four values for the variance computed
    %              over a given frequency band using rectangular (Riemann) 
    %              integration method: 
    %                   variance = [f_ob_0.01, f_in_0.01, f_ob_0.1, f_in_0.1]
    %             where: 
    %
    %             (1) f_ob_0.01 : Variance computed from integrating 
    %                             S_ob(f_ob)*df_ob from 0.01 to 1 Hz
    %             (2) f_in_0.01 : Total variance computed from integrating
    %                             S_in(f_in)*df_in from
    %                             f_in(f_ob = 0.01) to f_in(f_ob = 1) Hz
    %                             (head-seas case) and 1 Hz for
    %                             (following-seas case)
    %             (3) f_ob_0.1 : Variance computed from integrating 
    %                            S_ob(f_ob)*df_ob from 0.1 to 1 Hz
    %             (4) f_in_0.1 : Variance computed from integrating
    %                            S_in(f_in)*df_in from f_in(f_ob = 0.1) to 
    %                            f_in(f_ob = 1) Hz (head-seas case) and 
    %                            1 Hz for (following-seas case) 
    %
    %   Notes
    %   -----
    %   (1) Update variance calculation to include correct limits of
    %       integration. 
    %%%%
    
    % Set constants
    g = 9.81;                                                               % Gravitational Acceleration (units: ms^(-2))  
    df_obs = unique(round(diff(f_obs),4));                                  % Observed frequency resolution of input spectrum (units: Hz)
    f_obs_li = 0.01; f_obs_lf = 1;                                          % Limits of integration for the lower frequency domain integration (units: Hz)
    f_obs_hi = 0.1; f_obs_hf = 1;                                           % Limits of integration for the higher frequency domain integration (units: Hz)
    
    % Initialize variables
    S_fin = NaN(size(S_fob));
    f_b = NaN(1,2);
    
    % Compute the Doppler shift speed projected onto the ith wave direction
    c_p = U*cosd(theta_r); 

    %------- Branch Conditional statements -------%

    %---------- Branch 1 : Moving against waves ----------%
    if cosd(theta_r) < 0 && U > 0
        %disp('Branch 1 : Moving against waves');

        % Compute intrinsic frequency for branch 1
        f_int = (g - sqrt(g^2 - 8*pi*f_obs*g*c_p))/(4*pi*c_p);

        % Compute intrinsic frequency limits of integration
        f_int_li = (g - sqrt(g^2 - 8*pi*f_obs_li*g*c_p))/(4*pi*c_p);
        f_int_lf = (g - sqrt(g^2 - 8*pi*f_obs_lf*g*c_p))/(4*pi*c_p);
        f_int_hi = (g - sqrt(g^2 - 8*pi*f_obs_hi*g*c_p))/(4*pi*c_p);
        f_int_hf = (g - sqrt(g^2 - 8*pi*f_obs_hi*g*c_p))/(4*pi*c_p);

    %---------- Branch 2 : Moving noraml to waves ----------%
    elseif cosd(theta_r) == 0 || U == 0
        %disp('Branch 2 : Moving normal to waves ');
        
        % Compute intrinsic frequency for branch 2 
        f_int = f_obs; 

        % Compute intrinsic frequency limits of integration
        f_int_li = f_obs_li;
        f_int_lf = f_obs_lf;
        f_int_hi = f_obs_hi;
        f_int_hf = f_obs_hf;
        
    %---------- Branch 3-5 : Moving with wave ----------%
    elseif cosd(theta_r) > 0 && U > 0
        %disp('Branch 3-5 : Moving with wave');

        % Compute f_obs where df_in(f_obs)/df_ob tends towards infinity  
        f_ob_max = (g)/(8*pi*c_p);

        % Compute f_int at the f_obs value where df_in(f_obs)/df_ob
        % tends towards infinity and save in f_b array
        f_in_max = (g)/(4*pi*c_p);
        f_b(1,1) = f_in_max; 
        
        % Compute f_int when f_obs = 0
        f_in_zero = (g)/(2*pi*c_p);

        %---------- Platform moving slower than energy and crests ----------%
        if f_in_max > f_cut
            % disp('moving slower than energy and crests');
            
            % Compute the f_in value where f_int = f_cut for branch 3
            f_ob_cut = -((2*pi*c_p*f_cut^2)/(g)) + f_cut;

            % Obtain f_obs values that map to f_int less than or equal to f_ob_cut
            f_obs_map = f_obs(f_obs <= f_ob_cut);

            % Compute intrinsic frequency 
            f_int_s = (g - sqrt(g^2 - 8*pi*f_obs_map*g*c_p))/(4*pi*c_p);

            % Insert NaNs to extend the length of f_int_s to the length of
            % f_obs
            nNaNs = length(f_obs) - length(f_int_s); 
            f_int = [f_int_s, nan(nNaNs, 1)'];

            % Compute intrinsic frequency limits of integration
            f_int_li = (g - sqrt(g^2 - 8*pi*f_obs_li*g*c_p))/(4*pi*c_p);
            f_int_lf = f_obs_lf;
            f_int_hi = (g - sqrt(g^2 - 8*pi*f_obs_hi*g*c_p))/(4*pi*c_p);
            f_int_hf = f_obs_hf;

        %---------- Platform moving faster than energy but slower than crests ----------%
        elseif f_in_max  > 0 && f_in_max < f_cut && f_in_zero > f_cut
            % disp('moving faster than energy but slower than crests');
            
            % Obtain f_obs values for mapping to branch 3 (f_obs less than or equal to f_ob_max)
            f_obs_map = f_obs(f_obs <= f_ob_max);

            % Compute intrinsic frequency for branch 3
            f_int_s = (g - sqrt(g^2 - 8*pi*f_obs_map*g*c_p))/(4*pi*c_p);
            
            % Insert NaNs to extend the length of f_int to the length of
            % f_in
            nNaNs = length(f_obs) - length(f_int_s); 
            f_int = [f_int_s, nan(nNaNs, 1)'];

            % Compute intrinsic frequency limits of integration
            f_int_li = (g - sqrt(g^2 - 8*pi*f_obs_li*g*c_p))/(4*pi*c_p);
            f_int_lf = f_obs_lf;
            f_int_hi = (g - sqrt(g^2 - 8*pi*f_obs_hi*g*c_p))/(4*pi*c_p);
            f_int_hf = f_obs_hf;

        %---------- Platform moving faster than energy and crests ----------%
        elseif f_in_zero < f_cut
            % disp('moving faster than energy and crests');
            
            % Obtain f_obs values for mapping to branch 3 (f_obs less than or equal to f_ob_max)
            f_obs_map = f_obs(f_obs <= f_ob_max);

            % Compute intrinsic frequency for branch 3
            f_int_s = (g - sqrt(g^2 - 8*pi*f_obs_map*g*c_p))/(4*pi*c_p);
            
            % Insert NaNs to extend the length of f_in to the length of
            % f_in
            nNaNs = length(f_obs) - length(f_int_s); 
            f_int = [f_int_s, nan(nNaNs, 1)'];

            % Compute intrinsic frequency limits of integration
            f_int_li = f_obs_li;
            f_int_lf = f_obs_lf;
            f_int_hi = f_obs_hi;
            f_int_hf = f_obs_hf;
            
        end
    end 
    
    %-- Compute Intrinsic 1D Spectrum S(f_in) --%

    % Compute f_int differentials via finite differencing
    df_int = diff(f_int);

    % Compute Jacobian 
    J = df_obs./df_int;

    % Compute S_fin by multiplying S_fob by Jacobian 
    S_fin(1:length(J)) = S_fob(1:length(J)) .* J';

    % Map S_fin back onto the f_in regular grid (linear interpolation is not extrapolating to values where f_obs > f_int(end))
    idx_nan = ~isnan(S_fin); 
    S_fin = interp1(f_int(idx_nan),S_fin(idx_nan),f_obs)';

    % Find non-nan values in S_fin after mapping (mapping between f_int and
    % f_obs introduces new nans because ranges are incompatible)
    idx_nan = ~isnan(S_fin); 

    % Set the non-nan values for the intrinsic frequency and power
    % spectrum to truc variables
    f_truc = f_obs(idx_nan);
    S_truc = S_fin(idx_nan);

    % Set blocking frequency
    f_b(1,2) = f_truc(end);

    % Set frequency at which the tail will be attached at
    f_a = f_truc(end);
    
    % Obtain frequencies of tail
    f_tail = [f_a, f_obs(~idx_nan)];

    % Add a spectral tail if specified
    if tail(1) == true
        [S_fin, ~, ~, ~] = omnidir_spectral_tail(S_truc, f_truc, f_tail, tail(2), tail(3), f_cut);
    end

    % Obtain indices for integration (non-NaN values at frequencies higher the 0.01 and 0.1 Hz)
    idx_fob_l = ~isnan(S_fob) & (f_obs >= f_obs_li)' & (f_obs <= f_obs_lf)'; 
    idx_fin_l = ~isnan(S_fin) & (f_int >= f_int_li)' & (f_int <= f_int_lf)';
    idx_fob_h = ~isnan(S_fob) & (f_obs >= f_obs_hi)' & (f_obs <= f_obs_hf)'; 
    idx_fin_h = ~isnan(S_fin) & (f_int >= f_int_hi)' & (f_int <= f_int_hf)';

    % Compute variance by interagting over two bands of the spectrum
    variance = [sum(S_fob(idx_fob_l)*df_obs),...
                sum(S_fin(idx_fin_l)*df_obs),...
                sum(S_fob(idx_fob_h)*df_obs),...
                sum(S_fin(idx_fin_h)*df_obs)];
        
end
function [S_fin, f_int, f_b] = map_dir_spectrum(S_fob, f_obs, f_cut, U, theta_r)
    
    %%%%
    % [S_fin, f_int, f_b] = map_dir_spectrum(S_fob, f_obs, f_cut, U, theta_r)
    %
    % Function for mapping observed frequency (frequency measured in the 
    % moving reference frame of the platform) to intrinsic frequency 
    % (frequency measured in the absence of platform motion) for deep water
    % surface gravity waves. The function also maps the observed 2D 
    % spectrum S_ob(f_ob, theta) (i.e., the directional-frequency spectrum)
    % into intrinsic frequency space thus computing intrinsic 2D spectrum 
    % S_in(f_in, theta). The observed 2D spectrum may be estimated from 
    % vertical displacement of a platform and horizontal cartesian 
    % components of velocity using a maximum extropy or maximum likelyhood
    % method. However, other method exist. 
    %
    %   Parameters
    %   ----------
    %   S_fob : Observed 2D spectrum. Units: m^2 (Hz deg)^(-1).
    %   f_obs : Observed frequency.  Units: Hz. 
    %   f_cut : Noise cutoff frequency for the platform. For the SV3 wave
    %           gliders, f_cut is nominally ~1 Hz (corresponds to waves 
    %           with a wavelength half the length of the waveglider). Other
    %           more conservative cutoff frequencies may be used. This is 
    %           dependent on the frequency response and the sampling rate
    %           of the platform. Units Hz. 
    %   U : Propagation speed of the observer. Computed by time averaging 
    %       the speed of the platform projected onto the mean platform
    %       heading. Units: ms^(-1). 
    %   theta_r : Angle between the propagation velocity of the platform 
    %             and the direction of propagation of waves. This quantity 
    %             is a evaluated for all wave directions.
    %             Units: Degrees, CW, going towards, reference north.
    % 
    %   Returns
    %   -------
    %   S_fin : Intrinsic 2D spectrum. Units: m^2 (Hz deg)^(-1).
    %   f_int : Intrinsic frequency (computed for each direction). 
    %           Units: Hz.
    %   f_b : Bifurcation frequency and last frequency not an NaN (f_b 
    %         frequency). The bifircation frequency exists for branches 3-5
    %         (platform moving in direction of wave propagation) and is NaN
    %         for branches 1-2. The f_b frequency is set by the
    %         jacobian, the linear interpolation and the bifurcation
    %         frequency. The Bifurcation frequency is saved in first column
    %         f_bifurcation = f_b(:,1) and f_b frequency is saved in
    %         second column f_b = f_b(:,2). Units: Hz. 
    % 
    %   Notes
    %   -----  
    %   (1) ** A note on mapping to other reference frames ** 
    %       Mapping into another reference frame (other than the frame 
    %       absent of platform motion as suggested above) may be done by
    %       changing the variables U and theta_r. These quantities set the
    %       Doppler shift velocity. I recommend the researcher to compute 
    %       the Doppler shift velocity of their choice (nominally the 
    %       relative speed between frames) and then compute speed and
    %       direction from this velocity vector. This allows you to apply
    %       other mapping approaches such as from Hanson et al. 1997 
    %       which maps into the reference frame moving with the mean 
    %       currents.    
    % 
    %       IMPORTANT: Make sure to take the difference between the
    %       direction of the Doppler shift velocity and each wave
    %       direction.
    % 
    %   (1) ** Decision about adding a spectral tail **  
    %       Adding a spectral tail to the directional spectrum requires
    %       information about the directional dependence of the spectral
    %       slope in the equilibirum and saturation ranges. The
    %       conventional spectral forms of the spectra (i.e., f^-4 and 
    %       f^-5 from Phillips 1985 spectral model) do not include this 
    %       information. Therefore, the functional form
    %       must be either updated to include directional dependence or the
    %       spectral tail attachment must be preformed on the
    %       omni-directional spectrum. I choose the latter.
    % 
    %   (2) ** Discussion about conservation of variance ** 
    %       Comparing the total variance in the intrinsic and observed
    %       reference frames requires integrating the observed and
    %       intrinsic 2D spectrum over equivalent frequency bands given
    %       that the integrands S_ob and S_in = S_ob*df_ob/df_in are 
    %       equivalent (see equation 15 in main text). We use the mapping
    %       function (see equation 7 in main text) to obtain the limits of 
    %       integration after the change of variable to intrinsic 
    %       frequency.  In the head-sea and perpendicular-to-wave cases, 
    %       the mapping is bijective such that the limits of integration
    %       are well defined and the total variance before and after may be
    %       compared. However, in the following seas case, the mapping is
    %       no longer bijective. This is problematic because we loss 
    %       information about the high frequency wave above the bifurcation 
    %       frequency and the variance associated with these waves. 
    %       Therefore, the total variance is intrinsically not conserved
    %       for the following-seas mapping. We could constrain our
    %       integral to a low frequency band where the mapping is
    %       bijective. Note, this bijective mapping region is dependent
    %       on the speed and relative direction of the platform.
    %       I choose to not compute the total variance. In the future, I
    %       could implement code that works backwards from the integral of
    %       the intrinsic 2D spectrum setting the limits of intergation to
    %       the bijective region and mapping these limits back into
    %       observed frequency.
    % 
    %   Contact Information
    %   -------------------
    %   If you run into any diffculties using this code, please reach out
    %   to me via email at lcolosi@ucsd.edu. I'd be happy to help in any
    %   way I can! 
    %
    %%%%
    
    % Set constants
    g = 9.81;                                                               % Gravitational Acceleration (units: ms^(-2))
    df_obs = unique(round(diff(f_obs),4));                                  % Observed frequency resolution of input spectrum (units: Hz)
    
    % Initialize variables 
    S_fin = NaN(size(S_fob));
    f_int = zeros(size(S_fob)); 
    f_b = NaN(length(theta_r),2);

    %------- Map from observed to intrinsic frequency -------%
    
    % Loop through relative angle 
    for iangle = 1:length(theta_r)
        
        % Call ith relative angle
        itheta_r = theta_r(iangle);

        % Compute the Doppler shift speed projected onto the ith wave
        % direction
        c_p = U*cosd(itheta_r); 

        %---------- Branch 1 : Moving against waves ----------%
        if cosd(itheta_r) < 0 && U > 0
            %disp('Branch 1 : Moving against waves');

            % Compute intrinsic frequency 
            f_int(:,iangle) = (g - sqrt(g^2 - 8*pi*f_obs*g*c_p))/(4*pi*c_p);
            
        %---------- Branch 2 : Moving noraml to waves ----------%
        elseif cosd(itheta_r) == 0 || U == 0
            %disp('Branch 2 : Moving normal to waves');

            % Compute intrinsic frequency for branch 2 
            f_int(:,iangle) = f_obs;

        %---------- Branch 3-5 : Moving with wave ----------%
        elseif cosd(itheta_r) > 0 && U > 0
            %disp('Branch 3-5 : Moving with wave');

            % Compute f_obs where df_in(f_obs)/df_ob tends towards infinity  
            f_ob_max = (g)/(8*pi*c_p);

            % Compute f_int at the f_obs value where df_in(f_obs)/df_ob
            % tends towards infinity and save in f_b array
            f_in_max = (g)/(4*pi*c_p);
            f_b(iangle,1) = f_in_max; 

            % Compute f_int when f_obs = 0
            f_in_zero = (g)/(2*pi*c_p);

            %---------- Platform moving slower than energy and crests ----------%
            if f_in_max > f_cut
                % disp('moving slower than energy and crests');

                % Compute the f_ob value where f_int = f_cut for branch 3
                f_ob_cut = -((2*pi*c_p*f_cut^2)/(g)) + f_cut;

                % Obtain f_obs values that map to f_int less than or equal to f_ob_cut
                f_obs_map = f_obs(f_obs <= f_ob_cut);

                % Compute intrinsic frequency 
                f_int_s = (g - sqrt(g^2 - 8*pi*f_obs_map*g*c_p))/(4*pi*c_p);
                
                % Insert NaNs to extend the length of f_int_s to the length of
                % f_obs
                nNaNs = length(f_obs) - length(f_int_s); 
                f_int(:,iangle) = [f_int_s, nan(nNaNs, 1)'];

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
                f_int(:,iangle) = [f_int_s, nan(nNaNs, 1)'];

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
                f_int(:,iangle) = [f_int_s, nan(nNaNs, 1)'];
                
            end
        end
        
        %-- Compute Intrinsic Frequency Power Spectrum S(f_in, theta) --%

        % Compute f_int differentials via finite differencing
        df_int = diff(f_int(:,iangle));

        % Compute Jacobian 
        J = df_obs./df_int;

        % Compute S_fin by multiplying S_fob by Jacobian
        S_fin(1:length(J),iangle) = S_fob(1:length(J),iangle) .* J;

        % Map S_fin back onto the f_in regular grid (linear interpolation is not extrapolating to values where f_obs > f_int(end))
        idx_nan = ~isnan(S_fin(:,iangle)); 
        S_fin(:,iangle) = interp1(f_int(idx_nan,iangle),S_fin(idx_nan,iangle),f_obs)';
        
        % Find non-nan values in S_fin after mapping (mapping between f_int and
        % f_obs introduces new nans because ranges are incompatible)
        idx_nan = ~isnan(S_fin(:,iangle)); 

        % Set the non-nan values for the intrinsic frequency 
        f_truc = f_obs(idx_nan);

        % Set f_b frequency
        f_b(iangle,2) = f_truc(end);
             
    end  

end
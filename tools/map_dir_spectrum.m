function [S_fin, f_int, S_fob, f_obs, f_b, df_int, J, variance] = map_dir_spectrum(S_fob, f_obs, f_cut, df_obs, dtheta, U, theta_r, tail)
    
    %%%%
    %
    % [S_fin, f_int, S_fob, f_obs, f_b, df_int, J, variance] = map_dir_spectrum(S_fob, f_obs, f_cut, df_obs, dtheta, U, theta_r, tail)
    %
    % Function for mapping observed frequency (frequency measured in the 
    % moving reference frame of the platform) to intrinsic frequency 
    % (frequency measured in the absence of platform motion) for deep water
    % surface gravity waves. The function also transforms the directional 
    % wave power spectrum, which is a function of observed frequency
    % and direction S(f_obs, theta), to a power spectrum as a function of 
    % intrinsic frequency and direction S(f_int, theta). 
    %
    %   Parameters
    %   ----------
    %   S_fob : Wave power spectrum as function of observed frequency and 
    %          direction S(f_obs, theta). Defined nomially as the 
    %          directional-frequency power spectrum computed from vertical
    %          displacement or heave of a platform and horizontal cartesian 
    %          components of velocity. Units: m^2 (Hz deg)^(-1).
    %   f_obs : Observed frequency (i.e. frequency measured in the 
    %           moving reference frame of the platform). Units: Hz. 
    %   f_cut : Noise cutoff frequency for the platform. For the SV3 wave
    %           gliders, f_cut is nominally ~1 Hz (corresponds to waves 
    %           with a wavelength half the length of the waveglider).
    %           Units Hz. 
    %   df_obs : Observed frequency resolution of input spectrum. Units: Hz. 
    %   dtheta : Direction resolution of spectrum. Units: Degrees. 
    %   U : Propagation speed of the observer.  Units: ms^(-1). 
    %   theta_r : Angle between the propagation velocity of the platform 
    %             and the direction of propagation of waves. This quantity 
    %             is a evaluated at all directions of wave propagation.
    %             Units: Degrees, CW, going towards, reference north.
    % 
    %   Returns
    %   -------
    %   S_fin : Wave power spectrum as a function of intrinsic frequency and 
    %           theta. This quantity is multiplied by the jacobian. 
    %           Units: m^2 (Hz deg)^(-1).
    %   f_int : Intrinsic frequency (computed for each direction). 
    %           Units: Hz.
    %   S_fob : Wave power spectrum as a function of observed frequency and
    %           theta. This quantity is the same as the input.
    %           Units: m^2 (Hz deg)^(-1).
    %   f_obs : Observed Frequency. Units: Hz.
    %   f_b : Bifurcation frequency and last frequency not an NaN (blocking 
    %         frequency). The bifircation frequency exists for branches 3-5
    %         (platform moving in direction of wave propagation) and is NaN
    %         for branches 1-2. The bocking frequency is set by the 
    %         interpolation and the bifurcation freq. Bifurcation frequency
    %         is saved in first column f_bif = f_b(:,1) and blocking 
    %         frequency is saved in second column f_block = f_b(:,2). 
    %         Units: Hz. 
    %   df_int : Differential for intrinsic frequency. Units: Hz. 
    %   J : Jacobian. Units: unitless. 
    %   variance : Matrix containing four values for the variance computed
    %              over a given frequency band using rectangular (Riemann) 
    %             integration method: 
    %                   variance = [f_ob_0.01, f_in_0.01, f_ob_0.1, f_in_0.1]
    %             where: 
    %
    %             (1) f_ob_0.01 : Variance computed from integrating 
    %                             S(f_obs)*df_obs from 0.01 to 1 Hz
    %             (2) f_in_0.01 : Total variance computed from integrating
    %                             S(f_int)*df_obs/df_int*df_int from 0.01 to 1 Hz
    %             (3) f_ob_0.1 : Variance computed from integrating 
    %                            S(f_obs)*df_obs from 0.1 to 1 Hz
    %             (4) f_in_0.1 : Variance computed from integrating
    %                            S(f_int)*df_obs/df_int*df_int from 0.1 to 1 Hz
    % 
    %   Notes
    %   -----
    %   (1) Add a spectral tail to the directional spectrum requires
    %       information about the deirectional dependence of the spectral
    %       slope in the equilibirum and saturation ranges. The
    %       conventional spectral forms of the spectra (i.e., f^-4 and f^-5)
    %       do not include this information. Therefore, the functional form
    %       must be either updated to include directional dependence or the
    %       spectral tail attachment must be preformed on the
    %       omni-directional spectrum. I choose the later. 
    %
    %%%%
    
    % Set constants
    g = 9.81;  % Gravitational Acceleration (units: ms^(-2))
    
    % Initialize variables 
    S_fin = NaN(size(S_fob));
    f_int = zeros(size(S_fob)); 
    J = NaN(length(f_obs)-1,length(theta_r));
    f_b = NaN(length(theta_r),2);
    int_S_fob_l = NaN(length(theta_r),1); int_S_fin_l = NaN(length(theta_r),1);
    int_S_fob_h = NaN(length(theta_r),1); int_S_fin_h = NaN(length(theta_r),1);

    %------- Map from observed to intrinsic frequency -------%
    
    % Loop through relative angle 
    for iangle = 1:length(theta_r)
        
        % Call ith relative angle
        itheta_r = theta_r(iangle);
    
        %------- Branch Conditional statements -------%

        %---------- Branch 1 : Moving against waves ----------%
        if cosd(itheta_r) < 0 && U > 0
            %disp('Branch 1 : Moving against waves');

            % Compute intrinsic frequency 
            f_int(:,iangle) = f_obs + (g - 4*pi*f_obs*U*cosd(itheta_r) - sqrt(g^2 - 8*pi*f_obs*g*U*cosd(itheta_r)))/(4*pi* U * cosd(itheta_r));
            
        %---------- Branch 2 : Moving noraml to waves ----------%
        elseif cosd(itheta_r) == 0 || U == 0
            %disp('Branch 2 : Moving normal to waves');

            % Compute intrinsic frequency for branch 2 
            f_int(:,iangle) = f_obs;

        %---------- Branch 3-5 : Moving with wave ----------%
        elseif cosd(itheta_r) > 0 && U > 0
            %disp('Branch 3-5 : Moving with wave');

            % Compute f_obs where df_in(f_obs)/df_ob tends towards infinity  
            f_ob_max = (g)/(8*pi*U*cosd(itheta_r));

            % Compute f_int at the f_obs value where df_in(f_obs)/df_ob
            % tends towards infinity and save in f_b array
            f_in_max = (g)/(4*pi*U*cosd(itheta_r));
            f_b(iangle,1) = f_in_max; 

            % Compute f_int when f_obs = 0
            f_in_zero = (g)/(2*pi*U*cosd(itheta_r));

            %---------- Platform moving slower than energy and crests ----------%
            if f_in_max > f_cut
                % disp('moving slower than energy and crests');

                % Compute the f_ob value where f_int = f_cut for branch 3
                f_ob_cut = -((2*pi*U*cosd(itheta_r)*f_cut^2)/(g)) + f_cut;

                % Obtain f_obs values that map to f_int less than or equal to f_ob_cut
                f_obs_map = f_obs(f_obs <= f_ob_cut);

                % Compute intrinsic frequency 
                f_int_s = f_obs_map + (g - 4*pi*f_obs_map*U*cosd(itheta_r) - sqrt(g^2 - 8*pi*f_obs_map*g*U*cosd(itheta_r)))/(4*pi* U * cosd(itheta_r));
                
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
                f_int_s = f_obs_map + (g - 4*pi*f_obs_map*U*cosd(itheta_r) - sqrt(g^2 - 8*pi*f_obs_map*g*U*cosd(itheta_r)))/(4*pi* U * cosd(itheta_r));

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
                f_int_s = f_obs_map + (g - 4*pi*f_obs_map*U*cosd(itheta_r) - sqrt(g^2 - 8*pi*f_obs_map*g*U*cosd(itheta_r)))/(4*pi* U * cosd(itheta_r));

                % Insert NaNs to extend the length of f_tot to the length of
                % f_in
                nNaNs = length(f_obs) - length(f_int_s); 
                f_int(:,iangle) = [f_int_s, nan(nNaNs, 1)'];
                
            end
        end
        
        %-- Compute Intrinsic Frequency Power Spectrum S(f_int, theta) --%

        % Compute f_int differentials via finite differencing
        df_int = diff(f_int(:,iangle));

        % Compute Jacobian 
        J(:,iangle) = df_obs./df_int;

        % Compute S_fin by multiplying S_fob by Jacobian
        S_fin(1:length(J),iangle) = S_fob(1:length(J),iangle) .* J(:,iangle);

        % Map S_fin back onto the f_in regular grid (linear interpolation is not extrapolating to values where f_obs > f_int(end))
        idx_nan = ~isnan(S_fin(:,iangle)); 
        S_fin(:,iangle) = interp1(f_int(idx_nan,iangle),S_fin(idx_nan,iangle),f_obs)';
        
        % Find non-nan values in S_fin after mapping (mapping between f_int and
        % f_obs introduces new nans because ranges are incompatible)
        idx_nan = ~isnan(S_fin(:,iangle)); 

        % Set the non-nan values for the intrinsic frequency 
        f_truc = f_obs(idx_nan);

        % Set blocking frequency
        f_b(iangle,2) = f_truc(end);
        
        % Obtain indices for integration (non-NaN values at frequencies higher the 0.01 and 0.1 Hz)
        idx_fob_l = ~isnan(S_fob(:,iangle)) & (f_obs >= 0.01)'; idx_fin_l = ~isnan(S_fin(:,iangle)) & (f_int(:,iangle) >= 0.01);
        idx_fob_h = ~isnan(S_fob(:,iangle)) & (f_obs >= 0.1)'; idx_fin_h = ~isnan(S_fin(:,iangle)) & (f_int(:,iangle) >= 0.1);

        % Compute frequency integral for the ith angle excluding NaNs
        % and specified frequency bands from integral

        %----- Frequency band: 0.01 to 1 Hz -----% 
        int_S_fob_l(iangle,1) = sum(S_fob(idx_fob_l,iangle)*df_obs); 
        int_S_fin_l(iangle,1) = sum(S_fin(idx_fin_l,iangle)*df_obs);

        %----- Frequency band: 0.1 to 1 Hz -----% 
        int_S_fob_h(iangle,1) = sum(S_fob(idx_fob_h,iangle)*df_obs); 
        int_S_fin_h(iangle,1) = sum(S_fin(idx_fin_h,iangle)*df_obs);
             
    end  

    % Compute azimuthal integral for variance before and after mapping
    variance = [sum(int_S_fob_l*dtheta), sum(int_S_fin_l*dtheta),sum(int_S_fob_h*dtheta), sum(int_S_fin_h*dtheta)];  
    
end
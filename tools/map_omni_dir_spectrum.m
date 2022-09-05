function [S_fin, f_int, S_fob, f_obs, f_b, df_int, J, variance] = map_omni_dir_spectrum(S_fob, f_obs, f_cut, df_obs, U, theta_r, tail)
    
    %%%%
    % [S_fin, f_int, S_fob, f_obs, f_b, df_int, J, variance] = map_omni_dir_spectrum(S_fob, f_obs, f_cut, df_obs, U, theta_r, tail)
    %
    % Function for mapping observed frequency (frequency measured in the 
    % moving reference frame of the platform) to intrinsic frequency 
    % (frequency measured in the absence of platform motion) for deep water
    % surface gravity waves. Function also transforms the omni-directional 
    % wave power spectrum, which is a function of observed frequency
    % S(f_in) to a power spectrum as a function of intrinsic frequency.  
    %
    %   Parameters
    %   ----------
    %   S_fob : Wave power spectrum as function of observed frequency 
    %           S(f_obs). Defined nomially as the power spectrum of vertical
    %           displacement or heave of a platform. Units: m^2 Hz^(-1).
    %   f_obs : Observed Frequency (i.e. frequency measured in the 
    %           moving reference frame of the platform). Units: Hz. 
    %   f_cut : Noise cutoff frequency for the platform. For the SV3 wave
    %           gliders, f_cut is nominally ~1 Hz (corresponds to waves 
    %           with a wavelength half the length of the waveglider).
    %           Units Hz. 
    %   df_obs : Observed frequency resolution of input spectrum. Units: Hz. 
    %   U : Propagation speed of the observer.  Units: ms^(-1). 
    %   theta_r : Angle between the propagation velocity of the platform
    %             and the direction of propagation of waves. This is
    %             prescribed as a scalar value such that the mean wave
    %             direction is assumed to be aligned with the mean wind 
    %             direction for locally forced waves. Units: Degrees.
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
    %   S_fin : Wave power spectrum as a function of intrinsic frequency. 
    %           This quantity is multiplied by the jacobian. 
    %           Units: m^2 Hz^(-1).
    %   f_int : Intrinsic Frequency. Units: Hz.
    %   S_fob : Wave power spectrum as a function of observed frequency.
    %          This quantity is the same as the input. Units: m^2 Hz^(-1).
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
    %         (1) f_ob_0.01 : Variance computed from integrating 
    %                         S(f_obs)*df_obs from 0.01 to 1 Hz
    %         (2) f_in_0.01 : Total variance computed from integrating
    %                         S(f_int)*df_obs/df_int*df_int from 0.01 to 1 Hz
    %         (3) f_ob_0.1 : Variance computed from integrating 
    %                        S(f_obs)*df_obs from 0.1 to 1 Hz
    %         (4) f_in_0.1 : Variance computed from integrating
    %                        S(f_int)*df_obs/df_int*df_int from 0.1 to 1 Hz
    %
    %%%%
    
    % Set constants
    g = 9.81;  % Gravitational Acceleration (units: ms^(-2))  
    
    % Initialize variables
    S_fin = NaN(size(S_fob));
    f_b = NaN(1,2);
    
    %------- Branch Conditional statements -------%

    %---------- Branch 1 : Moving against waves ----------%
    if cosd(theta_r) < 0 && U > 0
        %disp('Branch 1 : Moving against waves');

        % Compute intrinsic frequency for branch 1
        f_int = f_obs + (g - 4*pi*f_obs*U*cosd(theta_r) - sqrt(g^2 - 8*pi*f_obs*g*U*cosd(theta_r)))/(4*pi* U * cosd(theta_r));

    %---------- Branch 2 : Moving noraml to waves ----------%
    elseif cosd(theta_r) == 0 || U == 0
        %disp('Branch 2 : Moving noraml to waves ');
        
        % Compute intrinsic frequency for branch 2 
        f_int = f_obs; 
        
    %---------- Branch 3-5 : Moving with wave ----------%
    elseif cosd(theta_r) > 0 && U > 0
        %disp('Branch 3-5 : Moving with wave');

        % Compute f_obs where df_in(f_obs)/df_ob tends towards infinity  
        f_ob_max = (g)/(8*pi*U*cosd(theta_r));

        % Compute f_int at the f_obs value where df_in(f_obs)/df_ob
        % tends towards infinity and save in f_b array
        f_in_max = (g)/(4*pi*U*cosd(theta_r));
        f_b(1,1) = f_in_max; 
        
        % Compute f_int when f_obs = 0
        f_in_zero = (g)/(2*pi*U*cosd(theta_r));

        %---------- Platform moving slower than energy and crests ----------%
        if f_in_max > f_cut
            % disp('moving slower than energy and crests');
            
            % Compute the f_in value where f_int = f_cut for branch 3
            f_ob_cut = -((2*pi*U*cosd(theta_r)*f_cut^2)/(g)) + f_cut;

            % Obtain f_obs values that map to f_int less than or equal to f_ob_cut
            f_obs_map = f_obs(f_obs <= f_ob_cut);

            % Compute total frequency 
            f_int_s = f_obs_map + (g - 4*pi*f_obs_map*U*cosd(theta_r) - sqrt(g^2 - 8*pi*f_obs_map*g*U*cosd(theta_r)))/(4*pi* U * cosd(theta_r));

            % Insert NaNs to extend the length of f_int_s to the length of
            % f_obs
            nNaNs = length(f_obs) - length(f_int_s); 
            f_int = [f_int_s, nan(nNaNs, 1)'];

        %---------- Platform moving faster than energy but slower than crests ----------%
        elseif f_tot_max  > 0 && f_tot_max < f_cut && f_in_zero > f_cut
            % disp('moving faster than energy but slower than crests');
            
            % Obtain f_obs values for mapping to branch 3 (f_obs less than or equal to f_ob_max)
            f_obs_map = f_obs(f_obs <= f_ob_max);

            % Compute intrinsic frequency for branch 3
            f_int_s = f_obs_map + (g - 4*pi*f_obs_map*U*cosd(theta_r) - sqrt(g^2 - 8*pi*f_obs_map*g*U*cosd(theta_r)))/(4*pi* U * cosd(theta_r));
            
            % Insert NaNs to extend the length of f_int to the length of
            % f_in
            nNaNs = length(f_obs) - length(f_int_s); 
            f_int = [f_int_s, nan(nNaNs, 1)'];

        %---------- Platform moving faster than energy and crests ----------%
        elseif f_in_zero < f_cut
            % disp('moving faster than energy and crests');
            
            % Obtain f_obs values for mapping to branch 3 (f_obs less than or equal to f_ob_max)
            f_obs_map = f_obs(f_obs <= f_ob_max);

            % Compute intrinsic frequency for branch 3
            f_int_s = f_obs_map + (g - 4*pi*f_obs_map*U*cosd(theta_r) - sqrt(g^2 - 8*pi*f_obs_map*g*U*cosd(theta_r)))/(4*pi*U * cosd(theta_r));
            
            % Insert NaNs to extend the length of f_tot to the length of
            % f_in
            nNaNs = length(f_obs) - length(f_int_s); 
            f_int = [f_int_s, nan(nNaNs, 1)'];
            
        end
    end 
    
    %-- Compute Intrinsic Frequency Power Spectrum S(f_int) --%

    % Compute f_int differentials via finite differencing
    df_int = diff(f_int);

    % Compute Jacobian 
    J = df_obs./df_int;

    % Compute S_fin by multiplying S_fob by Jacobian 
    S_fin(1:numel(J)) = S_fob(1:length(J)) .* J';

    % Map S_fin back onto the f_in regular grid (linear interpolation is not extrapolating to values where f_obs > f_int(end))
    idx_nan = ~isnan(S_fin); 
    S_fin(idx_nan) = interp1(f_int(idx_nan),S_fin(idx_nan),f_obs(idx_nan))';

    % Find non-nan values in S_fin after mapping (mapping between f_int and
    % f_obs introduces new nans because ranges are incompatible)
    idx_nan = ~isnan(S_fin); 

    % Set the non-nan values for the intrinsic frequency and power
    % spectrum to truc variables
    f_truc = f_obs(idx_nan);
    S_truc = S_fin(idx_nan);

    % Set blocking frequency
    f_b(1,2) = f_truc(end);

    %------ Add a spectral tail ------%
    if tail(1) == true 
        
        % Set frequency at which the tail will be attached at
        f_a = f_truc(end);
        
        % Obtain frequencies of tail
        f_tail = [f_a, f_obs(~idx_nan)];

        % disp(['Length of tail: ' num2str(length(f_tail))])
        
        %%%%%%%%%% Frequency Transition Conditional statements %%%%%%%%%%
        % Set transition frequencies
        f_eq = tail(2); % Spectral peak to equilibrium range 
        f_sat = tail(3); % Equilibrium to saturation range 
        
        %-------- f_sat less than noise cutoff frequency  --------%
        if f_cut > f_sat 
            disp('f_sat less than noise cutoff frequency')
            
            %-------- Blocking frequency in Equilibrium range --------%
            if f_a <= f_sat 
                disp('Blocking frequency in Equilibrium range')

                % Preform least squares fit within frequency range 
                % [f_eq f_b] and set the power spectral density at the
                % attachment frequency
                [~, ~, fit, ~] = spectral_slope(f_truc, S_truc, f_eq, f_a);
                S_c = fit(end);

                % Obtain frequency indices for transition frequency and 
                % spectral tails in equilibrium and saturation ranges
                idx_f_eq = f_tail >= f_a & f_tail <= f_sat; 
                idx_f_sat = f_tail >= f_sat; 

                % Obtain frequencies of tail from each frequency range
                f_tail_eq = f_tail(idx_f_eq);
                f_tail_sat = f_tail(idx_f_sat);

                % Compute Spectal tail in linear space
                S_tail_eq = S_c * (f_tail_eq/f_a).^(-4); 
                S_tail_sat = S_tail_eq(end) * (f_tail_sat/f_tail_sat(1)).^(-5); % Note: f_tail_sat(1) instead of f_t is used in the denominator because we want the first PSD level to be equal to PSD level at the last frequency in the equilibrium range.

                % Attach tail onto end of spectrum
                S_fin(:,iangle) = [S_truc; S_tail_eq(2:end)'; S_tail_sat(1:end)'];

            %-------- Blocking frequency in Saturation range --------%
            elseif f_a > f_sat
                disp('Blocking frequency in Saturation range')
                disp([num2str(length(find(f_truc >= f_sat)))])

                %-------- Enough frequencies in saturation range to preform least squares fit  --------%
                if length(find(f_truc >= f_sat)) >= 2

                    % Preform least squares fit within frequency range 
                    % [f_sat f_b] and set the power spectral density at the
                    % attachment frequency
                    [~, ~, fit, ~] = spectral_slope(f_truc, S_truc, f_sat, f_a);
                    S_c = fit(end);

                %-------- NOT Enough frequencies in saturation range to preform least squares fit  --------%
                else 
                
                    % Preform least squares fit within frequency range 
                    % [f_eq f_b] and set the power spectral density at the
                    % attachment frequency
                    [~, ~, fit, ~] = spectral_slope(f_truc, S_truc, f_eq, f_a);
                    S_c = fit(end);
                end 

                % Compute spectal tail 
                S_tail = S_c * (f_tail/f_a).^(-5);

                % Attach tail onto end of spectrum
                S_fin(:,iangle) = [S_truc; S_tail(2:end)'];

            end
            
        %-------- f_sat frequency more than cutoff frequency  --------%
        elseif f_cut < f_sat
            disp('f_sat frequency more than cutoff frequency')

            % Preform least squares fit within frequency range 
            % [f_eq f_b] and set the power spectral density at the
            % attachment frequency
            [~, ~, fit, ~] = spectral_slope(f_truc, S_truc, f_eq, f_a);
            S_c = fit(end);
            
            % Compute spectal tail in equilibrium range
            S_tail = S_c * (f_tail/f_a).^(-4);

            % Attach tail onto end of spectrum
            S_fin(:,iangle) = [S_truc; S_tail(2:end)'];
            
        end    
    end   

    % Obtain indices for integration (non-NaN values at frequencies higher the 0.01 and 0.1 Hz)
    idx_fob_l = ~isnan(S_fob) & (f_obs >= 0.01)'; idx_fin_l = ~isnan(S_fin) & (f_int >= 0.01)';
    idx_fob_h = ~isnan(S_fob) & (f_obs >= 0.1)'; idx_fin_h = ~isnan(S_fin) & (f_int >= 0.1)';

    % Compute variance by interagting over two bands of the spectrum
    variance = [sum(S_fob(idx_fob_l)*df_obs), sum(S_fin(idx_fin_l)*df_obs), sum(S_fob(idx_fob_h)*df_obs), sum(S_fin(idx_fin_h)*df_obs)];
        
end
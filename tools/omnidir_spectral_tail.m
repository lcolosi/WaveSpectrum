function [S_n, f_n, fit, f_fit] = omnidir_spectral_tail(S, f, f_tail, f_eq, f_sat, f_cut)

    %%%%
    % [S_n, f_n, fit, f_fit] = omnidir_spectral_tail(S, f, f_tail, f_eq, f_sat, f_cut)
    %
    % Function for computing the spectral tail for an omni-directional wave
    % spectrum. The functional form of the spectral tail depends on the
    % frequency range. Recall that the equilibirum frequency range of the
    % wave spectrum corresponds to waves where the terms
    % in the statistical equilibirum radiative transfer equation (i.e., 
    % wind-forcing, wave-breaking dissipation, and wave-wave interactions)
    % are balanced. In the saturation frequency range, the primary balance
    % is between the wind input and dissipation from wave breaking. This 
    % theory predicts an f^-4 slope for the equilibrium range and an f^-5 
    % slope for the saturation range. The limits of these ranges are
    % defined by the transition frequencies: 
    %
    %       f_eq = sqrt(2.25)*f_p
    %       f_sat = ((g*sqrt(r))/(2*pi*u))
    %
    % f_eq is the transition from the spectral peak to the equilibrium 
    % range, f_sat is the transition between the  equilibrium and
    % saturation ranges, and f_p is the peak frequency. Additionally, u is
    % the friction velocity (the mean over the period of time the spectrum 
    % is computed), g is the acceleration due to gravity, and r is
    % phillip's constant. 
    %
    %   Parameters
    %   ----------
    %   S : Omni-directional wave spectrum. Units: m^2 Hz^(-1).
    %   f : Frequencies below the spectral tail. Units: Hz. 
    %   f_tail : Frequencies of the spectral tail including the frequency 
    %            the tail will be attched at. Units: Hz. 
    %   f_eq : Spectral peak to equilibirum range transition frequency. 
    %           Units: Hz.
    %   f_sat : Equilibirum to saturation range transition frequency. 
    %           Units: Hz.
    %   f_cut : Noise cutoff frequency. Frequencies below this are not
    %           properly resolved and should be discarded. 
    % 
    %   Returns
    %   -------
    %   S_n : Omni-directional wave spectrum with the spectral tail
    %         attached. Units: m^2 Hz^(-1)
    %   f_n : Frequencies with spectral tail frequencies attached. 
    %         Units: Hz. 
    %   fit : Least squares fit for spectral tail. Units: m^2 Hz^(-1). 
    %   f_fit : Frequency subrange for least square fit. Units: Hz. 
    %   
    %   Notes
    %   -----
    %   (1) This function only fits a high frequency tail. That is, a
    %       spectral that extends the spectrum to higher frequencies. 
    %
    %   (2) This function accounts for the different cases of the
    %       attachment frequency's position in frequency space relative to
    %       the transition frequencies. This is accomplished in the 
    %       frequency transition conditional statements. 
    %
    %   (3) The spectral slope is attached to omni-directial wave spectrum
    %       at a power spectral density level set by least-squares fitting
    %       a line through a range of frequencies and grabbing value of the
    %       fit corresponding to the attachment frequency. This helps to
    %       lessen the impact the noise in the power spectrum at high 
    %       frequency has on the spectral tail.     
    %   
    %%%%
            
    % Make both frequency and omni-directional spectrum column vectors
    f = f(:); 
    f_tail = f_tail(:);
    S = S(:);

    % Set frequency at which the tail will be attached at
    f_a = f(end);

    % Attach frequencies on spectral tail to frequency vector
    f_n = [f; f_tail(2:end)];
    
    %-------- f_sat less than noise cutoff frequency  --------%
    if f_cut > f_sat 
        
        %-------- Attachment frequency in Equilibrium range --------%
        if f_a <= f_sat 
            % disp('Fit tail in both equilibrium and saturation ranges')

            % Preform least squares fit within frequency range 
            % [f_eq f_a] and set the power spectral density at the
            % attachment frequency
            [~, ~, fit, f_fit] = spectral_slope(f, S, f(end-10), f_a);
            S_c = fit(end);

            % Obtain frequency indices for transition frequency and 
            % spectral tails in equilibrium and saturation ranges
            idx_f_eq = f_tail >= f_a & f_tail <= f_sat; 
            idx_f_sat = f_tail >= f_sat; 

            % Obtain frequencies of tail from each frequency range
            f_tail_eq = f_tail(idx_f_eq);
            f_tail_sat = f_tail(idx_f_sat);

            % Compute spectal tail in linear space
            S_tail_eq = S_c * (f_tail_eq/f_a).^(-4); 
            S_tail_sat = S_tail_eq(end) * (f_tail_sat/f_tail_sat(1)).^(-5); % Note: f_tail_sat(1) instead of f_t is used in the denominator because we want the first PSD level to be equal to PSD level at the last frequency in the equilibrium range.

            % Attach tail onto end of spectrum
            S_n = [S; S_tail_eq(2:end); S_tail_sat(1:end)];

        %-------- Attachment frequency in Saturation range --------%
        elseif f_a > f_sat
            % disp('Fit tail in saturation ranges')

            %-------- Enough frequencies in saturation range to preform least squares fit  --------%
            if length(find(f >= f_sat)) >= 2

                % Preform least squares fit within frequency range 
                % [f_sat f_a] and set the power spectral density at the
                % attachment frequency
                [~, ~, fit, f_fit] = spectral_slope(f, S, f_sat, f_a);
                S_c = fit(end);

            %-------- NOT Enough frequencies in saturation range to preform least squares fit  --------%
            else 
            
                % Preform least squares fit within frequency range 
                % [f_eq f_a] and set the power spectral density at the
                % attachment frequency
                [~, ~, fit, f_fit] = spectral_slope(f, S, f(end-10), f_a);
                S_c = fit(end);

            end 

            % Compute spectal tail 
            S_tail = S_c * (f_tail/f_a).^(-5);

            % Attach tail onto end of spectrum
            S_n = [S; S_tail(2:end)];

        end
        
    %-------- f_sat frequency more than cutoff frequency  --------%
    elseif f_cut < f_sat
        % disp('Fit spectral tail in equilibrium range.')

        % Preform least squares fit within frequency range 
        % [f_eq f_b] and set the power spectral density at the
        % attachment frequency
        [~, ~, fit, f_fit] = spectral_slope(f, S, f(end-10), f_a);
        S_c = fit(end);
        
        % Compute spectal tail in equilibrium range
        S_tail = S_c * (f_tail/f_a).^(-4);

        % Attach tail onto end of spectrum
        S_n = [S; S_tail(2:end)];
        
    end
end
function [m, mm, yfit, f_range] = spectral_slope(f, E, fmin, fmax)
    %%%%
    % spectral_slope(f, E, fmin, fmax)
    %
    % Function for computing the spectral slope in log space of power spectral
    % density function E with uncertainty stdE in a frequency subrange
    % [fmin, fmax].
    %
    %   Parameters
    %   ----------
    %   f : Cyclical frequency for the power spectral density function E.
    %   E : Normalized power spectral density function
    %   fmin : Lower frequency limit for frequency subrange.
    %   fmax : Upper frequency limit for frequency subrange.
    % 
    %   Returns
    %   -------
    %   m : Spectral Slope.
    %   mm : Uncertainty of Spectral Slope.
    %   yfit : Unweighted least square fit model
    %   f_range : Frequency subrange out of log10 space
    %   
    %%%%
    
    % Find fequency indices in subrange
    idx_freq = f>=fmin & f<=fmax;
    idx_nan = ~isnan(E); 
    idx = logical(idx_freq.*idx_nan);
    
    % Extract f and E values within subrange in log10 space
    fi = log10(f(idx))';
    Ei = log10(E(idx))';
    
    % Initialize model matrix for computing least squares fit
    A = cat(2, ones(numel(fi),1) , fi);
    
    % Compute coefficients of least squares fit 
    coef = (inv(A'*A))*A'*Ei';
    
    % Grab coefficients for the spectral slope
    m = coef(2);
    
    %  Compute unweighted least square fit model
    yfit = A*coef;
    
    % Map frequency subrange and yfit back to linear space:
    f_range = 10.^(fi);
    yfit = 10.^yfit;
    
    % Initialize variables
    M = length(Ei);
    parameters = 2;
    
    % Compute uncertainty of data
    data_sigma = [];
    data_sigma(:,1) = sqrt((1/(M - parameters)) * sum((Ei - yfit).^2));  % Uncertainty in Log space
    data_sigma(:,2) = sqrt((1/(M - parameters)) * sum((10.^(Ei) - 10.^(yfit)).^2));  % Uncertainty in linear space

    % Compute uncertainty in spectral slope
    model_sig = [];
    delta = M * sum(f_range.^2) - (sum(f_range))^2;
    model_sig(:,1) = (data_sigma(1) * sqrt( sum(f_range.^2)/delta));  % Uncertainty in y-intercept
    model_sig(:,2) = (data_sigma(1) * sqrt(M/delta));  % Uncertainty in slope
    
    % Set slope uncertainty: 
    mm = model_sig(:,2);
    
end
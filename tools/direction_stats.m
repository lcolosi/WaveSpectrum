function [m_theta,std_theta, stdm_theta] = direction_stats(theta, dt, task)

    %%%%
    % [m_theta,std_theta, stdm_theta] = direction_stats(theta, dt, task)
    %
    % Function for computing the mean and standard deviation of an array of
    % direction measurements. 
    %
    %   Parameters
    %   ----------
    %   theta : A vector of angles in units of degrees. NaN may be present 
    %           in vector. Angles should be in the desired directional 
    %           convention.
    %   dt : Temporal or spatial separation between theta observations.
    %        Scalar quantity.
    %   task : Specifies whether the uncertainty of the mean should be
    %          computed. If true, uncertainty is computed. If false, then
    %          uncertainty is not computed. 
    % 
    %   Returns
    %   -------
    %   m_theta : Mean direction.
    %   std_theta : Standard deviation of direction. 
    %   stdm_theta : Standard error of the mean direction. 
    %   
    %%%%
    
    % Compute the cosine and sine components of the angle: 
    theta_cos = cosd(theta);
    theta_sin = sind(theta);

    % Compute the mean of theta_cos and theta_sin
    mtheta_cos = mean(theta_cos, 'omitnan');
    mtheta_sin = mean(theta_sin, 'omitnan');
    
    % Compute mean direction 
    m_theta = mod(atan2d(mtheta_sin,mtheta_cos),360); 
    
    % Compute eps
    mtheta_eps = sqrt(1-(mtheta_cos.^2 + mtheta_sin.^2));

    % Compute standard deviation of direction 
    std_theta = asind(mtheta_eps) * (1+(2/sqrt(3)-1)*mtheta_eps.^3);
    
    % Case 1: Compute uncertainty of the mean estimate 
    if task
        
        % Compute N effective (number of independent observations)
        [N_eff,~] = decor_scale(theta, dt);

        % Compute standard error of the mean (assume no serial correlation)
        stdm_theta = std_theta/sqrt(N_eff);
        
    % Case 2: Do not compute uncertainty of the mean. 
    else 
        
        % Set uncertainty to a flag
        stdm_theta = 'not computed';
        
    end

end


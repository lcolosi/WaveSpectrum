function [m_theta,std_theta, stdm_theta] = direction_stats(theta, dt, task)

    %%%%
    % [m_theta,std_theta, stdm_theta] = direction_stats(theta, dt, task)
    %
    % Function for computing the mean, standard deviation, and standard 
    % error of the mean for a time series of direction measurements. 
    %
    %   Parameters
    %   ----------
    %   theta : A vector of angles in units of degrees. NaN may be present 
    %           in vector. Angles should be in the desired directional 
    %           convention.
    %   dt : Temporal or spatial separation between theta observations.
    %        Scalar quantity.
    %   task : Specifies How the standard error of the mean is computed. If
    %          true, standard error is computed assuming data has 
    %          serial correlation. If false, then standard error is
    %          computed assuming all data points are independent. 
    % 
    %   Returns
    %   -------
    %   m_theta : Mean direction (Yamartino Method).
    %   std_theta : Standard deviation of direction (Yamartino Method).
    %   stdm_theta : Standard error of the mean direction. 
    % 
    %   Notes
    %   -----
    %   (1) The mean and standard deviation are computed using the 
    %       Yamartino Method. For more information, see:
    %
    %           (a) Wiki page - https://en.wikipedia.org/wiki/Yamartino_method
    %           (b) Publication - https://urldefense.com/v3/__https://journals.ametsoc.org/view/journals/apme/23/9/1520-0450_1984_023_1362_acospe_2_0_co_2.xml?tab_body=pdf__;!!Mih3wA!AMtkrURmRZGuK8lwOc1Zx2ixE7Zqo3CSEC3tVgpEC2X6gC9HD4UJu49dcrL_QgeyJutdnXi-DAqbcNE$
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
    
    % Case 1: Serial correlation
    if task
        
        % Compute N effective (number of independent observations)
        [N_eff,~] = decor_scale(theta, dt);

        % Compute standard error of the mean 
        stdm_theta = std_theta/sqrt(N_eff);
        
    % Case 2: No serial correlation
    else 
        
        % Compute N (number of observations assuming each are independent)
        N = length(theta);

        % Compute standard error of the mean 
        stdm_theta = std_theta/sqrt(N);
        
    end

end


function [ground_speed,true_course] = platform_speed_dir(ve, vn, dir_con)

    %%%%
    % [ground_speed,true_course] = platform_speed_dir(ve, vn, dir_con)
    %
    % Function for computing a platform's propagation speed and direction
    % given the velocity components of platform in the east-west and
    % north-south directions. 
    %
    %   Parameters
    %   ----------
    %   ve : Eastward component of velocity. May contain NaN values. 
    %        Units: ms^-1. 
    %   vn : Northward component of velocity. May contain NaN values.
    %        Units: ms^-1.
    %       %   dir_con : Specifies the directional convention for the 
    %             platform's progation direction (nov.true_course). Format 
    %             of the parameter is: dir_con = [rotation direction;
    %             coming from/going towards; zero angle reference direction]
    %             Options include: 
    %                 (1) rotation direction = 'CW' (Clockwise), 'CCW' (counter clockwise) 
    %                 (2) coming from/going towards = 'cf', 'gt'
    %                 (3) zero angle reference direction = 'rn' (refernece
    %                     north), 're' (reference east) 
    % 
    %   Returns
    %   -------
    %   ground_speed : Ground speed of the platform. Units: ms^-1. 
    %   true_course : Direction of propagation. Units: degrees.   
    %   
    %%%%
    
    % Compute ground speed
    ground_speed = sqrt(ve.^2 + vn.^2); 
    
    % Compute true course 
    %---------- CW, going towards, ref north ----------%
    if strcmp(dir_con(1), 'CW') && strcmp(dir_con(2), 'gt') && strcmp(dir_con(3), 'rn')
        true_course = mod(90 - atan2d(vn,ve), 360);

    %---------- CW, coming from, ref north ----------%
    elseif strcmp(dir_con(1), 'CW') && strcmp(dir_con(2), 'cf') && strcmp(dir_con(3), 'rn')
        true_course = mod(270 - atan2d(vn,ve), 360);
    end
    
end


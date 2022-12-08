function [ustar, z0, U10] = FUNC_ustar_Z(z, Uz, varargin)

    %%%%
    % [ustar, z0, U10] = FUNC_ustar_Z(z, Uz, varargin)
    %
    % Function for computing the friction velocity, surface roughness
    % length, and the wind speed 10 meters above the ocean surface by
    % assuming the wind speed profile, U(z), near the ocean surface is 
    % logarithmic. This logarithmic wind profile is based off the law of
    % the wall which states: 
    % 
    %   " The average velocity of a turbulent flow at a certain point is
    %     proportional to the lagarithm of the distance from that point to
    %     the "wall" or boundary of the fluid region." 
    % 
    % The Law of the Wall is applicable when "close" to the wall or < 20%
    % of the height of the flow. The mean wind speed (horizontally
    % averaged) at the height z above the fluid boundary is: 
    % 
    %   Uz = ustar/kappa * [ln(z-d/z0) + psi(z,z0,L)]
    % 
    % See Charnock 1955 "Wind Stress on a Water Surface" for more
    % information. 
    % 
    % !!! Important !!! : The solution of u_star, z0 and U10 assumes (1) a
    % logarithmic wind speed profile and (2) a neutrally stable atmosphere.
    % Also, the law of the wall only applies to mean wind speed profiles. 
    %
    %   Parameters
    %   ----------
    %   z : Altitude of wind measurement from the ocean surface. Units: m.
    %
    %   Uz : Mean wind speed measured at z. Uz can be a vector with time or
    %        horizontal distance as the non-unity dimension. Units: ms^{-1}. 
    % 
    %   Varargins (Optional parameters) :
    %       (1) No Charnock-wave relationship
    %               no varargin 
    %
    %       (2) Taylor & Yelland (2001) Model
    %               varargin(1) = 'TY' : Surface roughness model 
    %               varargin(2) = Hs : Significant wave height 
    %               varargin(3) = Fp : Peak frequency 
    %
    %       (3) Oost et al. (2002) Model 
    %               varargin(1) = 'Oo' : Surface roughness model
    %               varargin(2) = Fp : Peak frequency
    % 
    %   Returns
    %   -------
    %   ustar : Friction velocity in the air. Defined as the form by which
    %           shear stress (or surface drag) may be written in units of
    %           velocity. The firction velocity is dependent on the shear
    %           stress at the boundary of the flow and the fluid density. 
    %           Units: ms^(-1). 
    %
    %   z0 : Surface roughness length. Defined as a parameter to account
    %        for the effect of the surface roughness on the wind flow. 
    %        Therefore, the surface elevation of the fluid boundary. This
    %        dependence on the sea surface elevation leads to the
    %        charnock-wave relations specified in the varagin. Units: m. 
    % 
    %   U10 : Wind speed 10 meters above the ocean surface. Units: ms^(-1). 
    %   
    %   Notes
    %   -----
    %   (1) Varargin depends on the model of Charnock coefficient's
    %       dependence on the wave state. Two models are used here and
    %       their surface roughness lengths are defined as: 
    % 
    %           (1) Taylor & Yelland (2001) Model: 
    %                   z0 = 1200*Hs*(Hs/Lp)^4.5 + 0.11*nu/ustar
    %            
    %           (2) Oost et al. (2002) Model: 
    %                   z0 = 50/(2*pi)*Lp*(ustar/Cp)^4.5 + 0.11*nu/ustar
    %  
    %           (3) No charnock-wave Model: 
    %                   z0 = alpha*ustar^2/g + 0.11*nu
    %                   
    %       Here, Hs is significant wave height, Lp is the peak wavelength,
    %       Cp is the phase speed at the peak of the spectrum, g is
    %       gravitational acceleration, alpha is the charnock parameter 
    %       and nu is .
    % 
    %       The Charnock-wave models can produce weird result, so check 
    %       solution with no charnock-wave model output.     
    %   
    %   (2) In matlab, 1d-4 = 10^-4 and 15d-6 = 1.5*10^-5. The notation 
    %       nd-m where n and m are any positive integers provides a more 
    %       concise form compared to exponetial notation. 
    % 
    %   (3) u_star is computed using a difference equation because
    %       substitution leads to a transcendental equation for u_star. 
    %       This equation is NOT easy to solve. The difference equation 
    %       method provides a clever way to sovle this system of nonlinear
    %       equations. The logic behind the validity of the difference 
    %       equation approach is the following. The surface roughness 
    %       length can be isolated using the law of the wall and charnock 
    %       relation equation:        
    % 
    %           z0 = z/exp(kappa*Uz/ustar) --> Law of the Wall 
    %           z0 = alpha*ustar^2/g + 0.11*nu --> Charnock Relation
    %           
    %       Both of these estimates of z0 are correct, so nomially with the
    %       correct constants kappa, alpha, nu, and g and parameters z, Uz,
    %       and ustar, these z0 estimate should agree (i.e. equivalence). 
    %       If there is disagreement, one of the constants or parameters
    %       is incorrect. We know all constants and parameters except
    %       ustar. However, by computing z0 for a range of ustar values and
    %       finding where there is agreement between the estimates, we can
    %       find the proper ustar value such that the physical equations are
    %       statisfied. 
    %%%%

    % Set parameters
    x = 1d-4:1d-4:2;                                                        % Range of ustar values for computing z0 (Units: ms^(-2))
    nu = 15d-6;                                                             % Proportionality constant for the second term in the charnock relation (Units of 0.11*nu = m^2*s^(-1))
    alpha = 0.011;                                                          % Charnock Parameter (Units: dimensionless)
    g = 9.81;                                                               % Gravitational acceleration (Units: ms^(-2))
    kappa = 0.4;                                                            % Von Karman constant (Units: dimensionless)
    
    % Initialize output variables 
    ustar=NaN(length(Uz),1);                                                % Friction velocity (Units: ms^(-1))
    z0=NaN(length(Uz),1);                                                   % Surface roughness length (Units: m)
    U10=NaN(length(Uz),1);                                                  % Wind speed at 10 meters above the sea surface (Units: ms^(-1))

    % Loop through Uz vector
    for i=1:length(Uz)

        % Compute z0 using Law of the Wall for range of u_star values
        fustar1 = z./exp(kappa.*Uz(i)./x);

        %%--------- No charnock-wave Model ---------%%
        if isempty(varargin)

            % Compute z0 using charnock relation for range of u_star values
            fustar2 = alpha*x.^2/g + 0.11*nu./x;

        %%--------- Taylor & Yelland (2001) Model ---------%%
        elseif sum(varargin{1}=='TY')==2

            % Set wave state parameters
            Hs=varargin{2}(i);
            Fp=varargin{3}(i);

            % Compute peak wavelength
            Lp=g./(2*pi*Fp.^2);

            % Compute z0 using charnock-wave relation for range of u_star values
            fustar2 = 1200*Hs*(Hs./Lp)^4.5 + 0.11*nu./x;

        %%--------- Oost et al. (2002) Model ---------%%
        elseif sum(varargin{1}=='Oo')==2

            % Set wave state parameter
            Fp=varargin{2}(i);

            % Compute peak wavelength and phase speed
            Lp=g./(2*pi*Fp.^2);
            Cp=g./(2*pi*Fp);

            % Compute z0 using charnock-wave relation for range of u_star values
            fustar2=50/(2*pi)*Lp*(x./Cp).^4.5 + 0.11*nu./x;
       end
    
        % Compute the difference between z0 estimates
        fustar=fustar1-fustar2;

        % Recompute z0 using law of the wall estimate
        Z0=z./exp(kappa.*Uz(i)./x);
    
        % Find indices where difference between estimates vanish
        I=find(fustar==0);

        % Case 1: Only one u_star value has a difference between estimates that vanishes
        if length(I)==1

            % Set u_star and z0 (law of the wall) estimates at index 
            ustar(i)=x(I);
            z0(i)=Z0(I);

        % Case 1: Multiple u_star value have a difference between estimates that vanishes
        elseif isempty(I)

            % Find indices where the law of the wall estimates is greater
            % than the charnock relation estimate
            J=find(fustar>0);

            % Case 1: If no law of the wall estimates are greater
            if isempty(J)

                % Display loop variable
                disp(i)

            % Case 2: If there are multiple law of the wall estimates are greater
            else

                % Compute the mean of ustar and corresponding z0 (law of
                % the wall) of the first greater than zero fustar index and
                % one index below
                ustar(i)=mean([x(J(1)-1) x(J(1))]);
                z0(i)=mean([Z0(J(1)-1) Z0(J(1))]);
            end
        end

        % Compute the wind speed 10 meters above the ocean surface using
        % the u_star and z0 estimate above. 
        U10(i)=ustar(i)./kappa*log(10/z0(i));

        % Clear variables 
        clear Z0 fustar1 fustar2 fustar
    end
end

%% Developmental Code

% Questions for Luc
% -----------------
%   (1) How did you pick the range of u_star values in the x variable?
%       Because the size of the x array is so large (20,000 elements), this
%       is where most of the time expensive computing occurs
%       especially if you are iterating over multiple Uz values. 
%
%   (2) Where does the second term in the charnock relation for z0 come
%       from? What are parameters nu and 0.11? 
%
%   (3) Why do you find where fustar > 0 for the case where multiple u_star
%       values provide a difference between estimates equal to zero? This
%       is, why should the law of the wall estimate be greater than the
%       charnock relation estimate for z0? 
% 
%   (4) Why do you take the mean of u_star between the first index
%       in the J vector and first index in the J vector minus 1? I would
%       think you want to average all possible values? 


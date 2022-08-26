function [ustar,z0,U10]=FUNC_ustar_Z(z,Uz,varargin)

%  Input:
% z = altitude of wind measurement
% Uz = mean wind speed measured at z.

% Varargin depends on the model to rely Charnock coefficient to wave state
% Charnock-wave model can produce weird result --> check with no wave model

% No Charnock-wave relationship
% z0=alpha*ustar^2/g + 0.11*nu
% no varargin.

% Model Taylor & Yelland (2001): Required Hs and Fp.
% z0=1200*Hs*(Hs/Lp)^4.5 + 0.11*nu/ustar
% varargin(1)='TY'
% varargin(2)=Hs
% varargin(3)=Fp

% Model Oost et al. (2002): Required Fp.
% z0=50/(2*pi)*Lp*(ustar/Cp)^4.5 + 0.11*nu/ustar
% varargin(1)='Oo'
% varargin(2)=Fp

% Output
% Self-explanatory.

x=1d-4:1d-4:2;
nu=15d-6;
alpha=0.011;
g=9.81;
kappa=0.4;

ustar=NaN(length(Uz),1);
z0=NaN(length(Uz),1);
U10=NaN(length(Uz),1);
for i=1:length(Uz)
    fustar1=z./exp(kappa.*Uz(i)./x);
    if length(varargin)==0
        fustar2=alpha*x.^2/g+0.11*nu./x;
    elseif sum(varargin{1}=='TY')==2
        Hs=varargin{2}(i);
        Fp=varargin{3}(i);
        Lp=g./(2*pi*Fp.^2);
        fustar2=1200*Hs*(Hs./Lp)^4.5 + 0.11*nu./x;
    elseif sum(varargin{1}=='Oo')==2
        Fp=varargin{2}(i);
        Lp=g./(2*pi*Fp.^2);
        Cp=g./(2*pi*Fp);
        fustar2=50/(2*pi)*Lp*(x./Cp).^4.5 + 0.11*nu./x;
   end

    fustar=fustar1-fustar2;
    Z0=z./exp(kappa.*Uz(i)./x);

    I=find(fustar==0);
    if length(I)==1
        ustar(i)=x(I);
        z0(i)=Z0(I);
    elseif length(I)==0
        J=find(fustar>0);
        if length(J)==0
            disp(i)
        else
            ustar(i)=mean([x(J(1)-1) x(J(1))]);
            z0(i)=mean([Z0(J(1)-1) Z0(J(1))]);
        end
    end
    U10(i)=ustar(i)./kappa*log(10/z0(i));
    clear Z0 fustar1 fustar2 fustar
end
end

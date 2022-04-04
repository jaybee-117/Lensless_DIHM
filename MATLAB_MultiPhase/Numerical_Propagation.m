%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [Uz,I,Phi]  = Numerical_Propagation(U0,dz,Lambda)
% Inputs :
%	U0: Orignal complex field;
%	dz: Propagation distance(um);
%   Lambda: Wavelength of light source (nm);
% Outputs:
%   Uz: Complex field after propagation;
%   I: Intesity field after propagation;
%   Phi: Phase after propagation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uz,Iz,Phiz]  = Numerical_Propagation(U0,dz,Lambda)

[Ny,Nx]=size(U0);
Lambda = Lambda*10^-9;

k=2*pi/Lambda; % wave number

Fx = (-(Nx/2):1:(Nx/2-1))/Nx; % frequency coordinate X.
Fy = (-(Ny/2):1:(Ny/2-1))/Ny; % frequency coordinate Y.

[U,V] = meshgrid(Fx,Fy);

    Exp_term = sqrt(1-(Lambda*U).^2-(Lambda*V).^2);
    H = exp(1i*k*dz*Exp_term);
    H((1-(Lambda*U).^2-(Lambda*V).^2)<0) = 0; % neglect evanescent wave
    Uz=ifft2(ifftshift(U0.*H));
    Uz=ifftshift(Uz);
    
    Iz = abs(Uz);
    Phiz = angle(Uz);
end

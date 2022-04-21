%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [Uz,I,Phi]  = Numerical_Back_Propagation(U0,dz,Lambda)
% Inputs :
%	U0: Orignal complex field;
%	dz: Propagation distance(um);
%   Lambda: Wavelength of light source (nm);
% Outputs:
%   Uz: Complex field after propagation;
%   I: Intesity field after propagation;
%   Phi: Phase after propagation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uz,Iz,Phiz]  = Numerical_Back_Propagation(U0,dz,Lambda)

[Ny,Nx]=size(U0);
L = 1000; %FOV in microns
dx = L/Nx;
x = -L/2:dx:L/2 - dx;
dy = L/Ny;
y = -L/2:dy:L/2 - dy;

fx = (-(Nx/2):1:(Nx/2-1))/Nx*dx; % frequency coordinate X.
fy = (-(Ny/2):1:(Ny/2-1))/Ny*dy; % frequency coordinate Y.
[Fx,Fy] = meshgrid(fx,fy);

Lambda = Lambda*10^-3; %Wavelength of light in microns
k=2*pi/Lambda; % wave number

U1 = fftshift(fft2(U0));
Exp_term = sqrt(k^2-(2*pi*Fx).^2-(2*pi*Fy).^2);
H = exp(-1i*dz*Exp_term);
H((1-(Lambda*Fx).^2-(Lambda*Fy).^2)<0) = 0; % neglect evanescent wave
Uz = ifft2(ifftshift(U1.*H));
Iz = abs(Uz).^2;
Phiz = angle(Uz);
end

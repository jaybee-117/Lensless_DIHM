%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical propagation the complex field to another plane at a
% given distance using 'Angular Specturm' or 'Fresnel' method
% or 'TIE' method (TIE only predict intensity distribution).
% Version 1.0 - initial version;
% Version 1.5 - consider evanescent wave;
% Version 1.5 - rewrite for better readablility;
% Version 2.0 - add TIE method;
% last modified on 02/27/2020
% by Chao Zuo (zuochao@njust.edu.cn,surpasszuo@163.com)
% _________________________________________________________________________
% [Uz,I,Phi]  = Numerical_Propagation(U0,dz,Pixelsize,k,Method)
% Inputs :
%	U0: Orignal complex field;
%	dz: Propagation distance;
%	Pixelsize: Pixelsize (in the sample space);
%   k: Wave-number;
%   Method: type of method used ('Angular Spectrum' or 'Fresnel').
% Outputs:
%   Uz: Complex field after propagation;
%   I: Intesity field after propagation;
%   Phi: Phase after propagation.
% _________________________________________________________________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uz,Iz,Phiz]  = Numerical_Propagation(U0,dz,Pixelsize,Lambda,Method)

[Ny,Nx]=size(U0);

k=2*pi/Lambda; % wave number

Fx = (-(Nx/2):1:(Nx/2-1))/Nx/Pixelsize; % frequency coordinate X.
Fy = (-(Ny/2):1:(Ny/2-1))/Ny/Pixelsize; % frequency coordinate Y.

[U,V] = meshgrid(Fx,Fy);

FU0=fftshift(fft2(U0));

if(strcmp(Method,'Angular Spectrum'))   % Angular Specturm method
    Exp_term = sqrt(1-(Lambda*U).^2-(Lambda*V).^2);
    H = exp(1i*k*dz*Exp_term);
    H((1-(Lambda*U).^2-(Lambda*V).^2)<0) = 0; % neglect evanescent wave
    Uz=ifft2(ifftshift(FU0.*H));

elseif(strcmp(Method,'Fresnel'))        % Fresnel method
    H=exp(1i*k*dz*(1-((Lambda*U).^2+(Lambda*V).^2)/2));
    Uz=ifft2(ifftshift(FU0.*H));

elseif(strcmp(Method,'TIE'))            % TIE method
    I0 = U0.*(conj(U0));
    fx = fftshift(U);
    fy = fftshift(V);
    
    % Terms for derivative calculation in frequency domain
    Cx = 2*1i*pi*fx;
    Cy = 2*1i*pi*fy;
    
    % Calculate the estimate of dIdz using TIE
    Fphi = fft2(angle(U0));
    Fdphidx = Fphi.*Cx;
    Fdphidy = Fphi.*Cy;
    dphidx = real(ifft2(Fdphidx));
    dphidy = real(ifft2(Fdphidy));
    
    dPsidx = I0.*dphidx;
    dPsidy = I0.*dphidy;
    
    FdPsidx = fft2(dPsidx);
    FdPsidy = fft2(dPsidy);
    
    Fd2Psidx2 = FdPsidx.*Cx;
    Fd2Psidy2 = FdPsidy.*Cy;
    
    d2Psidx2 = real(ifft2(Fd2Psidx2));
    d2Psidy2 = real(ifft2(Fd2Psidy2));
    
    laplacePsi = d2Psidx2+d2Psidy2;
    
    dIdz = laplacePsi./(-k);
else
    errordlg('Type of transfer function must be ... <Angular Spectrum> or <Fresnel> or <TIE>','Error');
end

if(strcmp(Method,'TIE'))            % TIE method
    Iz = dz*dIdz+I0;
    Uz = NaN;
    Phiz = NaN;
else
    Iz = abs(Uz).^2;
    Phiz = angle(Uz);
end

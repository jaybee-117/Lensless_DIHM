%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Phi, dIdz_est] = TIE_MAX_solver(Phi_curr,dIdz_curr,dz,I0,Pixelsize,k,r,Method)
% Purpose: one-step solution to TIE using FFT-based Possion solver assuming
%          uniform maximum intensity
%'Inputs': 'Phi_curr', Current phase estimate;
%          'dIdz_curr', Current intensity derivative estimate;
%'         'dz',Propagation distance;
%          'I0', Infoucs intensity image
%          'pixelsize'��pixelsize('mm');
%          'k', Wavenumber
%          'r', regularzation parameter (to remove low-frequency noise)
%'Outputs': 'phi',object phase map;
%           'dIdz_est', estimated dIdz;
% ！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
% last modified on 02/27/2020
% by Chao Zuo (zuochao@njust.edu.cn,surpasszuo@163.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi, dIdz_est] = TIE_MAX_solver(Phi_curr,dIdz_curr,dz,I0,Pixelsize,k,r,Method)

J = -k*dIdz_curr;
[Ny,Nx] = size(J);
Fx = (-(Nx/2):1:(Nx/2-1))/Nx/Pixelsize; % frequency coordinate X.
Fy = (-(Ny/2):1:(Ny/2-1))/Ny/Pixelsize; % frequency coordinate Y.
[U,V] = meshgrid(Fx,Fy);
Lambda = 2*pi/k;

% Coordinates in frequency domain
fx = fftshift(U);
fy = fftshift(V);
% Terms for derivative calculation in frequency domain
Cx = 2*1i*pi*fx;
Cy = 2*1i*pi*fy;
% Calculate phi assuming uniform intensity.
FJ = fft2(J./max(max(I0)));
Fphi = FJ.*(Cx.*Cx+Cy.*Cy)./(r/Pixelsize.^4+(Cx.*Cx+Cy.*Cy).^2);
Phi = real(ifft2(Fphi));

if(strcmp(Method,'TIE'))
    % Calculate the estimate of dIdz using TIE
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
    % Estimate dIdz
    dIdz_est = laplacePsi./(-k);
    
elseif(strcmp(Method,'Angular Spectrum'))
    % Calculate the estimate of dIdz using Angular Spectrum
    FU=fftshift(fft2(sqrt(I0).*exp(1i*(Phi+Phi_curr))));
    Exp_term = sqrt(1-(Lambda*U).^2-(Lambda*V).^2);
    H = exp(1i*k*dz*Exp_term);
    H((1-(Lambda*U).^2-(Lambda*V).^2)<0) = 0; % neglect evanescent wave
    Uz=ifft2(ifftshift(FU.*H));
    Iz=Uz.*conj(Uz);
    dIdz_est = (Iz-I0)/dz;
    
elseif(strcmp(Method,'Fresnel'))
    % Calculate the estimate of dIdz using Fresnel method
    FU=fftshift(fft2(sqrt(I0).*exp(1i*(Phi+Phi_curr))));
    H=exp(1i*k*dz*(1-((Lambda*U).^2+(Lambda*V).^2)/2));
    Uz=ifft2(ifftshift(FU.*H));
    Iz=Uz.*conj(Uz);
    dIdz_est = (Iz-I0)/dz;
    
else
    errordlg('Type of transfer function must be <Angular Spectrum> or <Fresnel>','Error');
end

end

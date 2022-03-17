%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Phi, dIdz_est] = TIE_solver_one_step(Phi_curr,dIdz,dz,I0,Pixelsize,k,r,Method)
% Purpose: one-step solution to TIE using FFT-based Possion solver assuming
%          uniform maximum intensity
%'Inputs': 'Phi_curr', Current phase estimate;
%          'dIdz_curr', Current intensity derivative estimate;
%'         'dz',Propagation distance;
%          'I0', Infoucs intensity image
%          'pixelsize'��pixelsize('mm');
%          'k', Wavenumber
%          'RegPara', regularzation parameter (to remove low-frequency noise)
%'Outputs': 'phi',object phase map;
%           'dIdz_est', estimated dIdz;
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
% last modified on 02/27/2020
% by Chao Zuo (zuochao@njust.edu.cn,surpasszuo@163.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi, dIdz_est] = TIE_DCT_solver(Phi_curr,dIdz_curr,dz,I0,Pixelsize,k,RegPara,Method)

% Fast DCT based on mirror padding and FFT
dIdz_double = EvenFlip(dIdz_curr);
Phi_curr_double = EvenFlip(Phi_curr);
I0_double = EvenFlip(I0);
Lambda = 2*pi/k;

J = -k*dIdz_double;
[Ny,Nx] = size(J);
Fx = (-(Nx/2):1:(Nx/2-1))/Nx/Pixelsize; % frequency coordinate X.
Fy = (-(Ny/2):1:(Ny/2-1))/Ny/Pixelsize; % frequency coordinate Y.
[U,V] = meshgrid(Fx,Fy);

% Coordinates in frequency domain
fx = fftshift(U);
fy = fftshift(V);

% Terms for derivative calculation in frequency domain
Cx = 2*1i*pi*fx;
Cy = 2*1i*pi*fy;

% Calculate Psi by TIE with non-uniform intensity.
FJ = fft2(J);
Fpsi = FJ.*(Cx.*Cx+Cy.*Cy)./(RegPara/Pixelsize.^4+(Cx.*Cx+Cy.*Cy).^2);

% Calculate Phi by TIE with non-uniform intensity
Fdpsidx = Fpsi.*Cx;
Fdpsidy = Fpsi.*Cy;
dpsidx = real(ifft2(Fdpsidx));
dpsidy = real(ifft2(Fdpsidy));

% Intensity thresholding to prevent divide-by-zero error
% I0_double(find(I0_double<0.01*max(max(I0_double))))=inf;
I0_double(find(I0_double<0.01*max(max(I0_double))))=mean2(I0_double);


% Phase derivative
dPhidx = dpsidx./I0_double;
dPhidy = dpsidy./I0_double;

FdPsidx = fft2(dPhidx);
FdPsidy = fft2(dPhidy);

% Phase Laplacian
Fd2Phidx2 = FdPsidx.*Cx;
Fd2Phidy2 = FdPsidy.*Cy;

d2Phidx2 = real(ifft2(Fd2Phidx2));
d2Phidy2 = real(ifft2(Fd2Phidy2));
d2Phi = d2Phidx2 + d2Phidy2;

% Phase integration
Fd2Phi = fft2(d2Phi);
Fphi = Fd2Phi.*(Cx.*Cx+Cy.*Cy)./(RegPara/Pixelsize.^4+(Cx.*Cx+Cy.*Cy).^2);
Phi_double = real(ifft2(Fphi));

% Crop to orignal size
Phi = Phi_double(1:end/2,1:end/2);

if(strcmp(Method,'TIE'))
    % Calculate the estimate of dIdz using TIE
    Fdphidx = Fphi.*Cx;
    Fdphidy = Fphi.*Cy;
    dphidx = real(ifft2(Fdphidx));
    dphidy = real(ifft2(Fdphidy));
    dPsidx = I0_double.*dphidx;
    dPsidy = I0_double.*dphidy;
    FdPsidx = fft2(dPsidx);
    FdPsidy = fft2(dPsidy);
    Fd2Psidx2 = FdPsidx.*Cx;
    Fd2Psidy2 = FdPsidy.*Cy;
    d2Psidx2 = real(ifft2(Fd2Psidx2));
    d2Psidy2 = real(ifft2(Fd2Psidy2));
    laplacePsi = d2Psidx2+d2Psidy2;
    % Estimate dIdz
    dIdz_est = laplacePsi./(-k);
    dIdz_est = dIdz_est(1:end/2,1:end/2);
    
elseif(strcmp(Method,'Angular Spectrum'))
    % Calculate the estimate of dIdz using Angular Spectrum
    FU=fftshift(fft2(sqrt(I0_double).*exp(1i*(Phi_double+Phi_curr_double))));
    Exp_term = sqrt(1-(Lambda*U).^2-(Lambda*V).^2);
    H = exp(1i*k*dz*Exp_term);
    H((1-(Lambda*U).^2-(Lambda*V).^2)<0) = 0; % neglect evanescent wave
    Uz=ifft2(ifftshift(FU.*H));
    Iz=Uz.*conj(Uz);
    dIdz_est = (Iz-I0_double)/dz;
    dIdz_est = dIdz_est(1:end/2,1:end/2);

    
elseif(strcmp(Method,'Fresnel'))
    % Calculate the estimate of dIdz using Fresnel method
    FU=fftshift(fft2(sqrt(I0_double).*exp(1i*(Phi_double+Phi_curr_double))));
    H=exp(1i*k*dz*(1-((Lambda*U).^2+(Lambda*V).^2)/2));
    Uz=ifft2(ifftshift(FU.*H));
    Iz=Uz.*conj(Uz);
    dIdz_est = (Iz-I0_double)/dz;
    dIdz_est = dIdz_est(1:end/2,1:end/2);
else
    errordlg('Type of transfer function must be <Angular Spectrum> or <Fresnel>','Error');
end

end


% Subfunction for fast DCT based on mirror padding and FFT
function AA = EvenFlip(A)

temp = [A,fliplr(A)];
AA = [temp',(flipud(temp))']';

end
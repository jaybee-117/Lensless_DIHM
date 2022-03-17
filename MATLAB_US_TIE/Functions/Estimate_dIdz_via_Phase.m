%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate intensity derivative via current phase
% last modified on 02/27/2020
% by Chao Zuo (zuochao@njust.edu.cn,surpasszuo@163.com)
% _________________________________________________________________________
% dIdz_est = Estimate_dIdz_via_Phase(Phi,I0,Pixelsize,k,Method)
% Inputs :
%	Phi: In-focus Phase;
%	I0: In-focus intensity;
%	Pixelsize: Pixelsize (in the sample space);
%   dz: Propagation distance (used for 'Angular Spectrum' or 'Fresnel')
%   k: Wave-number;
%   Method: type of method used ('Angular Spectrum' or 'Fresnel' or 'TIE').
% Output:
%   dIdz_est: Estimated intensity derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dIdz_est = Estimate_dIdz_via_Phase(Phi,I0,Pixelsize,dz,Lambda,Method)

[Ny,Nx]=size(Phi);

k=2*pi/Lambda; % wave number

Fx = (-(Nx/2):1:(Nx/2-1))/Nx/Pixelsize; % frequency coordinate X.
Fy = (-(Ny/2):1:(Ny/2-1))/Ny/Pixelsize; % frequency coordinate Y.

[U,V] = meshgrid(Fx,Fy);


if(strcmp(method,'Angular Spectrum'))
    H=exp(1i*k*z*sqrt(1-uu.^2-vv.^2)); % Angular Specturm method
elseif(strcmp(method,'Fresnel'))
    H=exp(1i*k*z*(1-(uu.^2+vv.^2)/2)); % Fresnel method
    
    
elseif(strcmp(method,'TIE'))  % TIE method
    fx = fftshift(U);
    fy = fftshift(V);
    
    Cx = 2*1i*pi*fx;
    Cy = 2*1i*pi*fy;
    
    Fphi = fft2(Phi);
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
    
    dIdz_est = laplacePsi./(-k);
    
else
    errordlg('Type of Method must be <Angular Spectrum> or <Fresnel> or <TIE>','Error');
end

end

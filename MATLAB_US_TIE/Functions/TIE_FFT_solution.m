%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function for FFT-TIE solver.
% Version 2.0 -
% Related Reference:
% D. Paganin and K. A. Nugent, Phys. Rev. Lett. 80, 2586 (1998).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Phi,RMSE,Time] = TIE_FFT_solution(dIdz,I0,Pixelsize,k,RegPara,IntThr,True_phase)
% Purpose: Solution to TIE using FFT-based Possion solver
% inputs :
%	dIdz: Intensity derivative dI/dz;
%   dz: propagation distance;
%   I0: In-focus intensity;
%	Pixelsize: Pixelsize (in the sample space);
%   RegPara: Regularization parameter (to remove low-frequency noise);
%   k: Wave-number;
%   Iter: Maximum iteration times;
%   Method: type of method used ('FFT' or 'Direct');
%   IntThr, Intensity threshold to prevent divide-by-zero error
%   True_phase�� Ture phase.
% Outputs:
%   Phi: Phase recovered;
%   RMSE: Root mean squre error;
%   Time: Computation time;
%��������������������������������������������������������������������������������
% last modified on 02/27/2020
% by Chao Zuo (zuochao@njust.edu.cn,surpasszuo@163.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [Phi,RMSE,Time] = TIE_FFT_solution(dIdz,I0,Pixelsize,k,RegPara,IntThr,True_phase)

t0 = clock;
J = -k*dIdz;
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

% Calculate Psi by TIE with non-uniform intensity.
FJ = fft2(J);
Fpsi = FJ.*(Cx.*Cx+Cy.*Cy)./(RegPara/Pixelsize.^4+(Cx.*Cx+Cy.*Cy).^2);

% Calculate Phi by TIE with non-uniform intensity
Fdpsidx = Fpsi.*Cx;
Fdpsidy = Fpsi.*Cy;
dpsidx = real(ifft2(Fdpsidx));
dpsidy = real(ifft2(Fdpsidy));

% Intensity thresholding to prevent divide-by-zero error
I0(find(I0<IntThr*max(max(I0))))...
    =IntThr*max(max(I0));

% Phase derivative
dPhidx = dpsidx./I0;
dPhidy = dpsidy./I0;

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
Phi = real(ifft2(Fphi));

% Computation time
t1 = clock;
t2 = etime(t1,t0);
Time = t2;

% Calculate RMSE
if exist('True_phase')
    err = Phi-True_phase;
    err = err-mean(err(~isnan(err)));
    RMSE=sqrt(sum(sum(err.^2,'includenan'),'includenan')/(numel(err)-numel(find(isnan(err)))));
end

end

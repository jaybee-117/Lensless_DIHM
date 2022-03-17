%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function for FFT-TIE solver-.
% Version 3.0 -
% Related Reference:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  C. Zuo, Q. Chen, and A. Asundi,
% "Boundary-artifact-free phase retrieval with the transport of intensity
%  equation: fast solution with use of discrete cosine transform,"
% Optics Express 22, 9220-9244 (2014).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Phi,RMSE,Time] = TIE_DCT_solution(dIdz,I0,Pixelsize,k,RegPara,IntThr,True_phase)
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
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
% Version 3.0 -
% Coded by Chao Zuo - 2013-11-19
% Lastest edited by Chao Zuo - 2014-3-9
% Last modified on 02/27/2020
% by Chao Zuo (zuochao@njust.edu.cn,surpasszuo@163.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [Phi,RMSE,Time] = TIE_DCT_solution(dIdz,I0,Pixelsize,k,RegPara,IntThr,True_phase)

t0 = clock;

% find the apecture region
Aperture = false(size(I0));
bw = im2bw(I0,0.5*graythresh(I0));
stats = [regionprops(bw); regionprops(not(bw))];
rect = round(stats(1).BoundingBox);

% Determine the BC_Region
if([rect(4),rect(3)] == size(I0))
    Aperture = true(size(I0));
    BC_Region = 0;
else
    for BC_Region = 1:20
        if( (rect(2)-BC_Region)<1 || (rect(2)+rect(4)-1)+BC_Region> size(I0,1)...
         ||(rect(1)-BC_Region)<1 || (rect(1)+rect(3)-1)+BC_Region> size(I0,2))
            BC_Region = BC_Region-1;
            Aperture(rect(2)-BC_Region:(rect(2)+rect(4)-1)+BC_Region,...
            rect(1)-BC_Region:(rect(1)+rect(3)-1)+BC_Region) = 1;
            break;
        end
        Aperture(rect(2)-BC_Region:(rect(2)+rect(4)-1)+BC_Region,...
            rect(1)-BC_Region:(rect(1)+rect(3)-1)+BC_Region) = 1;
        if(max(max(dIdz(~Aperture)))<max(max(dIdz(Aperture)))/1e3)
            break;
        end
    end
end

% crop the image
dIdz = reshape(dIdz(Aperture),[rect(4)+2*BC_Region,rect(3)+2*BC_Region]);
I0 = reshape(I0(Aperture),[rect(4)+2*BC_Region,rect(3)+2*BC_Region]);

% Fast DCT based on mirror padding and FFT
dIdz = EvenFlip(dIdz);

J = -k*dIdz;
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
Psi=real(ifft2(Fpsi));

% Intensity thresholding to prevent divide-by-zero error
I0(find(I0<IntThr*max(max(I0))))=IntThr*max(max(I0));

if(BC_Region~=0)
    % Exclude BC region
    Psi = Psi(1:end/2,1:end/2);
    Psi = Psi(1+BC_Region:end-BC_Region,1+BC_Region:end-BC_Region);
    Psi = EvenFlip(Psi);
    Fpsi = fft2(Psi);
    
    % Phase integration on the aperture region (exclude the BC region)
    I0 = I0(1+BC_Region:end-BC_Region,1+BC_Region:end-BC_Region);
    Aperture = false(size(Aperture));
    Aperture(rect(2):(rect(2)+rect(4)-1),rect(1):(rect(1)+rect(3)-1)) = 1;
    
    % Fast DCT based on mirror padding and FFT
    I0 = EvenFlip(I0);
    [Ny,Nx] = size(I0);
    Fx = (-(Nx/2):1:(Nx/2-1))/Nx/Pixelsize; % frequency coordinate X.
    Fy = (-(Ny/2):1:(Ny/2-1))/Ny/Pixelsize; % frequency coordinate Y.
    [U,V] = meshgrid(Fx,Fy);
    
    % Coordinates in frequency domain
    fx = fftshift(U);
    fy = fftshift(V);
    
    % Terms for derivative calculation in frequency domain
    Cx = 2*1i*pi*fx;
    Cy = 2*1i*pi*fy;
    
    % Calculate Phi by TIE with non-uniform intensity
    Fdpsidx = Fpsi.*Cx;
    Fdpsidy = Fpsi.*Cy;
    dpsidx = real(ifft2(Fdpsidx));
    dpsidy = real(ifft2(Fdpsidy));
    
else
    % Calculate Phi by TIE with non-uniform intensity
    Fdpsidx = Fpsi.*Cx;
    Fdpsidy = Fpsi.*Cy;
    dpsidx = real(ifft2(Fdpsidx));
    dpsidy = real(ifft2(Fdpsidy));
    
    I0 = EvenFlip(I0);
end

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
Phi_ROI = real(ifft2(Fphi));

% Crop to orignal size
Phi_ROI = Phi_ROI(1:end/2,1:end/2);

% Put back to the original position
Phi = zeros(size(Aperture));
Phi(~Aperture) = NaN;
Phi(Aperture) = Phi_ROI;
Phi = reshape(Phi,size(Aperture));

% Computation time
t1 = clock;
t2 = etime(t1,t0);
Time = t2;

% Calculate RMSE
if exist('True_phase')
    err = Phi-True_phase;
    err = err-mean(err(~isnan(err)));
    RMSE=sqrt(nansum(nansum(err.^2))/(numel(err)-numel(find(isnan(err)))));
end

end


% Subfunction for fast DCT based on mirror padding and FFT
function AA = EvenFlip(A)

temp = [A,fliplr(A)];
AA = [temp',(flipud(temp))']';

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file for Simulation 2 for Universal Solution to the TIE.
% Version 2.0 -
% Related Reference:
% "On a universal solution to the transport-of-intensity equation"
% Jialin Zhang, Qian Chen, Jiasong Sun, Long Tian, Chao Zuo
%
% last modified on 02/27/2020
% by Chao Zuo (zuochao@njust.edu.cn,surpasszuo@163.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clear everything existing.
clc
close all
clearvars
addpath('.\Functions');
addpath('.\Data')

%% Defined all the parameters
% Set the pixel resolution of image.
Nx = 256;
Ny = Nx;

Pixelsize = 2e-6;   % pixel size (m)
Lambda = 0.633e-6;  % wavelength (m)
k = 2*pi/Lambda;    % wave number
dz = 1e-6;          % defocus distance (m)

Flag_Intensity = 0;	% 1: uniform intensity
                    % 0: nonuniform intensity distribution

Flag_Hann = 0;      % 1: Use Hann distribution phase object.
                    % 0: Do not use Hann phase object.

Flag_AperShape = 'circle';
% 'rect': Rect.
% 'circle': R=0.65*Rmax Aperture.
% 'circle-': circle-.
% 'ellipse': ellipse.
% 'ellipse-': ellipse-.
% 'telescope-': telescope-.
% 'telescope': telescope.
% 'full': Do not use Aperture.

sigma = 1.1e-4;    % sigma of Gaussian intensity distribution (m)

% Symbolize expression.
syms Phase_sym Inten_sym x y
Inten_sym = exp((-x.^2-y.^2)/(2*sigma.^2));
Phase_sym = 0.82+(-0.7*(x*1e3)+(x*1e3).^2*10-(y*1e3).^2*10+2*(y*1e3));

%% Generate intensity and phase distribution.
eval Generate_Optical_Field;

%% Numerical propagation to generate defocused images.
% Complex amplitude of optical field in foucs.
U0z = sqrt(I0).*exp(1i.*Phi);

% Cut the centrial region which is able to be captured by CCD.
Indx = Nx_LF/4+1:Nx_LF/4+Nx;
Indy = Ny_LF/4+1:Ny_LF/4+Ny;
I0 = I0(Indy,Indx);
U0z = U0z(Indy,Indx);
Phi = Phi(Indy,Indx);

% Numerical propagation for defocused images.
Method = 'Angular Spectrum';
[Uz,Iz] = Numerical_Propagation(U0z,dz,Pixelsize,Lambda,Method);
Aperture = double(Aperture(Indy,Indx));
Phi(~Aperture) = NaN;
Aperture(Aperture==0)=NaN;

% Axial intensity derivative.
dIdz = (Iz-I0)/dz;

% Check if satisfying Energy Conservation .
disp(['mean2(Iz-I0)) = ' num2str(mean2(Iz-I0))]);
if abs(mean2(Iz-I0))>max(max(I0))/255    % Need to check the threshold 
         warndlg('Conservation of Energy may not be satisfied!');
end

%% Show the images.
figure
subplot(1,3,1)
imshow(I0,[0,1]);
title('In-focus intensity');
subplot(1,3,2)
imshow(Aperture,[]);
title('Aperture');
subplot(1,3,3)
imshow(Phi,[0,1.5]);
title('Phase');
colormap(gca,jet);

figure
subplot(1,3,1)
imshow(I0,[]);
title('In-focus intensity');
subplot(1,3,2)
imshow(Iz,[]);
title('Defocused intensity');
subplot(1,3,3)
imshow(dIdz,[]);
title('Intensity Derivative');


%% Solve TIE with iteration parameters.
RegPara = eps;
IntThr = 0.01;
JudgeFlag_DCT = 1;
JudgeFlag_US = 1;
MaxIterNum = 500;
%% Solve TIE with FFT-TIE.
[Phi_FFT,RMSE_FFT,Time_FFT] = TIE_FFT_solution(dIdz,I0,Pixelsize,k,RegPara,IntThr,Phi);
% Error Calculation.
phiErr_FFT = Remove_Piston(Phi_FFT - Phi);

%% Solve TIE with IterDCT-TIE.
[Phi_IterDCT,RMSE_IterDCT,Time_IterDCT] = Iter_DCT_solution(dIdz,dz,I0,Pixelsize,k,RegPara,MaxIterNum,JudgeFlag_DCT,Method,Phi);
% Error Calculation.
phiErr_IterDCT = Remove_Piston(Phi_IterDCT - Phi);

%% Solve TIE with US-TIE.
[Phi_US,RMSE_US,Time_US] = Universal_Solution(dIdz,dz,I0,Pixelsize,k,RegPara,MaxIterNum,JudgeFlag_US,Method,Phi);
% Error Calculation.
phiErr_US = Remove_Piston(Phi_US - Phi);

%% Show Reconstruction Result
figure;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',scrsz);

subplot(2,3,1)
imshow(Phi_FFT.*Aperture,[]);
colormap(gca,jet);colorbar
title('Reconstructed Phase (FFT-TIE)');
subplot(2,3,2)
imshow(Phi_IterDCT.*Aperture,[]);
colormap(gca,jet);colorbar
title('Reconstructed Phase (IterDCT-TIE)');
subplot(2,3,3)
imshow(Phi_US.*Aperture,[]);
colormap(gca,jet);colorbar
title('Reconstructed Phase (US-TIE)');

subplot(2,3,4)
imshow(phiErr_FFT,[ ]);
colormap(gca,jet);colorbar
title(['Error (FFT-TIE) ','RMSE =',num2str(RMSE_FFT(end))]);
subplot(2,3,5)
imshow(phiErr_IterDCT,[ ]);
colormap(gca,jet);colorbar
title(['Error (IterDCT-TIE) ','RMSE = ', num2str(RMSE_IterDCT(end))]);
subplot(2,3,6)
imshow(phiErr_US,[-0.05,0.05]);
colormap(gca,jet);colorbar
title(['Error (US-TIE) ','RMSE = ', num2str(RMSE_US(end))]);

figure
plot(0:size(RMSE_US,2)-1,RMSE_US,'color',[0.8,0.2,0.2],'LineWidth',3);
hold on
plot(0:size(RMSE_IterDCT,2)-1,RMSE_IterDCT,'color',[0.2,0.2,0.8],'LineWidth',3);
title('RMSE Comparison')
ylim([0,0.3]);
legend('US-TIE','IterDCT');
xlabel('Iteration NO.')

%% Computation Time
figure
subplot(1,2,1)
plot(Time_IterDCT,'color',[0.8,0.2,0.2],'LineWidth',3);
hold on
plot(Time_US,'color',[0.2,0.2,0.8],'LineWidth',3);
title('Time per iteration (s)')
legend('IterDCT','US-TIE');
xlabel('Iteration NO.')
subplot(1,2,2)
TimeUS_plot=cumsum(Time_US);
Time_IterDCT_plot=cumsum(Time_IterDCT);
plot(Time_IterDCT_plot,'color',[0.8,0.2,0.2],'LineWidth',3);
hold on
plot(TimeUS_plot,'color',[0.2,0.2,0.8],'LineWidth',3);
title('Overall time (s)')
legend('IterDCT','US-TIE');
xlabel('Iteration NO.')
disp(['FFT-TIE computation time:', num2str(Time_FFT(end)),'s']);
disp(['IterDCT computation time:', num2str(Time_IterDCT_plot(end)),'s']);
disp(['US--TIE computation time:', num2str(TimeUS_plot(end)),'s']);



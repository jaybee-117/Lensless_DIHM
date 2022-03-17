%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation using DIHM module in ImageJ and phase retrieval using TIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear everything existing.
clc
close all
clearvars
addpath('Functions/');
addpath('Data/Microlens');
addpath('DIHM_ImageJ/');
rmpath('Data/HeLaCell');

%% Setting up system parameters
I0 = im2double(imread('I_300.tif'));
Iz = im2double(imread('I_299.tif'));

Mag = 20;                  % Magnification
Pixelsize = 4.65e-6./Mag;  % Pixelszie
lambda = 550e-9;           % Wavelength (m)
k = 2*pi/lambda;           % Wave number
dz = -1e-6;       	   % Defocus distance (m)

% Axial intensity derivative
dIdz = (Iz-I0)/(dz);

% Valid domain within the apecture
Aperture = ones(size(I0));
Aperture(I0<max(max(I0))/10)=0;

% Check if satisfying Energy Conservation
disp(['mean2(Iz-I0)) = ' num2str(mean2(Iz-I0))]);
if abs(mean2(Iz-I0))>max(max(I0))/255    % Need to check the threshold 
         warndlg('Conservation of Energy may not be satisfied!');
end

%% Show the images.
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

%% Set parameters.
RegPara = 5e-9; % remove low-freq noise and aberration 
IntThr = 0.01;
MaxIterNum = 30; % maximum iteration
JudgeFlag_DCT = 1;
JudgeFlag_US = 1;

%% Solve TIE with FFT-TIE
phi_FFT = TIE_FFT_solution(dIdz,I0,Pixelsize,k,RegPara,IntThr);

%% Solve TIE with Iter-DCT-TIE
phi_DCT = Iter_DCT_solution(dIdz,dz,I0,Pixelsize,k,RegPara,MaxIterNum,JudgeFlag_DCT,'Angular Spectrum');

%% Solve TIE with US-TIE
phi_US = Universal_Solution(dIdz,dz,I0,Pixelsize,k,RegPara,MaxIterNum,JudgeFlag_US,'Angular Spectrum');

%% Show Reconstruction Result
figure
subplot(1,3,1);
imshow(phi_FFT.*Aperture,[]);
title('FFT-TIE');
colorbar;
subplot(1,3,2);
imshow(phi_DCT.*Aperture,[]);
title('Iter-DCT-TIE');
colorbar;
subplot(1,3,3);
imshow(phi_US.*Aperture,[]);
title('US-TIE');
colorbar;




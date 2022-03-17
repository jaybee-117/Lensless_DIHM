%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file for Experiment 2 for Universal Solution to the TIE.
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
addpath('.\Data\HeLaCell')
rmpath('.\Data\Microlens')

%% Setting up system parameters
I0 = double(imread('I0.tif'));
Iz = double(imread('Iz.tif'));

Mag = 20;                  % Magnification
Pixelsize = 4.65e-6./Mag;  % Pixelszie
lambda = 550e-9;           % Wavelength (m)
k = 2*pi/lambda;           % Wave number
dz = 2.5e-6;       	   % Defocus distance (m)

% Axial intensity derivative
dIdz = (Iz-I0)/(dz);

% Valid domain within the apecture
Aperture = ones(size(I0));
Aperture(I0<max(max(I0))/10)=NaN;

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
MaxIterNum = 20; % maximum iteration
JudgeFlag_DCT = 1; % avoid overshooting during the later stage of iteration 
JudgeFlag_US = 1;
%% Solve TIE with FFT-TIE
phi_FFT = TIE_FFT_solution(dIdz,I0,Pixelsize,k,RegPara,IntThr);

%% Solve TIE with Iter-DCT-TIE
phi_DCT = Iter_DCT_solution(dIdz,dz,I0,Pixelsize,k,RegPara,MaxIterNum,JudgeFlag_DCT,'Angular Spectrum');

%% Solve TIE with US-TIE
phi_US = Universal_Solution(dIdz,dz,I0,Pixelsize,k,RegPara,MaxIterNum,JudgeFlag_US,'TIE');

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
imshow(phi_US.*Aperture,[-1.5,2]);
title('US-TIE');
colorbar;



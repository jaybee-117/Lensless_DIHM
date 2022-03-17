%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function for universal solution to TIE.
% Version 2.0 -
% Related Reference:
% "On a universal solution to the transport-of-intensity equation"
% Jialin Zhang, Qian Chen, Jiasong Sun, Long Tian, Chao Zuo
%
% last modified on 02/27/2020
% by Chao Zuo (zuochao@njust.edu.cn,surpasszuo@163.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _________________________________________________________________________
% [Phi,RMSE,Time]  = Universal_Solution(dIdz,dz,I0,Pixelsize,k,r,Iter, Method, True_phase)
% Inputs :
%	dIdz: Intensity derivative dI/dz;
%   dz: propagation distance;
%   I0: In-focus intensity;
%	Pixelsize: Pixelsize (in the sample space);
%   r: Regularization parameter (to remove low-frequency noise);
%   k: Wave-number;
%   Iter: Maximum iteration times;
%   JudgeFlag: Check for convergence (1-Yes 0-No)
%   Method: type of method used ('FFT' or 'Direct');
%   Ture_phase�� Ture phase.
% Outputs:
%   Phi: Phase recovered;
%   RMSE: Root mean squre error;
%   Time: Computation time;
% _________________________________________________________________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi,RMSE,Time] = Universal_Solution(dIdz,dz,I0,Pixelsize,k,r,Iter, JudgeFlag, Method, True_phase)

%% Check nargin
if strcmpi(Method,'TIE')
elseif strcmpi(Method,'Angular Spectrum')
elseif strcmpi(Method,'Fresnel')
else
    error('Currently ''Method'' only allow ''TIE'' or ''Angular Spectrum'' or ''Fresnel''.');
end

%% Initalize dIdz and phinow
RMSE=[];
Time=[];

dIdz_curr = dIdz;
dIdz_temp = dIdz;

Phi_curr = zeros(size(dIdz));

if exist('True_phase')
    err = Phi_curr-True_phase;
    err = err-mean(err(~isnan(err)));
    RMSE_Temp = sqrt(sum(sum(err.^2,'includenan'),'includenan')/(numel(err)-numel(find(isnan(err)))));
    RMSE = [RMSE RMSE_Temp];
end


%% iterative US TIE solver
for n = 1:Iter
    
    % Estimate Phase with TIE_MAX_solver
    t0 = clock;
    [Phi_est, dIdz_est] = TIE_MAX_solver(Phi_curr,dIdz_curr,dz,I0,Pixelsize,k,r,Method);
    t1 = clock;
    t2 = etime(t1,t0);
    dIdz_temp = dIdz_curr;
    % Calculate dIdz Error
    if strcmpi(Method,'TIE')
        dIdz_curr = dIdz_curr - dIdz_est;
    elseif strcmpi(Method,'Fresnel') || strcmpi(Method,'Angular Spectrum')
        dIdz_curr = dIdz - dIdz_est;
    end
    
    % Current Phase Estimate
    Phi_curr = Phi_curr + Phi_est;
    
    % Calculate RMSE
    if exist('True_phase')
        err = Phi_curr-True_phase;
        err = err-mean(err(~isnan(err)));
        RMSE_Temp=sqrt(sum(sum(err.^2,'includenan'),'includenan')/(numel(err)-numel(find(isnan(err)))));
        RMSE = [RMSE RMSE_Temp];
        Time = [Time t2];
    end
    
    % Check Convergence
    if(JudgeFlag)
        if(max(max(dIdz_curr))<max(max(dIdz))*1e-3 || max(max(dIdz_curr))>1.05*max(max(dIdz_temp)) )
            disp(['US-TIE iteration time:', int2str(n)]);
            break
        end
    end
    if(n == Iter)
        disp(['US-TIE iteration time:', int2str(Iter)]);
    end
end
Phi = Phi_curr;
end
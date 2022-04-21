%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [Phi_Retrieved]  = Phase_Retrieve(U_Object,lambda_min,d_lambda,lambda_max)
% Inputs :
%	U_Object: Wavefront at the object;
%	lambda_min: Starting Laser Wavelength;
%   d_lambda: Difference between two consecutive wavelengths;
%   lambda_max: Ending Laser Wavelength;
% Outputs:
%   Phi_Retrieved: Retrieved Phase;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Phi_Retrieved] = Phase_Retrieve(U_Object,lambda_min,d_lambda,lambda_max)

%%%%%%%%%%%%%%%%%%%%%%% Generating Recordings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Obj_Cam_Dist = 300;

i = 1;
for x = lambda_min:d_lambda:lambda_max
    [~,temp,~] = Numerical_Front_Propagation(U_Object,Obj_Cam_Dist,x);
    if (x == lambda_min)
        I = temp;
    else
        I = cat(3,I,temp);
    end
%     i = i+1;
%     figure (i)
%     imagesc(temp);
end
[a,b,c] = size(I);
fprintf("Size of intensity stack: %dx%dx%d", a,b,c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Phase Retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi_t = I(:,:,1)./256 * 2 * pi;
% phi_q = I(:,:,2)./256 * 2 * pi;
% phi_t = zeros(256,256);
% phi_q = zeros(256,256);
phi_t = rand(256,256);
phi_q = rand(256,256);
MSE = 1000;
Minimum_error = 1000;
it = 10;
MSE_trend = 0;
a = 1;
b = 1;
next_lambda = lambda_min;
h = zeros(256,256);
while MSE > 0 && it > 0
    lambda = lambda_min;
    it = it - 1;
    for t = 1:c
            q = mod(t,c) + 1;
            %% Calculate Wavefront
            U_t = sqrt(I(:,:,t)).*exp(1i*phi_t);
            U_q = sqrt(I(:,:,q)).*exp(1i*phi_q);
            %% Back Propagate
            [U_back_t,I_back_t,Phi_back_t] = Numerical_Back_Propagation(U_t,Obj_Cam_Dist,lambda);
            next_lambda = lambda + d_lambda;
            if next_lambda > lambda_max
                next_lambda = lambda_min;
            end
            [U_back_q,I_back_q,Phi_back_q] = Numerical_Back_Propagation(U_q,Obj_Cam_Dist,next_lambda);
            %% Calculate Object Profile
            synth_lambda = lambda * next_lambda / abs(lambda - next_lambda);
            h1 = h;
            h = (synth_lambda/(1.5 - 1.0)) .* abs(Phi_back_t - Phi_back_q)/ ( 2* pi) * 1e-9;
            h2 = h;
            h_error = abs(sqrt(h2)-sqrt(h1))^2;
            H_Error = sum(h_error(:))/numel(sqrt(h1));
            next_next_lambda = next_lambda + d_lambda;
            if next_next_lambda > lambda_max
                next_next_lambda = lambda_min;
            end
            %% Calculate Wavefront at Object Plane
            U_tt = (sqrt(I_back_t) + sqrt(I_back_q))/2 .* exp(2*pi*i*h/next_lambda);
            U_qq = (sqrt(I_back_t) + sqrt(I_back_q))/2 .* exp(2*pi*i*h/next_next_lambda);
            %% Front Propagation
            [U_front_t,I_front_t,phi_t] = Numerical_Front_Propagation(U_tt,Obj_Cam_Dist,next_lambda);
            [U_front_q,I_front_q,phi_q] = Numerical_Front_Propagation(U_qq,Obj_Cam_Dist,next_next_lambda);
            %% Error Calculation
            D = abs((sqrt(I_front_t) - sqrt(I(:,:,q))))^2;
            MSE = sum(D(:))/numel(sqrt(I_front_t));
%             fprintf("H_Error = %d\n",H_Error);
%             fprintf("MSE = %d\n",MSE);
            if Minimum_error > MSE
                Minimum_error = MSE;
                Phi_Retrieved = Phi_back_t;
            end
            MSE_trend = [MSE_trend; MSE];
            lambda = next_lambda;
    end
end
figure (10)
plot(MSE_trend);
figure(11)
imagesc(h);
colormap('hot');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
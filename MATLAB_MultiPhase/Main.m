%%%%%%%%%%%%%%%%%%%%  Parameters  %%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
Obj_Cam_Dist = 3;
r = 10; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  Test Object  %%%%%%%%%%%%%%%%%%%%%%%%
% x = 0:255;
% y = 0:255;
% [X,Y] = meshgrid(x,y);
% I = zeros(length(x),length(y));
% Phi = zeros(length(x),length(y));
% for i = 1:256
%     for j = 1:256
%         if ((i-64)^2 + (j-64)^2 <= r^2)
%             I(i,j) = 255;
%         end
%     end
% end
% 
% Phi = imgaussfilt(I,5);
% figure ( 10)
% mesh(X,Y,Phi);
% U = sqrt(I) .*exp(Phi*1i);

I = double(rgb2gray(imread("Dots.png")));
% I = double(imread("Dots.png"));
Phi = (I/256)*(pi - 1); %normalising phase

I_add = zeros(256,256);
Phi_add = zeros(256,256);
for i = 1:256
    for j = 1:256
        if (Phi(i,j) == 0)
            I_add(i,j) = (2*128^2 - (i-128)^2 - (j-128)^2)* (1/(2*128^2)) * 128;
            Phi_add(i,j) = (2*128^2 - (i-128)^2 - (j-128)^2)* (1/(2*128^2)) * 2;
        end
    end
end

I = imgaussfilt(I + I_add,2);
Phi = imgaussfilt(Phi + Phi_add,2);
[a,b] = size(I);
x = 0:a-1;
y = 0:b-1;

[X,Y] = meshgrid(x,y);
U = sqrt(I) .*exp(Phi*1i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  Display the Object  %%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
% mesh(X,Y,abs(U));
imagesc(abs(U).^2);
colorbar;
title("Object Intensity");
ax = gcf;
exportgraphics(ax,'Object_Intensity.png','Resolution',300);

figure (2)
% mesh(X,Y,angle(U));
imagesc(angle(U));
colorbar;
title("Object Phase");
ax = gcf;
exportgraphics(ax,'Object_Phase.png','Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  Generate propagated wavefront  %%%%%%%%%%%%%%%%%%%%%%
[U_n,I_n,Phi_n] = Numerical_Front_Propagation(U,Obj_Cam_Dist,840);

figure(3)
% mesh(X,Y,I_n);
imagesc(I_n);
colorbar;
title("Recorded Intensity pattern");
ax = gcf;
exportgraphics(ax,'Recorded_Intensity.png','Resolution',300);

figure(4)
% mesh(X,Y,Phi_n);
imagesc(Phi_n);
colorbar;
title("Recorded Phase patern");
ax = gcf;
exportgraphics(ax,'Recorded_Phase.png','Resolution',300);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  Back propagated wavefront  %%%%%%%%%%%%%%%%%%%%%%%%
[U_back,I_back,Phi_back] = Numerical_Back_Propagation(U_n,Obj_Cam_Dist,840);

figure(5)
% mesh(X,Y,I_back);
imagesc(I_back);
colorbar;
title("Back propagated Intensity");
ax = gcf;
exportgraphics(ax,'Back_Propagated_Intensity.png','Resolution',300);

figure(6)
% mesh(X,Y,Phi_back);
imagesc(Phi_back);
colorbar;
title("Back propagated Phase");
ax = gcf;
exportgraphics(ax,'Back_Propagated_Phase.png','Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  Phase Retrieved Wavefront  %%%%%%%%%%%%%%%%%%%%%%%%
[~,Retrieved_Intensity,Retrieved_Phase] = Phase_Retrieve(U,800,5,850);

figure(7)
% mesh(X,Y,Retrieved_Phase);
imagesc(Retrieved_Phase);
colorbar;
title("Retrieved Phase");
ax = gcf;
exportgraphics(ax,'Retrieved_Phase.png','Resolution',300);

figure(8)
% mesh(X,Y,Retrieved_Phase);
imagesc(Retrieved_Intensity);
colorbar;
title("Reconstructed Intensity");
ax = gcf;
exportgraphics(ax,'Retrieved_Intensity.png','Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
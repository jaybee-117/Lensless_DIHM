%%%%%%%%%%%%%%%%%%%%  Parameters  %%%%%%%%%%%%%%%%%%%%%%%%

Obj_Cam_Dist = 300;
r = 10; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  Test Object  %%%%%%%%%%%%%%%%%%%%%%%%
x = 0:255;
y = 0:255;

[X,Y] = meshgrid(x,y);

I = zeros(length(x),length(y));
Phi = zeros(length(x),length(y));
for i = 0:256
    for j = 0:256
        if ((i-128)^2 + (j-128)^2 <= r^2)
            I(i,j) = 255;
            Phi(i,j) = 2*(pi/r^2)*sqrt(r^2 - (i-128)^2 - (j-128)^2);
        end
    end
end
U = I .* exp(1i*Phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  Display the Object  %%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(2,1,1)
mesh(x,y,I);
title("Object Intensity");
imwrite(I,"Object_Intensity.png");

subplot(2,1,2)
mesh(x,y,Phi);
title("Object Phase");
imwrite(Phi,"Object_Phase.png");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  Generate propagated wavefront  %%%%%%%%%%%%%%%%%%%%%%
[U_n,I_n,Phi_n] = Numerical_Propagation(U,Obj_Cam_Dist,840);

figure(3)
mesh(x,y,I_n);
title("Recorded Intensity pattern");
imwrite(I_n,"Recorded_Intensity.png");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  Back propagated wavefront  %%%%%%%%%%%%%%%%%%%%%%%%
[U_back,I_back,Phi_back] = Numerical_Propagation(U_n,Obj_Cam_Dist,840);

figure(4)
mesh(x,y,I_back);
title("Back Propagated Intensity pattern");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
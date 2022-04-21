I = double(rgb2gray(imread("Squiggly.jpeg")));
[a,b]=size(I);
x=0:a-1;
y=0:b-1;
[X,Y] = meshgrid(x,y);

I_fft = (fft2(I));
I_ifft = ifft2(I_fft);
mesh(X,Y,abs(I_fft));
mesh(X,Y,I_ifft);
Fx = 0:a-1;
Fy = 0:b-1;
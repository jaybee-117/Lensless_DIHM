%% Generate Optical Field for TIE.

% Amplitude: Uniform(Flag_Intensity==1) or Gaussian(Flag_Intensity==0)
% phi: Hann2 (Flag_Hann==1) or Arbitrary(Flag_Hann==0) 
% x2: light fieldis two times larger than CCD.

%% Set field of view smaller than light field in numerical propagration.

% Set the size of the Light Field.
Nx_LF = Nx*2;
Ny_LF = Ny*2;

% Cartesian coordinates.
[x, y] = meshgrid ...
    ( (-Nx_LF/2:Nx_LF/2-1)*Pixelsize ...
    , (-Ny_LF/2:Ny_LF/2-1)*Pixelsize);

% Set the Aperture.........................................................
if strcmpi(Flag_AperShape,'rect')
    Bd = 16;    % Boundary
    Lt = Nx_LF/4+1+Bd;  % Left
    Rt = Nx_LF/4+Nx-Bd; % Right
    Tp = Ny_LF/4+1+Bd;  % Top
    Bm = Ny_LF/4+Ny-Bd; % Bottom
    Aperture = false(Ny_LF,Nx_LF);
    Aperture(Tp:Bm,Lt:Rt) = true;

elseif strcmpi(Flag_AperShape,'circle')
    [Theta,R] = cart2pol(x,y);
    Rmax = max(R(~isnan(R)));
    Rmax_CCD = Rmax/2; % Note: Light field is larger than CCD, so Rmax/2.
    Aperture = R<0.65*Rmax_CCD; 

elseif strcmpi(Flag_AperShape,'circle-')
    [Theta,R] = cart2pol(x,y);
    Rmax = max(R(~isnan(R)));
    Rmax_CCD = Rmax/2; % Note: Light field is larger than CCD, so Rmax/2.
    Aperture = R<0.65*Rmax_CCD; 
    Aperture(y>0.6*Rmax_CCD) = 0;

elseif strcmpi(Flag_AperShape,'ellipse')
    Aperture = x.^2+y.^2-0.3*x.*y<max(x(~isnan(x)))*8e-5;

elseif strcmpi(Flag_AperShape,'ellipse-')
    Aperture = x.^2+y.^2-0.3*x.*y<max(x(~isnan(x)))*8e-5;
    [Theta,R] = cart2pol(x,y);
    Rmax = max(R(~isnan(R)));
    Rmax_CCD = Rmax/2; % Note: Light field is larger than CCD, so Rmax/2.
    Aperture(y>0.4*Rmax_CCD) = 0;

elseif strcmpi(Flag_AperShape,'telescope-')
    [Theta,R] = cart2pol(x,y);
    Rmax = max(R(~isnan(R)));
    Aperture = x.^2+y.^2+0.5*x.*y<0.8*Rmax*8e-5;
    Rmax_CCD = Rmax/2; % Note: Light field is larger than CCD, so Rmax/2.
    Aperture(y>0.4*Rmax_CCD) = 0;
    Aperture(R<0.2*Rmax_CCD) = 0;

elseif strcmpi(Flag_AperShape,'telescope')
    [Theta,R] = cart2pol(x,y);
    Rmax = max(R(~isnan(R)));
    Aperture = x.^2+y.^2-0.3*x.*y<Rmax*8e-5;
    Rmax_CCD = Rmax/2; % Note: Light field is larger than CCD, so Rmax/2.
%     Aperture(y>0.4*Rmax_CCD) = 0;
    Aperture(R<0.2*Rmax_CCD) = 0;
    Aperture(x>-4e-6 & x<4e-6) = 0;
    Aperture(y>-4e-6 & y<4e-6) = 0;
    Aperture(sqrt((x).^2+(y-1.6e-4).^2)<2e-5) = 0;
    Aperture(sqrt((x-1.3e-4).^2+(y+1e-4).^2)<2e-5) = 0;
    Aperture(sqrt((x+1.5e-4).^2+(y+1.2e-4).^2)<2e-5) = 0;

elseif strcmpi(Flag_AperShape,'full')
    Aperture = true(Ny_LF,Nx_LF);

elseif strcmpi(Flag_AperShape,'import')
    mask = double(rgb2gray(imread('mask.png')));
    mask = mask==max(mask(:));
    [Ny_mask, Nx_mask] = size(mask);
    Aperture = false(Ny_LF,Nx_LF);
    Aperture ...
        ( (Ny_LF-Ny_mask)/2+1 : (Ny_LF-Ny_mask)/2+Ny_mask ...
        , (Nx_LF-Nx_mask)/2+1 : (Nx_LF-Nx_mask)/2+Nx_mask ...
        ) = mask;
    
else
    
    error('Flag_AperShape mismatches.')  % Send error info.
end

% Assign phi.............................................................
if Flag_Hann == 0
    Phi = eval(vectorize(char(Phase_sym)));
    if numel(Phi)==1
        Phi = Phi*ones(Ny_LF, Nx_LF);
    end
elseif Flag_Hann == 1
%     phi = zeros(Ny_LF,Nx_LF);
    Phi = eval(vectorize(char(Phase_sym)));
    if numel(Phi)==1
        Phi = Phi*ones(Ny_LF, Nx_LF);
    end

    Hann2 = (hann(Ny/2)*hann(Nx/2)').^2;
    % 1:
    Phi(Ny_LF/4+Ny/4+1 : Ny_LF/4+Ny/4+Ny/2 ...
        , Nx_LF/4+Nx/4+1 : Nx_LF/4+Nx/4+Nx/2) ...
        =    Phi...
        ( Ny_LF/4+Ny/4+1 : Ny_LF/4+Ny/4+Ny/2 ...
        , Nx_LF/4+Nx/4+1 : Nx_LF/4+Nx/4+Nx/2) ...
        + Hann2;
    % 2:
    Shiftx = -20;
    Shifty = -50;
    Phi(Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        = Phi ...
        ( Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        + Hann2;
    
	% 3:
    Shiftx = -30;
    Shifty = 60;
    Phi(Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        = Phi ...
        ( Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        + Hann2;

elseif Flag_Hann == 2
%     phi = zeros(Ny_LF,Nx_LF);
    Phi = eval(vectorize(char(Phase_sym)));
    if numel(Phi)==1
        Phi = Phi*ones(Ny_LF, Nx_LF);
    end

    Hann2 = (hann(Ny/2)*hann(Nx/2)').^2;
    % 1:
    Phi(Ny_LF/4+Ny/4+1 : Ny_LF/4+Ny/4+Ny/2 ...
        , Nx_LF/4+Nx/4+1 : Nx_LF/4+Nx/4+Nx/2) ...
        =    Phi...
        ( Ny_LF/4+Ny/4+1 : Ny_LF/4+Ny/4+Ny/2 ...
        , Nx_LF/4+Nx/4+1 : Nx_LF/4+Nx/4+Nx/2) ...
        + Hann2;
    % 2:
    Shiftx = -20;
    Shifty = -50;
    Phi(Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        = Phi ...
        ( Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        + Hann2;

	% 3:
    Shiftx = -30;
    Shifty = 60;
    Phi(Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        = Phi ...
        ( Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        + Hann2;


	% 4:
    Shiftx = -110;
    Shifty = 80;
    Phi(Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        = Phi ...
        ( Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        + Hann2;


    
	% 5:
    Shiftx = 10;
    Shifty = 120;
    Phi(Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        = Phi ...
        ( Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        + Hann2;
    
	% 6:
    Shiftx = 110;
    Shifty = 20;
    Phi(Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        = Phi ...
        ( Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        + Hann2;

    % 7:
    Shiftx = 70;
    Shifty = -80;
    Phi(Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        = Phi ...
        ( Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        + Hann2;


    
    % 8:
    Shiftx = -90;
    Shifty = -135;
    Phi(Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        = Phi ...
        ( Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        + Hann2;


	% 9:
    Shiftx = -85;
    Shifty = -5;
    Phi(Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        = Phi ...
        ( Ny_LF/4+Ny/4+1+Shifty : Ny_LF/4+Ny/4+Ny/2+Shifty ...
        , Nx_LF/4+Nx/4+1+Shiftx : Nx_LF/4+Nx/4+Nx/2+Shiftx) ...
        + Hann2;

end

% Shift the phase bias to 0, which helps to check the phase error.
% phi = phi-mean(phi(~isnan(phi)));
% phi(~Aperture)=0;

% Assign Intensity according to the Flag_Intensity.........................
% Initialize the intensity matrix.
I0 = zeros(Ny_LF,Nx_LF);
if Flag_Intensity == 1
    I0(Aperture) = 1;                          % Uniform Intensity
elseif Flag_Intensity == 0
    Inten = eval(vectorize(char(Inten_sym)));
    I0(Aperture) = Inten(Aperture);            % Nonuniform Intensity
else
    error('Flag_Intensity should only be 0 or 1.')  % Send error info.
end



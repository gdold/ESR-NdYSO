% --- Spin system ---
I_nuc = 3.5;

Sys = struct();
Sys.Nucs = '145Nd';
%Sys.I = 3.5; %Specified by Nucs(?)
Sys.S = 0.5;
%Sys.g = [-1.47 -1.02 4.20];


%% g and A tensors as principal values %%

%Sys.g = [-0.96 1.48 -4.14];
%Sys.A = [-0.146 -495.5 -808.0]; %MHz

%gFrame_zyz = [12.4 38.6 -86.4]*(pi/180);
%AFrame_zyz = [106.4 37.5 -90.7]*(pi/180);

%Sys.gFrame = eulang(inv(erot(gFrame_zyz))); % try inverting matrix
%Sys.AFrame = eulang(inv(erot(AFrame_zyz)));

%Sys.gFrame = gFrame_zyz;
%Sys.AFrame = AFrame_zyz;

%% g and A tensors as full 3x3 tensor %%

Sys.g = [1.30 0.62 0.22;
        0.62 -2.07 1.62;
        0.22 1.62 -2.86];
Sys.A = [-37.1 -99.9 -83.4;
        -99.9 -589.2 169.4;
        -83.4 169.4 -678.4];

% The eigenvaues of these matrices
% do not match up with the above principal values
% Is this right?!

% --- Experiment properties ---

% Mag field close to D1 axis
Exp = struct();
Exp.mwFreq = 9.385; %GHz
Exp.Range = [0 600]; %mT % [350 600]
Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 

% Calculate euler angles given xyz rotations
phi = 45*(pi/180);
theta = 45*(pi/180);
rho = 90*(pi/180);
cryst_rot = eulang(rotaxi2mat(ang2vec(phi,theta),rho));

Exp.CrystalOrientation = cryst_rot; %Euler angles
Exp.Temperature = 20; %Kelvin

% --- Pepper ---
pepper(Sys,Exp)


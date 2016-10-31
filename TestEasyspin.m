%% Spin system %%
I_nuc = 3.5;

Sys = struct();
Sys.Nucs = '145Nd';
%Sys.I = 3.5; %Specified by Nucs(?)
Sys.S = 0.5;


%% g and A tensors from source %%
% Allowed sources:
% 'Maier-FlaigTensor'
% 'Maier-FlaigPrincipal'
% 'Wolfowicz' - corrected Euler convention
parameter_source = 'Wolfowicz';
Sys = NdYSOparams(Sys,parameter_source); % Appends chosen parameters to Sys

% The eigenvaues M-F Tensor do not match up with
% M-F Principal - is this right?!

%% Experiment properties %%

% Mag field close to D1 axis
Exp = struct();
Exp.mwFreq = 9.385; %GHz
Exp.Range = [0 600]; %mT % [350 600]
Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 


%% Crystal rotation %%
% Calculate euler angles given xyz rotations
phi = 45*(pi/180);
theta = 45*(pi/180);
rho = 90*(pi/180);
cryst_rot = eulang(rotaxi2mat(ang2vec(phi,theta),rho));

Exp.CrystalOrientation = cryst_rot; %Euler angles
Exp.Temperature = 20; %Kelvin

% --- Pepper ---
pepper(Sys,Exp)


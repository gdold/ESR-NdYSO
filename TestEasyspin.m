%% Spin system %%
I_nuc = 3.5;

Sys = struct();
Sys.Nucs = '145Nd';
%Sys.I = 3.5; %Specified by Nucs(?)
Sys.S = 0.5;


%% g and A tensors from source %%
% Allowed sources:
% 'Maier-FlaigTensor'
% 'Maier-FlaigPrincipal' % gives vastly different results to other two
% 'Wolfowicz' - corrected Euler convention
parameter_source = 'Wolfowicz';
Sys = NdYSOparams(Sys,parameter_source); % Appends chosen parameters to Sys

% The eigenvaues M-F Tensor do not match up with
% M-F Principal - is this right?!

%% Experiment properties %%

% Mag field close to D1 axis
Exp = struct();
Exp.mwFreq = 9.7; %GHz
Exp.Range = [0 600]; %mT % [350 600]
Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 


%% Crystal rotation %%
% Calculate euler angles given xyz rotations
phi = 69.83;
theta = 3.75;
rho = 80;
cryst_rot = eulang( rotz(phi) * roty(theta) * rotx(rho) );



Exp.CrystalOrientation = cryst_rot; %Euler angles
Exp.Temperature = 20; %Kelvin

% --- Pepper ---
pepper(Sys,Exp)


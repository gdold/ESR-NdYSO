deg = pi/180;

%% Spin system %%
Sys = struct();
Sys.Nucs = '145Nd';
%Sys.I = 3.5; %Specified by Nucs(?)
Sys.S = 0.5;


%% g and A tensors from source %%
% Allowed sources:
% 'Maier-FlaigTensor'
% 'Maier-FlaigPrincipal' % gives vastly different results to other two
% 'Wolfowicz' - corrected Euler convention
%parameter_source = 'Maier-FlaigPrincipal';
parameter_source = 'Maier-FlaigTensor';
Sys = NdYSOparams(Sys,parameter_source); % Appends chosen parameters to Sys

% The eigenvaues M-F Tensor do not match up with
% M-F Principal - is this right?!

%% Experiment properties %%

Exp = struct();
Exp.mwFreq = 9.385; %GHz % 9.7 perhaps?
Exp.Range = [0 1000]; %mT % [350 600]
Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 
Exp.Temperature = 20; %Kelvin

%% Crystal rotation %%
% Axes defined in crystal basis [D1 D2 b]

% CHOOSE WHETHER TO REPRODUCE FIG. 4.4a, 4.4b, 4.4c
% 4.4a
% rot_axis = ang2vec(69.83*deg,3.75*deg); % Fig 4.4a
% init_mag_vect = [0.999745;-0.001384;-0.022508];
% 4.4b
rot_axis = ang2vec(189.13*deg,96.21*deg); % Fig 4.4b
init_mag_vect = [-0.156794;0.987480;-0.017280];
% 4.4c
% rot_axis = ang2vec(89.72*deg,-92.77*deg); % Fig 4.4c
% init_mag_vect = [0.99998808;-0.004875;-2.358e-04];

% The initial rotation from the crystal frame to the lab frame is then
init_mag_vect = init_mag_vect/norm(init_mag_vect);
y_axis_in_cryst = cross(init_mag_vect,rot_axis); % true if rot,mag orthogonal
init_rotm = [rot_axis,y_axis_in_cryst,init_mag_vect]'; % orthonormal?


% Want a threshold intensity below which eigfields will ignore transition?
Opt = struct();
%Opt.Threshold = 100;

% Number of steps in rotation simulation
rot_steps = 72;
total_angle = 360*deg;

%% Calculate spectrum %%

Rot_inc = rotaxi2mat(rot_axis,total_angle/rot_steps);
Rot_inc_lab = rotaxi2mat([1,0,0],-total_angle/rot_steps);
Rot_inc_lab = eye(3);

cryst_rot = init_rotm;
cryst_rot_lab = init_rotm;
x = [];
y1 = [];
y2 = [];
for n = 0:rot_steps
    disp(n);
    angle = n*total_angle/rot_steps;
    %mag_vect_crystal = cryst_rot*init_mag_vect;
    %Exp.CrystalOrientation = eulang(alignMagRot(mag_vect_crystal));
    Exp.CrystalOrientation = eulang(cryst_rot);
    out = eigfields(Sys,Exp,Opt);
    fields1 = out{1}';
    fields2 = out{2}';
    x_temp = linspace(angle,angle,length(fields1))/deg;
    
    x = [x x_temp];
    y1 = [y1 fields1];
    y2 = [y2 fields2];
    
    cryst_rot = Rot_inc_lab*cryst_rot*Rot_inc'; %Inverse crystal rotation
end

figure
hold off
scatter(x,y1,'.')
hold on
scatter(x,y2,'.')
xlabel('Angle (degrees)')
ylabel('B (mT)')
text_label = {['Source: ',parameter_source],['Rotation axis: ',num2str(rot_axis')]};
annotation('textbox',[.2 .5 .3 .3],'string',text_label,'FitBoxToText','on');
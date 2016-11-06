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
%rot_axis = ang2vec(69.83*deg,3.75*deg); % Fig 4.4a
%axis_wrt = [1,0,0]';
% 4.4b
rot_axis = ang2vec(189.13*deg,96.21*deg); % Fig 4.4b
axis_wrt = [0,1,0]';
% 4.4c
%rot_axis = ang2vec(89.72*deg,-92.77*deg);% Fig 4.4c
%axis_wrt = [1,0,0]';

% Can define initial B field and rotation axis in crystal frame
%init_mag_vect = [1,0,0]; % a,c: [1,0,0], b:[0,1,0]
%rot_axis = [0 0 1]; % D1 D2 b

% Relevant rotation axes for reproducing Maier-Flaig's figures
%init_mag_approx = [1,0,0]; % a,c: [1,0,0], b:[0,1,0]

% To reproduce Maier-Flaig's figures, need B-field orthogonal to rot_axis
% as rot_axis along zL and B along xL
orth_axis = cross([0,0,1],rot_axis);
init_mag_vect = orth_axis/norm(orth_axis);

% Rotate init_mag_vect to be close to reference axis
proj = axis_wrt - dot(axis_wrt,rot_axis)*rot_axis;% projection onto normal plane
theta_offset = atan2(norm(cross(proj,init_mag_vect)),dot(proj,init_mag_vect));
init_mag_vect = rotaxi2mat(rot_axis,theta_offset)*init_mag_vect';

% Want a threshold intensity below which eigfields will ignore transition?
Opt = struct();
%Opt.Threshold = 100;

% Number of steps in rotation simulation
rot_steps = 72;
total_angle = 360*deg;

%% Calculate spectrum %%

Rot_inc = rotaxi2mat(rot_axis,total_angle/rot_steps);

cryst_rot = eye(3);
x = [];
y1 = [];
y2 = [];
for n = 0:rot_steps
    disp(n);
    angle = n*total_angle/rot_steps;
    mag_vect_crystal = cryst_rot*init_mag_vect;
    Exp.CrystalOrientation = eulang(alignMagRot(mag_vect_crystal));
    out = eigfields(Sys,Exp,Opt);
    fields1 = out{1}';
    fields2 = out{2}';
    x_temp = linspace(angle,angle,length(fields1))/deg;
    
    x = [x x_temp];
    y1 = [y1 fields1];
    y2 = [y2 fields2];
    
    cryst_rot = Rot_inc*cryst_rot;
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
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
Exp.mwFreq = 9.7; %GHz
Exp.Range = [0 1000]; %mT % [350 600]
Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 
Exp.Temperature = 20; %Kelvin

%% Crystal rotation %%
% [D1 D2 b]
init_mag_vect = [1,0,0]; % a,c: [1,0,0], b:[0,1,0]
%rot_axis = [0 0 1]; % D1 D2 b
rot_axis = ang2vec(69.83*deg,3.75*deg); % Fig 4.4a
%rot_axis = ang2vec(189.13*deg,96.21*deg); % Fig 4.4b % Not quite right
%rot_axis = ang2vec(89.72*deg,-92.77*deg);% Fig 4.4c

%Find some axis orthogonal to rot_axis
orth_axis = cross([0,0,1],rot_axis);
init_mag_vect = orth_axis/norm(orth_axis);

init_rot = alignMagRot(init_mag_vect); % matrix


Exp.CrystalOrientation = eulang(init_rot); %Euler angles


Opt = struct();
%Opt.Threshold = 100;



% Calculate rotation increment for each step
% given rotation axis in crystal frame,
% total angle, steps
rot_steps = 72;


total_angle = 360*deg;

Rot_inc = rotaxi2mat(rot_axis,total_angle/rot_steps);

cryst_rot = eye(3);
x = [];
y1 = [];
y2 = [];
for n = 0:rot_steps
    disp(n);
    angle = n*total_angle/rot_steps;
    mag_vect_crystal = cryst_rot*init_mag_vect';
    Exp.CrystalOrientation = eulang(alignMagRot(mag_vect_crystal));
    out = eigfields(Sys,Exp,Opt);
    fields1 = out{1}';
    fields2 = out{2}';
    x_temp = linspace(angle,angle,length(fields1));
    
    x = [x x_temp];
    y1 = [y1 fields1];
    y2 = [y2 fields2];
    
    cryst_rot = Rot_inc*cryst_rot;
    % perform rotation in crystal frame before converting to lab frame
end

figure
hold off
scatter(x,y1,'.')
hold on
scatter(x,y2,'.')
xlabel('Angle (radians)')
ylabel('B (mT)')
text_label = {['Source: ',parameter_source],['Rotation axis: ',num2str(rot_axis)]};
annotation('textbox',[.2 .5 .3 .3],'string',text_label,'FitBoxToText','on');
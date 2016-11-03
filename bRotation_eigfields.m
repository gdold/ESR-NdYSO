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

% Mag field close to D1 axis
Exp = struct();
Exp.mwFreq = 9.7; %GHz
Exp.Range = [0 1000]; %mT % [350 600]
Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 
Exp.Temperature = 20; %Kelvin

%% Crystal rotation %%
%cryst_to_lab = rotz(135*deg)*roty(-90*deg)*rotz(180*deg);
%cryst_to_lab = eye(3); %B along b
cryst_to_lab = roty(pi/2); %B along D1
%cryst_to_lab = rotx(-pi/2); % B along D2
% Could issues be coming from imperfect B-field direction?

cryst_rot_mat = cryst_to_lab;
Exp.CrystalOrientation = eulang(cryst_rot_mat); %Euler angles


Opt = struct();
%Opt.Threshold = 100;



% Calculate rotation increment for each step
% given rotation axis in crystal frame,
% total angle, steps
steps = 12;
%axis = [0 0 1]; % x y z
axis = ang2vec(69.83*deg,3.75*deg)'; % Fig 4.4a
%axis = ang2vec(189.13*deg,96.21*deg)'; % Fig 4.4b
%axis = ang2vec(89.72*deg,-92.77*deg)';% Fig 4.4c

total_angle = 360*deg;
Rot_inc = rotaxi2mat(axis,total_angle/steps);

x = [];
y = [];
for n = 0:steps
    disp(n);
    angle = n*total_angle/steps;
    Exp.CrystalOrientation = eulang(cryst_rot_mat);
    out = eigfields(Sys,Exp,Opt);
    fields = out{1}';
    x_temp = linspace(angle,angle,length(fields));
    
    x = [x x_temp];
    y = [y fields];
    
    cryst_rot_mat = cryst_rot_mat*Rot_inc; % perform rotation in crystal frame before converting to lab frame
end

scatter(x,y,'.')
xlabel('Angle (radians)')
ylabel('B (mT)')
text_label = {['Source: ',parameter_source],['Rotation axis: ',num2str(axis)]};
annotation('textbox',[.2 .5 .3 .3],'string',text_label,'FitBoxToText','on');
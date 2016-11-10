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
%Exp.mwFreq = 9.385; %GHz % 9.7 perhaps?
Exp.Range = [0 3]; %mT % [350 600]
%Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 
Exp.Temperature = 20; %Kelvin

%% Crystal rotation %%
% Axes defined in crystal basis [D1 D2 b]

% DIRECT MAGNETIC FIELD ALONG D1
init_rotm = [0,0,1; %  b->x
             0,-1,0;% D2->-y
             1,0,0];% D1->z
rt = 1.0/sqrt(2);
init_rotm = rotaxi2mat([rt,0,rt],180*deg);

init_rotm = alignMagRot([0,1,0]);

% The initial rotation from the crystal frame to the lab frame is then


% Want a threshold intensity below which eigfields will ignore transition?
Opt = struct();
%Opt.Threshold = 0.1;

% Number of steps in rotation simulation
rot_steps = 144;
total_angle = 360*deg;

%% Calculate spectrum %%

%Rot_inc = rotaxi2mat(rot_axis,total_angle/rot_steps);
Rot_inc_lab = rotaxi2mat([1,0,0],-total_angle/rot_steps);
Rot_inc = eye(3);

cryst_rot = init_rotm;
cryst_rot_lab = init_rotm;
x = [];
y1 = [];
y2 = [];
y = [];
freqs = [];

Exp.Field = 20; % mT % WANT TO SWEEP MW FIELD
Exp.mwFreq = 2.0; % GHz
threshold = 0.1*0.01;
Opt.Threshold = 0;

max_field = 50;
field_steps = 500;
for n = 0:field_steps
    freqs = [];
    disp(n);
    field = n*max_field/field_steps;
    Exp.Field = field;
    Exp.CrystalOrientation = eulang(cryst_rot);
    [Pos,Amp,Wid,Trans] = resfreqs_matrix(Sys,Exp,Opt);
    out = [Pos,Amp]';
    for i = 1:length(out)
        if out(2,i) >= threshold
            freqs = [freqs out(1,i)];
        end
    end
    %fields1 = out{1}';
    %fields2 = out{2}';
    %freqs = out(1,:);
    x_temp = linspace(field,field,length(freqs));
    
    x = [x x_temp];
    %y1 = [y1 fields1];
    %y2 = [y2 fields2];
    y = [y freqs];
    
end

% for n = 0:rot_steps
%     freqs = [];
%     disp(n);
%     angle = n*total_angle/rot_steps;
%     %mag_vect_crystal = cryst_rot*init_mag_vect;
%     %Exp.CrystalOrientation = eulang(alignMagRot(mag_vect_crystal));
%     Exp.CrystalOrientation = eulang(cryst_rot);
%     [Pos,Amp,Wid,Trans] = resfreqs_matrix(Sys,Exp,Opt);
%     out = [Pos,Amp]';
%     for i = 1:length(out)
%         if out(2,i) >= threshold
%             freqs = [freqs out(1,i)];
%         end
%     end
%     %fields1 = out{1}';
%     %fields2 = out{2}';
%     %freqs = out(1,:);
%     x_temp = linspace(angle,angle,length(freqs))/deg;
%     
%     x = [x x_temp];
%     %y1 = [y1 fields1];
%     %y2 = [y2 fields2];
%     y = [y freqs];
%     
%     cryst_rot = Rot_inc_lab*cryst_rot*Rot_inc'; %Inverse crystal rotation
% end

figure
hold off
%scatter(x,y1,'.')
scatter(x,y,'.')
hold on
%scatter(x,y2,'.')
xlabel('Angle (degrees)')
ylabel('B (mT)')
text_label = {['Source: ',parameter_source],['Rotation axis: ',num2str(rot_axis')]};
annotation('textbox',[.2 .5 .3 .3],'string',text_label,'FitBoxToText','on');
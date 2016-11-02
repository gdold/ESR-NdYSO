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
parameter_source = 'Wolfowicz';
Sys = NdYSOparams(Sys,parameter_source); % Appends chosen parameters to Sys

% The eigenvaues M-F Tensor do not match up with
% M-F Principal - is this right?!

%% Experiment properties %%

% Mag field close to D1 axis
Exp = struct();
Exp.mwFreq = 9.7; %GHz
Exp.Range = [0 1000]; %mT % [350 600]
Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 


%% Crystal rotation %%
% Calculate euler angles given xyz rotations
%phi = 69.83;
%theta = 3.75;
%rho = 0;
phi = 69.83;
theta = 86;
rho = 86;
cryst_ang = [phi, theta, rho];
cryst_rot = eulang( rotz(cryst_ang(1)) * roty(cryst_ang(2)) * rotx(cryst_ang(3)) );



Exp.CrystalOrientation = cryst_rot; %Euler angles
Exp.Temperature = 20; %Kelvin

Opt = struct();
%Opt.Threshold = 100;

% --- Eigenfields ---
%out = eigfields(Sys,Exp,Opt);

num_of_points = 3;
angles = linspace(0,360,num_of_points)';
x = [];
y = [];
angle_choice = 1; % phi, theta, rho

for n = 1:length(angles)
    disp(n);
    cryst_ang(angle_choice) = angles(n);
    Exp.CrystalOrientation = eulang( rotz(cryst_ang(1)) * roty(cryst_ang(2)) * rotx(cryst_ang(3)) );
    out = eigfields(Sys,Exp,Opt);
    fields = out{1}';
    %fields = fields((length(fields)-16):length(fields)); % exclude all but 16 strongest transitions
    %x_temp = linspace(rho,rho,length(fields));
    x_temp = linspace(cryst_ang(angle_choice),cryst_ang(angle_choice),length(fields));
    
    x = [x x_temp];
    y = [y fields];
end

scatter(x,y,'.')
xaxis_texts = {'Phi','Theta','Rho'};
xlabel(xaxis_texts(angle_choice))
ylabel('B (mT)')
text_label = {['Source: ',parameter_source],
              ['Phi: ',num2str(phi),' Theta: ',num2str(theta),' Rho: ',num2str(rho)]};
annotation('textbox',[.2 .5 .3 .3],'string',text_label,'FitBoxToText','on');
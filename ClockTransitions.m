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
Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 
Exp.Temperature = 20; %Kelvin

%% Crystal rotation %%
% Axes defined in crystal basis [D1 D2 b]

% DIRECT MAGNETIC FIELD ALONG D1
init_rotm = [0,0,1; %  b->x
             0,-1,0;% D2->-y
             1,0,0];% D1->z
rt = 1.0/sqrt(2);
init_rotm = rotaxi2mat([rt,0,rt],180*deg);

init_axis = [1,0,0];
init_rotm = alignMagRot(init_axis);

rot_axis = [0,0,1];

% The initial rotation from the crystal frame to the lab frame is then


% Want a threshold intensity below which eigfields will ignore transition?
Opt = struct();
Opt.Threshold = 0; % Get all transitions (incl forbidden)
Opt.Sites = [1,2];

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
y = [];
freqs = [];

Exp.Field = 20; % mT % WANT TO SWEEP MW FIELD
Exp.mwFreq = 2.0; % GHz
threshold = 0.1*0.01;
Opt.Threshold = 0;

max_field = 50;
field_steps = 500;


for step = 0:rot_steps
    disp(['Step ',int2str(step),' of ',int2str(rot_steps)]);
    angle = step*total_angle/rot_steps;

    
    output = [];
    
    for n = 0:field_steps
        freqs = [];
        %disp(n);
        mag_field = n*max_field/field_steps;
        Exp.Field = mag_field;
        Exp.CrystalOrientation = eulang(cryst_rot);
        [Pos,Amp,Wid,Trans] = resfreqs_matrix(Sys,Exp,Opt);
        Site = repelem(Opt.Sites(1),length(Trans))';
        if length(Trans)<length(Pos)
            Trans = [Trans;Trans]; % Add extra labels for simulation of site 2
            Site = [repelem(Opt.Sites(1),length(Site))';repelem(Opt.Sites(2),length(Site))'];
        end
        transition_label = Trans(:,1)*100 + Trans(:,2);% label x->y by int xxyy, works only if <100 transitions
        field = repelem(mag_field,length(Pos))';
        out = [field,Pos,Amp,transition_label,Site]; % label is cast to double
        output = [output; out];
        for i = 1:1:length(out)
            %if out(i,3) >= threshold
            %if (out(i,4) == 15) && (out(i,5) == 16)
            if 1
                freqs = [freqs out(i,2)];
            end
        end
        %fields1 = out{1}';
        %fields2 = out{2}';
        %freqs = out(1,:);
        %x_temp = linspace(mag_field,mag_field,length(freqs));
        
        %x = [x x_temp];
        %y1 = [y1 fields1];
        %y2 = [y2 fields2];
        %y = [y freqs];
        
    end
    
    % Find transitions with large amplitude
    largeamp = [];
    threshold = 0.1;
    output(:,3) = output(:,3)/max(output(:,3));
    for i = 1:length(output(:,3))
        if output(i,3) > threshold
            largeamp = [largeamp; output(i,4)];
        end
    end
    largeamp = unique(largeamp);
    
    % Create array containing only these transitions
    transitions = [];
    for i = 1:length(output(:,4))
        if any(abs( output(i,4) - largeamp )<1e-10) % transition has a large amplitude at some point
            transitions = [transitions; output(i,:)];
        end
    end
    
    x = transitions(:,1); % field
    y = transitions(:,2); % freq
    z = transitions(:,3); % intensity
    
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
    
    %figure
    init_axis_str = 'D1';
    rot_axis_str = 'b';
    
    %hold off
    %scatter(x,y1,'.')
    scatter(x,y,[],z,'.')
    colormap(flipud(hot))
    cbar = colorbar();
    cbar.Label.String = 'Relative amplitude';
    caxis([0.0,1.0])
    ylim([0,3000])
    %hold on
    %scatter(x,y2,'.')
    xlabel('B (mT)')
    ylabel('Transition frequency (MHz)')
    title(['Init axis: ',init_axis_str,'; Rot axis: ',rot_axis_str,'; angle: ',num2str(angle/deg)])
    saveas(gcf,['figure',int2str(step),'.png'])
    %text_label = {['Source: ',parameter_source],['Rotation axis: ',num2str(rot_axis')]};
    %annotation('textbox',[.2 .5 .3 .3],'string',text_label,'FitBoxToText','on');
    
    cryst_rot = Rot_inc_lab*cryst_rot*Rot_inc'; %Inverse crystal rotation
end



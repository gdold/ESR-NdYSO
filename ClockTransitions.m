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

%% Experimental parameters %%

Exp = struct();
%Exp.mwFreq = 9.385; %GHz % 9.7 perhaps?
Exp.Range = [0 3]; %mT % [350 600]
Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 
Exp.Temperature = 20; %Kelvin
Exp.Mode = 'perpendicular'; % direction of MW field, 'perpendicular' or 'parallel'

%% Crystal rotation %%
% Axes defined in crystal basis [D1 D2 b]

magaxis = 'b';
MWaxis = 'b';
[init_rotm, Exp.Mode] = setInitialAxes(magaxis,MWaxis); % currently only crystal axes

rot_axis = [0,1,0];

% The initial rotation from the crystal frame to the lab frame is then


% Want a threshold intensity below which eigfields will ignore transition?
Opt = struct();
Opt.Threshold = 0; % Get all transitions (incl forbidden)
%Opt.Sites = [1,2];
Opt.Sites = [1];



%% Calculate spectrum %%



cryst_rot = init_rotm;
cryst_rot_lab = init_rotm;
x = [];
y = [];
freqs = [];

Exp.Field = 20; % mT % WANT TO SWEEP MW FIELD
Exp.mwFreq = 2.0; % GHz
threshold = 0.1*0.01;
Opt.Threshold = 0;

% Set up field sweep
max_field = 50;
field_steps = 500;

% Set up rotation
total_angle = 360*deg;
rot_steps = 0; % 0 for no angular sweep

Rot_inc = rotaxi2mat(rot_axis,total_angle/rot_steps);
%Rot_inc_lab = rotaxi2mat([1,0,0],-total_angle/rot_steps);
Rot_inc_lab = eye(3);


full_empty = struct();
full_empty.data = struct();
full_empty.angle = [];


dat_empty = struct();
NaNarray = repelem(NaN,field_steps+1,1); % preallocate memory for vects
dat_empty.transition = [];
dat_empty.label = [];
dat_empty.site = [];
dat_empty.frequency = NaNarray;
dat_empty.field = NaNarray;
dat_empty.amplitude = NaNarray;


full = full_empty();

for step = 0:rot_steps
    disp(['Step ',int2str(step),' of ',int2str(rot_steps)]);
    angle = step*total_angle/rot_steps;
    
    full(step+1) = full_empty;
    full(step+1).angle = angle;
    
    
    Pos = [];
    Amp = [];
    Trans = [];
    Site = [];
    
    dat = dat_empty();
    %dat.angle = angle;
    
    
    
    for n = 0:field_steps
        freqs = [];
        %disp(n);
        mag_field = n*max_field/field_steps;
        Exp.Field = mag_field;
        Exp.CrystalOrientation = eulang(cryst_rot);
        [Pos,Amp,~,Trans] = resfreqs_matrix(Sys,Exp,Opt);
        Site = repelem(Opt.Sites(1),length(Trans))';
        if length(Trans)<length(Pos)
            Trans = [Trans;Trans]; % Add extra labels for simulation of site 2
            Site = [repelem(Opt.Sites(1),length(Site))';repelem(Opt.Sites(2),length(Site))'];
        end
        transition_label = Site(:)*10000 + Trans(:,1)*100 + Trans(:,2);% label x->y by int xxyy, works only if <100 transitions
        field = repelem(mag_field,length(Pos))';
        %out = [field,Pos,Amp,transition_label,Site];% label is cast to double
        
        dat(1).transition = [1 2]; % prevent duplicate...
        dat(1).label = 10102; % lazy hack...
        
        for i = 1:length(Trans)
            [in,loc] = ismember(transition_label(i),[dat(:).label]);
            if ~in
                dat(end+1) = dat_empty;
                loc = length(dat);
            end
            dat(loc).transition = Trans(i,:);
            dat(loc).label = transition_label(i);
            dat(loc).site = Site(i);
            dat(loc).frequency(n+1) = Pos(i)';
            dat(loc).field(n+1) = mag_field';
            dat(loc).amplitude(n+1) = Amp(i)';
        end
        
    end
    
    threshold = 0.1;
    strong_transitions = findStrongTransitions(dat,threshold);

    x = [dat(strong_transitions).field];
    x = x(:); % make column vector
    y = [dat(strong_transitions).frequency];
    y = y(:);
    z = [dat(strong_transitions).amplitude];
    z = z(:);
    
    %figure
    rot_axis_str = 'D2';
    
    %hold off
    %scatter(x,y,'.')
    scatter(x,y,[],z,'.')
    colormap(flipud(hot))
    cbar = colorbar();
    cbar.Label.String = 'Amplitude';
    %caxis([0.0,1.0])
    ylim([0,3000])
    %hold on
    %scatter(x,y2,'.')
    xlabel('B (mT)')
    ylabel('Transition frequency (MHz)')
    %title(['Mag vector: ',num2str(init_mag_vect')])
    title(['Mag axis: ',magaxis,'; MW axis: ',MWaxis,'; Rot axis: ',rot_axis_str,'; angle: ',num2str(angle/deg)])
    %tit = 'TEST';
    %title(tit);
    findClockTransitions;
    %saveas(gcf,['figure',int2str(5*step),'.png'])
    %saveas(gcf,[tit,'.png'])
    %save(['output',int2str(5*step),'.mat'],'output')
    %text_label = {['Source: ',parameter_source],['Rotation axis: ',num2str(rot_axis')]};
    %annotation('textbox',[.2 .5 .3 .3],'string',text_label,'FitBoxToText','on');
    
    full(step+1).data = dat;
    
    cryst_rot = Rot_inc_lab*cryst_rot*Rot_inc'; %Inverse crystal rotation
end

save(['full_',datestr(now,'yyyy-mm-ddTHH-MM-SS'),'.mat'],'full')

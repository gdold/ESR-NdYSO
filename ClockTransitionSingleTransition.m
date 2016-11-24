deg = pi/180;
rt = 1.0/sqrt(2);

%% Spin system %%
Sys = struct();
% g and A tensors from source %
% Allowed sources:
% 'Maier-FlaigTensor'
% 'Maier-FlaigPrincipal' % gives vastly different results to other two
% 'Wolfowicz' - corrected Euler convention
parameter_source = 'Maier-FlaigTensor';
Sys = NdYSOparams(Sys,parameter_source); % Appends chosen parameters to Sys

%% Experimental parameters %%

Exp = struct();
Exp.Range = [0 3]; %GHz - must start at 0 to number transitions consistently
Exp.CrystalSymmetry = 'C2h'; %monoclinic C^6_2h spacegroup 
Exp.Temperature = 20; %Kelvin

%% Other simulation parameters %%
% Want a threshold intensity below which eigfields will ignore transition?
Opt = struct();
Opt.Threshold = 0; % Get all transitions (incl forbidden)
Opt.Sites = 1; % [1,2] for both YSO sites

%% Sweep parameters %%
% Axes defined in crystal basis [D1 D2 b]

% Set up field sweep parameters
max_field = 50;
field_steps = 500;
magaxis = 'D2';
MWaxis = 'b';
[init_rotm, Exp.Mode] = setInitialAxes(magaxis,MWaxis); % currently only crystal axes
%Exp.Mode = [90*deg,0*deg];

% Set up rotation parameters
rotaxis = 'D1'; % set this manually below
rot_axis = [1,2,0]; % rotate mag field in crystal frame
start_angle = 0*deg;
total_angle = 90*deg;
rot_steps = 2; % 0 for no angular sweep
rot_points = linspace(start_angle,start_angle+total_angle,rot_steps+1);

%rot_axis_lab = [1,0,0]; % lab frame

cryst_rot = init_rotm;
cryst_rot_lab = init_rotm;
Rot_inc = rotaxi2mat(rot_axis,total_angle/rot_steps);
%Rot_inc_lab = rotaxi2mat(rot_axis_lab,-total_angle/rot_steps);
Rot_inc_lab = eye(3);

%% Run sweep %%

% Full struct contains all dat structs from rotation sweep
full_empty = struct();
full_empty.data = struct();
full_empty.clocks = struct();
full_empty.angle = [];

% Dat struct contains data from each field sweep
dat_empty = struct();
NaNarray = repelem(NaN,field_steps+1,1); % preallocate memory for vects
dat_empty.transition = [];
dat_empty.label = [];
dat_empty.site = [];
dat_empty.frequency = NaNarray;
dat_empty.field = NaNarray;
dat_empty.amplitude = NaNarray;
dat_empty.angle = NaN;


full = full_empty();

for step = 0:rot_steps
    angle = rot_points(step+1);
    %disp(['Step ',int2str(step),' of ',int2str(rot_steps)]);
    fprintf('\nStep %i of %i\t',step,rot_steps)
    full(step+1) = full_empty;
    full(step+1).angle = angle;
    
    
    Pos = [];
    Amp = [];
    Trans = [];
    Site = [];
    
    dat = dat_empty();
    
    for n = 0:field_steps
        mag_field = n*max_field/field_steps;
        fprintf('%3.f%% - %2.0f mT',100*n/field_steps,mag_field) % display percentage complete and current field
        %disp([int2str(mag_field),'mT'])
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
            dat(loc).angle = angle;
        end
        
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8)%delete last 12 characters
    end
    fprintf('\n')
    
    threshold = 0.1;
    strong_transitions = findStrongTransitions(dat,threshold);

    x = [dat(strong_transitions).field];
    x = x(:); % make column vector
    y = [dat(strong_transitions).frequency];
    y = y(:);
    z = [dat(strong_transitions).amplitude];
    z = z(:);
    
    %figure
    

    scatter(x,y,[],z,'.')
    colormap(flipud(hot))
    cbar = colorbar();
    cbar.Label.String = 'Amplitude';
    %caxis([0.0,1.0])
    ylim(Exp.Range*1000)
    xlabel('B (mT)')
    ylabel('Transition frequency (MHz)')
    title(['Mag axis: ',magaxis,'; MW axis: ',MWaxis,'; Rot axis: ',rotaxis,'; angle: ',num2str(angle/deg)])
    clocks = findClockTransitions(dat,threshold);
    saveas(gcf,['figure',int2str(step),'.png'])
    %text_label = {['Source: ',parameter_source],['Rotation axis: ',num2str(rot_axis')]};
    %annotation('textbox',[.2 .5 .3 .3],'string',text_label,'FitBoxToText','on');
    
    full(step+1).data = dat;
    full(step+1).clocks = clocks;
    
    cryst_rot = Rot_inc_lab*cryst_rot*Rot_inc'; %Inverse crystal rotation
end

%% Look for clock transitions that meet criteria %%
promising_clocks = struct([]);
freqrange = [1000, 3000]; % MHz
max_deriv2 = 10; % MHz/mT^2...
min_amplitude_relative = 0.1; % normalised to 1
num_displayed_clocks = 20;



spec = [7,12]; % specify index interesting
num_of_levels = 0.5*(1+sqrt(1+8*length(full(1).data(:))));
transition_indices = sum(num_of_levels-(1:spec(1)))-num_of_levels+spec(2);
transition_str = [int2str(spec(1)),'-->',int2str(spec(2))];



full_peak_amplitude = 0;
for i = 1:length(full) % neater way to do this?
    full_peak_amplitude = max([[full(i).clocks(:).amplitude], full_peak_amplitude]);
end
min_amplitude_absolute = min_amplitude_relative*full_peak_amplitude;

x = [];
y = [];
z = [];

if ~isempty(transition_index)
    for i = 1:length(full)
        for transition = 1:length(transition_index)
            idx = find([full(i).clocks(:).index] == transition_index(transition));
            x = [x, repelem(full(i).angle/deg,length([full(i).clocks(idx).frequency]))];
            y = [y, [full(i).clocks(idx).field]];
            %y = [y, [full(i).clocks(idx).field]];
            z = [z, [full(i).clocks(idx).amplitude]];
        end
    end
end
scatter(x,y,[],z,'.')
colormap(flipud(hot))
cbar = colorbar();
cbar.Label.String = 'Amplitude';
%caxis([0.0,1.0])
%ylim(Exp.Range*1000)
xlabel('Angle (degrees)')
ylabel('Clock transition field (mT)')
title(['Mag axis: ',magaxis,'; MW axis: ',MWaxis,'; Rot axis: ',rotaxis,'; transition: ',transition_str])
saveas(gcf,['transitions',transition_str,'.png'])


% for i = 1:length(full)
%     if isempty(fieldnames(clocks))
%         continue
%     end
%     for j = 1:length(full(i).clocks)
%         if full(i).clocks(j).amplitude > min_amplitude_absolute ...
%         && abs(full(i).clocks(j).deriv2) < max_deriv2...
%         && full(i).clocks(j).frequency > freqrange(1) && full(i).clocks(j).frequency < freqrange(2)
%             promising_clocks = [promising_clocks, full(i).clocks(j)];
%         end
%     end
% end

% % prioritise low second derivative over amplitude
% best_clocks = nestedSortStruct(promising_clocks,{'deriv2mag','amplitude'},[1,-1]);
% 
% fprintf('Top %i clock transitions:\n',num_displayed_clocks)
% %disp(['Min amplitude (relative): ',num2str(min_amplitude_relative)]);
% %disp(['Max 2nd deriv: ',num2str(max_deriv2),' MHz/mT^2']);
% for i = 1:num_displayed_clocks
%     fprintf('f=%.f\tB=%.1f\tamp=%.3f\tderiv2=%.3f\ttransition=%i-->%i\tangle=%.0f\n',...
%                 best_clocks(i).frequency,...
%                 best_clocks(i).field,...
%                 best_clocks(i).amplitude/full_peak_amplitude,...
%                 best_clocks(i).deriv2,...
%                 best_clocks(i).transition(1),...
%                 best_clocks(i).transition(2),...
%                 best_clocks(i).angle/deg...
%     );
% end



save(['full_',datestr(now,'yyyy-mm-ddTHH-MM-SS'),'.mat'],'full')

%% Display freq and second derivative vs field for one transition
% Requires full.data struct output from ClockTransitions

deg = pi/180;
rt = 1.0/sqrt(2);

% IF NEED TO LOAD full STRUCT
filename = ''; % file containing full().data() struct
if ~isempty(filename)
    load(filename);
end

rotation_step = 1;
transition = [7,12]; % specify transition - lower level first
threshold = 0.1;

field_step_size = full(1).data(1).field(2) - full(1).data(1).field(1); % mT
num_of_levels = 0.5*(1+sqrt(1+8*length(full(1).data(:)))); % inverse of 0.5*n*(n-1)
transition_index = sum(num_of_levels-(1:transition(1)))-num_of_levels+transition(2);
transition_str = [int2str(transition(1)),'-->',int2str(transition(2))];

x = [];
y = [];
z = [];

x = [full(rotation_step).data(transition_index).field];
x = x(:); % make column vector
y = [full(rotation_step).data(transition_index).frequency];
y = y(:);
z = [full(rotation_step).data(transition_index).amplitude];
z = z(:);

% Find clock transitions
clocks = findClockTransitions(full(rotation_step).data(transition_index),0.1);
clockfreqs = [clocks.frequency]';
clockfields = [clocks.field]';
clockderiv2 = [clocks.deriv2]';

% Deriv2 for specified transition
deriv1 = diff(y)/field_step_size;
deriv2 = diff(deriv1)/field_step_size; % approximate second derivative
deriv2 = [deriv2(1);deriv2;deriv2(end)]; % MHz/mT^2

figure
subplot(3,1,1)

% Freq for specified transition
hold on
scatter(x,y,[],z,'.')
colormap(flipud(hot))
cbar = colorbar('east');
cbar.Label.String = 'Amplitude';
scatter(clockfields,clockfreqs);
hold off
%caxis([0.0,1.0])
%ylim(Exp.Range*1000)
xlabel('B (mT)')
ylabel('Transition frequency (MHz)')
title('Frequency')
%title(['Mag axis: ',full(1).magaxis,'; MW axis: ',full(1).MWaxis,'; Rot axis: ',full(1).rotaxis,'; transition: ',transition_str])
%saveas(gcf,['transitions',transition_str,'.png'])

subplot(3,1,2)
hold on
scatter(x(1:end-1),deriv1,'.')
xL = get(gca,'XLim');
line(xL,[0 0],'Color',[0.8 0.8 0.8])
%colormap(flipud(hot))
%cbar = colorbar();
%cbar.Label.String = 'Amplitude';
%caxis([0.0,1.0])
scatter(clockfields,repelem(0,length(clockfields)));
hold off
%ylim(Exp.Range*1000)
%ax = gca;
%ax.XAxisLocation = 'origin';
xlabel('B (mT)')
ylabel('df/dB (MHz/mT)')
title('First derivative');
%title(['Mag axis: ',full(1).magaxis,'; MW axis: ',full(1).MWaxis,'; Rot axis: ',full(1).rotaxis,'; transition: ',transition_str])


subplot(3,1,3)

hold on
scatter(x,deriv2,'.')
xL = get(gca,'XLim');
line(xL,[0 0],'Color',[0.8 0.8 0.8])
%colormap(flipud(hot))
%cbar = colorbar();
%cbar.Label.String = 'Amplitude';
%caxis([0.0,1.0])
scatter(clockfields,clockderiv2);
hold off
%ylim(Exp.Range*1000)
%ax = gca;
%ax.XAxisLocation = 'origin';
xlabel('B (mT)')
ylabel('d^2f/dB^2 (MHz/mT^2)')
title('Second derivative');
%title(['Mag axis: ',full(1).magaxis,'; MW axis: ',full(1).MWaxis,'; Rot axis: ',full(1).rotaxis,'; transition: ',transition_str])

text_label = {['Mag axis: ',full(1).magaxis,'; MW axis: ',full(1).MWaxis],['Rot axis: ',full(1).rotaxis,', angle: ',num2str(full(rotation_step).angle)]};
annotation('textbox',[.55 .25 .0 .0],'string',text_label,'FitBoxToText','on');


saveas(gcf,['singletransition',transition_str,'.png'])
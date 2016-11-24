%% Display clock transition freqs vs angle for one transition

deg = pi/180;
rt = 1.0/sqrt(2);

% IF NEED TO LOAD full STRUCT
filename = 'full_2016-11-24T13-46-07.mat'; % file containing full().data() struct
load(filename);

spec = [7,12]; % specify transition - lower level first
num_of_levels = 0.5*(1+sqrt(1+8*length(full(1).data(:)))); % inverse of 0.5*n*(n-1)
transition_index = sum(num_of_levels-(1:spec(1)))-num_of_levels+spec(2);
transition_str = [int2str(spec(1)),'-->',int2str(spec(2))];

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
title(['Mag axis: ',full(1).magaxis,'; MW axis: ',full(1).MWaxis,'; Rot axis: ',full(1).rotaxis,'; transition: ',transition_str])
saveas(gcf,['transitions',transition_str,'.png'])

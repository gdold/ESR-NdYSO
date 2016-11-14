%output[field,Pos,Amp,transition_label,Site]
%transition_label x->y xxyy

selected_transitions = strong_transitions;

% condition = transitions(:,4)~=103;
% transitions(condition,:) = []

% if both sites look at site 1 only
% CURRENTLY BROKEN - can only simulate 1 site anyway
if any([dat(:).site] == 1) && any([dat(:).site] == 2)
    disp('Simulation contains both sites, stripping out site 2...')
    selected_transitions([dat(selected_transitions).site] == 2) = [];
    %match = clock_data(:,5)~=1;
    %clock_data(match,:) = []; % strip out site 2
end

AllClocks = [];

for transition = selected_transitions'
    data = dat(transition).frequency;
    field = dat(transition).field;
    
    %data = data(:,2);
    [Maxima,MaxIndex] = findpeaks(data); % pass this field as 2nd arg
    [Minima,MinIndex] = findpeaks(-data); % to yeild correct field idx
    Minima = -Minima;
    
    ClockFreqs = [Maxima;Minima];
    ClockIndices = [MaxIndex;MinIndex];
    ClockFields = field(ClockIndices);

    
    
    deriv1 = diff(data);
    deriv2 = diff(deriv1); % 2nd derivative, lazy version
    deriv2 = [deriv2(1);deriv2;deriv2(end)];
    % this won't work for max/min at start or end of data
    
    ClockDeriv2 = deriv2(ClockIndices);
    
    if ~isempty(ClockFreqs)
        AllClocks = [AllClocks;ClockFreqs,ClockFields,ClockDeriv2];
    end
    
end

if ~isempty(AllClocks)
    hold on
    scatter(AllClocks(:,2),AllClocks(:,1))
    hold off
end
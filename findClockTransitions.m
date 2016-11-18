function clocks = findClockTransitions(dat,threshold)

%output[field,Pos,Amp,transition_label,Site]
%transition_label x->y xxyy



% condition = transitions(:,4)~=103;
% transitions(condition,:) = []

% if both sites look at site 1 only
% CURRENTLY BROKEN - can only simulate 1 site anyway
bothsites = false;
if any([dat(:).site] == 1) && any([dat(:).site] == 2)
    disp('Clock transitions for site 1 only')
    bothsites = true;
end

AllClocks = [];
clocks = struct();

field_step_size = dat(1).field(2) - dat(1).field(1); % mT


for transition = 1:length(dat)
    if bothsites && (dat(transition).site == 2)
        continue;
    end
    frequency = dat(transition).frequency;
    field = dat(transition).field;
    amplitude = dat(transition).amplitude;
    
    
    %data = data(:,2);
    [Maxima,MaxIndex] = findpeaks(frequency); % pass this field as 2nd arg
    [Minima,MinIndex] = findpeaks(-frequency); % to yeild correct field idx
    Minima = -Minima;
    
    ClockFreqs = [Maxima;Minima];
    ClockIndices = [MaxIndex;MinIndex];
    ClockFields = field(ClockIndices);
    ClockAmplitudes = amplitude(ClockIndices);
    
    
    deriv1 = diff(frequency)/field_step_size;
    deriv2 = diff(deriv1)/field_step_size; % approximate second derivative
    deriv2 = [deriv2(1);deriv2;deriv2(end)];
    % this won't work for max/min at start or end of data
    
    ClockDeriv2 = deriv2(ClockIndices);
    ClockDeriv2mag = abs(ClockDeriv2);
    
    for clock = 1:length(ClockFreqs)
        % Perhaps want to look for degenerate transitions here?
        clocks(end+1).frequency = ClockFreqs(clock);
        clocks(end).field = ClockFields(clock);
        clocks(end).amplitude = ClockAmplitudes(clock);
        clocks(end).deriv2 = ClockDeriv2(clock);
        clocks(end).deriv2mag = ClockDeriv2mag(clock);
        clocks(end).transition = dat(transition).transition;
        clocks(end).index = transition;
        clocks(end).angle = dat.angle;
    end
    
    if ~isempty(ClockFreqs)
        AllClocks = [AllClocks;ClockFreqs,ClockFields,ClockDeriv2];
    end
    
end

% remove empty first entry
clocks(1) = [];

% check there were some clock transitions
if isempty(fieldnames(clocks))
    disp('No clock transitions found')
    return
else
    % find maximum amplitude
    threshold_amplitude = threshold*max([clocks(:).amplitude]);
    strong_transitions = [clocks(:).amplitude] >= threshold_amplitude;
    clockfreqs = [clocks(strong_transitions).frequency];
    clockfreqs = clockfreqs(:);
    clockfields = [clocks(strong_transitions).field];
    clockfields = clockfields(:);
    disp(['Found ',num2str(length(clockfreqs)),' strong transitions'])
end

if isempty(clockfreqs)
    % This should never happen if threshold <=1
    disp(['No clock transitions with relative amplitude >',num2str(threshold)])
else
    hold on
    scatter(clockfields,clockfreqs)
    hold off
end

%clocks = [];

end
%output[field,Pos,Amp,transition_label,Site]
%transition_label x->y xxyy

clock_data = transitions;

% condition = transitions(:,4)~=103;
% transitions(condition,:) = []

% if both sites look at site 1 only
if any(abs(1-clock_data(:,5))<1e-10) && any(abs(2-clock_data(:,5))<1e-10)
    disp('Simulation contains both sites, stripping out site 2...')
    match = clock_data(:,5)~=1;
    clock_data(match,:) = []; % strip out site 2
end

all_transitions = clock_data;

AllMaxima = [];
AllMinima = [];

for transition = unique(clock_data(:,4))'
    data = all_transitions;
    match = data(:,4)~=transition;
    data(match,:) = []; % strip out data that doesn't match transition
    
    %data = data(:,2);
    [Maxima,MaxIndex] = findpeaks(data(:,2)); % pass this freq as 2nd arg
    [Minima,MinIndex] = findpeaks(-data(:,2)); % to yeild correct freq idx
    
    deriv1 = diff(data(:,2));
    deriv2 = [0;diff(deriv1);0]; % 2nd derivative, lazy version
    % this won't work for max/min at start or end of data
    
    MaxDeriv2 = deriv2(MaxIndex);
    MinDeriv2 = deriv2(MinIndex);
    
    %delta = 0.5;
    %[maxtab,mintab] = peakdet(data(:,2),1e-50,data(:,1));
    
    Minima = -Minima;
    
    MaxIndex = data(MaxIndex,1);
    MinIndex = data(MinIndex,1);
    
    if ~isempty(Maxima)
        AllMaxima = [AllMaxima;Maxima,MaxIndex,MaxDeriv2];
    end

    
    if ~isempty(Minima)
        AllMinima = [AllMinima;Minima,MinIndex,MinDeriv2];
    end
    
end



if ~isempty(AllMaxima)
    hold on
    scatter(AllMaxima(:,2),AllMaxima(:,1))%,[],AllMaxima(:,3))
    hold off
end
if ~isempty(AllMinima)
    hold on
    scatter(AllMinima(:,2),AllMinima(:,1))%,[],AllMinima(:,3))
    hold off
end
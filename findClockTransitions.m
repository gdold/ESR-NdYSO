%output[field,Pos,Amp,transition_label,Site]
%transition label x->y xxyy

transitions = output;

% condition = transitions(:,4)~=103;
% transitions(condition,:) = []

% if both sites look at site 1 only
if any(abs(1-transitions(:,5))<1e-10) && any(abs(2-transitions(:,5))<1e-10)
    disp('Simulation contains both sites, stripping out site 2...')
    match = output(:,5)~=1;
    transitions(match,:) = []; % strip out site 2
end

all_transitions = transitions;

AllMaxima = [];
AllMinima = [];

for transition = unique(output(:,4))'
    data = all_transitions;
    match = data(:,4)~=transition;
    data(match,:) = [];
    
    data = data(:,2);
    [Maxima,MaxIndex] = findpeaks(data);
    [Minima,MinIndex] = findpeaks(-data);
    Minima = -Minima;
    
    MaxIndex = MaxIndex*50/500;
    MinIndex = MinIndex*50/500;
    
    if ~isempty(Maxima)
        AllMaxima = [AllMaxima;Maxima,MaxIndex];
    end
    
    if ~isempty(Minima)
        AllMinima = [AllMinima;Minima,MinIndex];
    end
    
end

hold on
scatter(AllMaxima(:,2),AllMaxima(:,1))
scatter(AllMinima(:,2),AllMinima(:,1))
hold off
function strongTransitions = findStrongTransitions(dat,threshold)

    strongTransitions = [];
%    threshold = 0.1;
    for i = 1:length(dat)
        dat(i).peak_amplitude = max(dat(i).amplitude);
    end
    max_amplitude = max([dat(:).peak_amplitude]);
    for i = 1:length(dat)
        if dat(i).peak_amplitude > threshold*max_amplitude
            strongTransitions = [strongTransitions; i];
        end
    end
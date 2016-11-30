function best_clocks = findBestClocks(full,freqrange,max_deriv2,min_amplitude_relative,num_displayed_clocks,require_levels)

deg = pi/180;

promising_clocks = struct([]);

full_peak_amplitude = 0;
for i = 1:length(full) % neater way to do this?
    full_peak_amplitude = max([[full(i).clocks(:).amplitude], full_peak_amplitude]);
end
min_amplitude_absolute = min_amplitude_relative*full_peak_amplitude;



for i = 1:length(full)
    if isempty(fieldnames(full(i).clocks))
        continue
    end
    for j = 1:length(full(i).clocks)
        if full(i).clocks(j).amplitude > min_amplitude_absolute ...
        && abs(full(i).clocks(j).deriv2) < max_deriv2...
        && full(i).clocks(j).frequency > freqrange(1) && full(i).clocks(j).frequency < freqrange(2)
            promising_clocks = [promising_clocks, full(i).clocks(j)];
        end
    end
end

clocks_matching_transition = struct([]);

% remove transitions that don't match require_levels
if isempty(require_levels)%~exist('require_levels','var')
    clocks_matching_transition = promising_clocks;
else
    for i = 1:length(promising_clocks)
        if any(require_levels(1) == promising_clocks(i).transition)...
        && any(require_levels(end) == promising_clocks(i).transition)
            clocks_matching_transition = [clocks_matching_transition, promising_clocks(i)];
        end
    end
end

if isempty(fieldnames(clocks_matching_transition))
    disp('No clock transitions matching criteria');
    return
end

% prioritise low second derivative over amplitude
best_clocks = nestedSortStruct(clocks_matching_transition,{'deriv2mag','amplitude'},[1,-1]);
%best_clocks = nestedSortStruct(clocks_matching_transition,{'amplitude','deriv2mag'},[-1,1]);

num_displayed_clocks = min(num_displayed_clocks,length(best_clocks));

fprintf('Top %i clock transitions:\n',num_displayed_clocks)
%disp(['Min amplitude (relative): ',num2str(min_amplitude_relative)]);
%disp(['Max 2nd deriv: ',num2str(max_deriv2),' MHz/mT^2']);
for i = 1:num_displayed_clocks
    fprintf('f=%.f\tB=%.1f\tamp=%.4f\tderiv2=%.3f\ttransition=%i-->%i\tangle=%.2f\n',...
                best_clocks(i).frequency,...
                best_clocks(i).field,...
                best_clocks(i).amplitude/full_peak_amplitude,...
                best_clocks(i).deriv2,...
                best_clocks(i).transition(1),...
                best_clocks(i).transition(2),...
                best_clocks(i).angle/deg...
    );
end

function [filtered_LFP, clean_events] = pre_process(LFP, events,k,SWRLTDIdx)
load(events)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Per Ripple Analysis %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess signal
%%%%%%

filtered_LFP = BPfilter(LFP, 1250, 30, 60); % isolate gamma frequency band
LTD_events = SWRevents(SWRLTDIdx(k).R, 1:2); % grab SWRs that occur during this period
LTD_events(LTD_events <= 625) = NaN;
LTD_events(LTD_events >= length(LFP)-1250) = NaN;
clean_events = rmmissing(LTD_events, 1);
if isempty(clean_events)
    disp(['No LTD events for animal: ', num2str(cur_animal), ' recording: ', num2str(k)])
end

end
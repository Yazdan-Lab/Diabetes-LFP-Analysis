function [ripple, label] = label_ripples(prev_rip, prev_label, events,idx,rec,c)
ripple = [prev_rip; events(idx(rec).R,:)]; % DB+ 200D
temp_label = zeros(size(events(idx(rec).R,:),1),1);
temp_label(:) = c;
label = [prev_label; temp_label];
end

function [d] = calculate_theta_state(times,d,side)
% This function takes in a struct with two fields, start and end
% and calculates length of state (HTD), ratio of HTD to LTD, and
% number of states

% initialize state length storage variables
htd = [];
ltd = [];

% Calculate state length
for iT = 1:length(times.start)
    % subtract end time (i) from start time (i) to find each length
    % of on state
    htd = [htd (times.end(iT)-times.start(iT))];
    
    % subtract start time (i) from end time (i-1) to find each length
    % of off state. At i = 1, start at 0.
    if iT == 1
        ltd = ltd + times.start(iT);
    else
        ltd = [ltd (times.start(iT) - times.end(iT-1))];
    end
end

perc = sum(htd)/(sum(htd)+sum(ltd));
if length(times.start) == length(times.end)
    s = length(times.start);
else
    disp('Start and end lengths differ');
end
switch side
    case 'L'
        d.lNum = [d.lNum s];
        d.lHtd = [d.lHtd htd];
        d.lLtd = [d.lLtd ltd];
        d.lR   = [d.lR   perc];
    case 'R'
        d.rNum = [d.rNum s];
        d.rHtd = [d.rHtd htd];
        d.rLtd = [d.rLtd ltd];
        d.rR   = [d.rR   perc];
end
end
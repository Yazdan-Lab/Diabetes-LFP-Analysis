function [temp_gamma_pyr, temp_avg_rip, temp_CSD] = process_events(LFP, gamma_LFP, LTD_events, cur_animal, chans, Fs)

% initialize temp storage variables
temp_gamma_pyr = zeros(1, size(LTD_events,1));
temp_avg_rip = zeros(1876, size(LTD_events,1));
temp_CSD = zeros(1876,12, size(LTD_events,1));
for r = 1:size(LTD_events, 1)
    % gamma power of each ripple
    temp_gamma_pyr(r) = SignalPower(gamma_LFP(LTD_events(r, 1):LTD_events(r, 2), chans(2, cur_animal)), 1250);
    temp_avg_rip(:, r) = gamma_LFP(LTD_events(r)-625:LTD_events(r)+1250, chans(2, cur_animal));
    
    % Create CSD of each ripple
    temp_CSD(:, :, r) = CSDlite(LFP(LTD_events(r)-625:LTD_events(r)+1250, chans(2, cur_animal)-5:chans(2, cur_animal)+6), Fs, 1e-4);
    %csd = save_check(csd, temp_CSD);
end

end
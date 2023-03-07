% This notebook will take the dbdb Data and detect ripples using
% thresholding for both pyramidal ripples and radiatum sharp waves
clear all
%



Fs = 1250;
bval = 0;
USER = 'Z';
if USER == 'S'
    addpath('C:\COM\ePhy\dbdb\code\utils-toolbox\utils-toolbox')
    addpath('C:\COM\ePhy\dbdb\code\spectral-analysis-tools')
    data_path = 'C:\COM\ePhy\dbdb\Data\dbdb electrophy';
else
    data_path = 'C:\Users\ipzach\Documents\MATLAB\Data\dbdb electrophy';
    addpath('C:\Users\ipzach\Documents\MATLAB\spectral-analysis-tools')
end
volt_conv_factor  = 0.000000015624999960550667;
smoothing_width = 0.01; % 300 ms
kernel = gaussian(smoothing_width*Fs, ceil(8*smoothing_width*Fs));
cd(data_path)
animals = dir;
%wbar = waitbar(0,'Initializing');
load('SpkInfo.mat'); %MS: Contains
load('PyramChans.mat');
viz = 1;
for i = 1:4
    if i ==1
        grouping = 4:9; % DB+ 200D %Is based on
        message = 'Processing: DB+ 200 ';
    elseif i == 2
        grouping = 10:14; % DB+ 400D
        message = 'Processing: DB+ 400 ';
    elseif i ==3
        grouping = [15:18 20 21]; % DBDB 200D
        message = 'Processing: DBDB 200 ';
    elseif i ==4
        grouping = [22 24:27]; % DBDB 400D
        message = 'Processing: DBDB 400 ';
    end % if i
    group_spwrs = [];
    group_ripples = [];
    for j = grouping
        cd(animals(j).name)
        load('REM.mat'); %ME: loads Theta_Delta info >  figure;plot(rem(1).R.Theta_Delta(:))
        SWR_files = dir('SWR_R_*');
        SWR_files = {SWR_files.name};
        
        LFP_files = dir('LFP*');
        LFP_files = {LFP_files.name};
        bval = bval + (.25/length(grouping));
        %waitbar(bval, wbar, [message num2str(round(bval*100)) '%']);
        for k = 1:size(rem,2)
            if ~isempty(rem(k).R) % makes sure ripple occured during this period
                disp('loading')
                load(char(LFP_files(k)));
                LFP = LFPs{1,2} .*volt_conv_factor;
                % view TD states
                raw_pyr = LFP(:,chans(j));
                theta = BPfilter(raw_pyr, Fs, 4, 7);
                delta = BPfilter(raw_pyr, Fs, 0.1, 3);
                [idx, TD] = calculate_theta_state(theta, delta);
                
                if viz
                   view_data(raw_pyr,TD)
                end
                disp('filtering')
                [pyr, rad] = process_LFP(LFP, j, chans);
                disp('detecting')
                
                
                ripples = detect_events(pyr, 6);
                %waves = detect_events(-rad, 4);
                
                %SPWRs = SPWR_filter(ripples, waves);
                
                %                LTD_SPWRs = TD_check(SPWRs, rem, k);
                LTD_ripples = TD_check(ripples, rem, k);
                
                
                
                % Get intuition for SPWR quality
                if viz
                    view_ripples(LTD_ripples, LFP, chans, j);
                end % if viz
                %group_ripples = cat(1,group_ripples,LTD_ripples);
                %group_spwrs = cat(1,group_spwrs,LTD_SPWRs);
                
            else
                group_spwrs = [];
                group_ripples = [];
            end % if empty
            % save with consistent file name for future use
            SWRevents = group_spwrs;
            save(['SPWR_R_' num2str(k)],'pyr','rad','SWRevents')
        end % For k
        cd ..
    end % for j
    
    if i ==1
        DB200_spwrs = group_spwrs;
    elseif i ==2
        DB400_spwrs = group_spwrs;
    elseif i ==3
        DBDB200_spwrs = group_spwrs;
    elseif i ==4
        DBDB400_spwrs = group_spwrs;
    end % if i
    
    
end % for i
close(wbar)
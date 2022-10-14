% This notebook will take the dbdb Data and detect ripples using
% thresholding for both pyramidal ripples and radiatum sharp waves
clear all
addpath('C:\Users\ipzach\Documents\MATLAB\Toolbox Zach',...
    'C:\Users\ipzach\Documents\MATLAB\spectral-analysis-tools')

load('SpkInfo.mat');
load('PyramChans.mat');
Fs = 1250;
bval = 0;
data_path = 'C:\Users\ipzach\Documents\MATLAB\Data\dbdb electrophy';
volt_conv_factor  = 0.000000015624999960550667;
smoothing_width = 0.01; % 300 ms
kernel = gaussian(smoothing_width*Fs, ceil(8*smoothing_width*Fs));
cd(data_path)
animals = dir;
wbar = waitbar(0,'Initializing');



viz = 1;
for i = 1:4
    if i ==1
        grouping = 3:9; % DB+ 200D
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
        load('REM.mat');
        SWR_files = dir('SWR_R_*');
        SWR_files = {SWR_files.name};
        
        LFP_files = dir('LFP*');
        LFP_files = {LFP_files.name};
        bval = bval + (.25/length(grouping));
        waitbar(bval, wbar, [message num2str(round(bval*100)) '%']);
        for k = 1:size(rem,2)
            if ~isempty(rem(k).R) % makes sure ripple occured during this period
                load(char(LFP_files(k)));
                LFP = LFPs{1,2} .*volt_conv_factor;
                
                [pyr, rad] = process_LFP(LFP, j, chans);
                
                ripples = detect_events(pyr, 4);
                waves = detect_events(-rad, 4);
                
                SPWRs = SPWR_filter(ripples, waves);
                
                LTD_SPWRs = TD_check(SPWRs, rem, k);
                LTD_ripples = TD_check(ripples, rem, k);
                % Get intuition for SPWR quality
                if viz
                    figure
                    for viz = 1:size(LTD_SPWRs,1)
                        ripple_visualize(LTD_SPWRs(viz,1),...
                            LTD_SPWRs(viz,2),...
                            LFP,...
                            chans,...
                            j);
                        pause(1)
                    end % for ripple viz
                end % if viz
                group_ripples = cat(1,group_ripples,LTD_ripples);
                group_spwrs = cat(1,group_spwrs,LTD_SPWRs);
                
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
% This notebook will take the dbdb Data and detect ripples using
% thresholding for both pyramidal ripples and radiatum sharp waves
clear all
addpath('C:\Users\ipzach\Documents\MATLAB\Toolbox Zach',...
    'C:\Users\ipzach\Documents\MATLAB\spectral-analysis-tools')

load('SpkInfo.mat');
load('PyramChans.mat');
Fs = 1250;
data_path = 'C:\Users\ipzach\Documents\MATLAB\Data\dbdb electrophy';
volt_conv_factor  = 0.000000015624999960550667;

smoothing_width = 0.01; % 300 ms
kernel = gaussian(smoothing_width*Fs, ceil(8*smoothing_width*Fs));

% SpwrStats
% SWR Occurance
% SPW Occurance
% Overlap Pct
% SWR During LTD
% SPW During LTD
% Overlap during LTD
cd(data_path)
animals = dir;

for i = 1:4
    disp(num2str(i))
    if i ==1
        grouping = 3:9; % DB+ 200D
    elseif i ==2
        grouping = 10:14; % DB+ 400D
    elseif i ==3
        grouping = [15:18 20 21]; % DBDB 200D
    elseif i ==4
        grouping = [22 24:27]; % DBDB 400D
    end % if i
    group_spwrs = [];
    for j = grouping
        cd(animals(j).name)
        load('REM.mat');
        SWR_files = dir('SWR_R_*');
        SWR_files = {SWR_files.name};
        
        LFP_files = dir('LFP*');
        LFP_files = {LFP_files.name};
        
        for k = 1:size(rem,2)
            if ~isempty(rem(k).R) % makes sure ripple occured during this period
                load(char(LFP_files(k)));
                disp('Detecting')
                LFP = LFPs{1,2} .*volt_conv_factor;
                
                [pyr, rad] = process_LFP(LFP, j, chans);
                
                ripples = detect_events(pyr, 5);
                waves = detect_events(-rad, 3);
                
                SPWRs = SPWR_filter(ripples, waves);
                
                LTD_SPWRs = TD_check(SPWRs, rem, k);
                
                group_spwrs = cat(1,group_spwrs,LTD_SPWRs);
            end % if empty
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
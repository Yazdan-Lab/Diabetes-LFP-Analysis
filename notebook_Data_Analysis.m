% Last edited by Zachary Ip 3/7/23
% This script loads in LFP files from folder structure and saves isolated
% measures of the LFP for further statistical processing.
% This code was designed for use in both the UW NERD lab and Jialing Liu's
% UCSF, as designated by `USER` variable. 


%% Global variables
% For LFP processing
VOLTCONV = 0.000000091555527603759401; % Neurolynx saves data in a unitless value, we need this to convert it to volts
FS = 1250;

% For signal smoothing (Ripple detection)
KERNEL = gaussian(FS, ceil(8*FS));

% For CSD, after centering on pyramidal channel
HIGH_CHAN = 7; %pyramidal channel
LOW_CHAN = 11; % Radiatum channel
PRE_WIN = 1:550; % pre indices
WIN = 650:750; % ripple indices
POST_WIN = 850:1300; % post indices
%%%%%%
USER = 'Z'; %'Z' for Zach or 'S' for Shahram for path stuff

switch USER
    case 'Z'
        cd('C:\Users\ipzach\Documents\MATLAB\Diabetes-Data-Analysis')
        addpath('C:\Users\ipzach\Documents\MATLAB\Toolbox Zach', ...
            'C:\Users\ipzach\Documents\MATLAB\spectral-analysis-tools',...
            'C:\Users\ipzach\Documents\MATLAB\Diabetes-Data-Analysis')
    case 'S'
        cd('C:\COM\ePhy\dbdb\code\Diabetes-LFP-Analysis')
end

% Load reference 
load('SpkInfo.mat')
load('chans.mat')
switch USER
    case 'Z'
        cd('C:\Users\ipzach\Documents\MATLAB\Data\dbdb electrophy')
    case 'S'
        cd('C:\COM\ePhy\dbdb\Data\dbdb electrophy'); % here is the data
end

% Find animal LFP Data
animal_list = dir; % create a list of every folder (each one is one animal)

% Initialize storage variable for final measures
rip.DB2 = []; rip.DB4 = []; rip.DBDB2 = []; rip.DBDB4 = [];
label.DB2 = []; label.DB4 = []; label.DBDB2 = []; label.DBDB4 = [];
intSlo_Store = cell(4, 7); intSlo0_Store = cell(4, 7);
Co = NaN(4, 7, 7, 3); PLI = NaN(4, 7, 7, 3);
powers = NaN(4,7,7,3);
GR_powers = NaN(4,7,7,3);
% Data structure Group, Animal, freq_band, Layer_comb
slowing_score = NaN(4, 7, 3);
state_changes = NaN(4, 7);
% Group,Animal,Layer_comb

% Begin parsing data
%%%%%%
counter = 0;
% Separate animals by group first
for group = 1:4
    % Grab indices of animals in a particular group
    if group == 1
        grouping = 3:9; % DB+ 200D #7
        gn = 'DB+ 200D';
    elseif group == 2
        grouping = 10:14; % DB+ 400D #5
        gn = 'DB+ 400D';
    elseif group == 3
        grouping = [15:18, 20, 21]; % DBDB 200D #6
        gn = 'DBDB 200D';
    elseif group == 4
        grouping = [22, 24:27]; % DBDB 400D #5
        gn = 'DBDB 400D';
    end
    % start iterator since group number is based on place in directory
    counter = 0;
    
    % Initialize temp storage variables for group
    gamma_pyr = [];
    gamma_slm = [];
    gamma_ctx = [];
    avg_rip = [];
    csd = [];
    dur = [];
    iri = [];
    % initialize struct for data validation and saving for UCSF 
    single_animal_measures.gamma = [];
    single_animal_measures.CSD = [];
    single_animal_measures.CSD_pre = [];
    single_animal_measures.CSD_post = [];
    single_animal_measures.events = [];
    single_animal_measures.events = [];
    single_animal_measures.dur = [];
    single_animal_measures.IRI = [];
    
    % Start Processing single animals
    for cur_animal = grouping
        disp(['Animal: ', num2str(cur_animal)])
        [SWR_files, LFP_files, SWRLTDIdx, rem] = initalize_animal(cur_animal, animal_list);
       
        counter = counter + 1;
        full_LFP = [];
        state_changes(group, counter) = 0;
        % Each animal has multiple sequential recordings, they must first
        % be concatenated for full LFP measures (signal power, coherence,
        % etc), however, ripple timings are saved per recording, so they
        % are processed first, and then added together
        for k = 1:size(SWRLTDIdx,2) % run through each of the LTD periods (where SWRs occur)
            if ~isempty(rem(k).R) % Make sure recording exists
                state_changes(group, counter) = state_changes(group, counter) + length(rem(k).R.start);
            else
                disp(['Rem file missing for animal ', num2str(cur_animal)])
            end % if REM not empty
            load(char(LFP_files(k))); % Load LFP events
            
            %create consistent LFP variable
            LFP = LFPs{1, 2} .* VOLTCONV *1000; % load LFP, convert to mV
            full_LFP = [full_LFP; LFP];
            
            % To detect ripples, must get theta state index first, and look
            % during non-theta. Theta states were recorded by Gratianne and
            % given as part of meta data for each animal
            if ~isempty(SWRLTDIdx(k).R) % makes sure ripple occured during this period
                % Load in animal data
                %%%%%%
                [gamma_LFP, LTD_events] = pre_process(LFP, char(SWR_files(k)), k, SWRLTDIdx); % Load SWR events
                
                [temp_gamma_pyr, temp_avg_rip, temp_CSD] = process_events(LFP,gamma_LFP, LTD_events, cur_animal,chans, FS);
               
                
                temp_dipole_pre = calculate_CSD_dipole(temp_CSD, HIGH_CHAN, LOW_CHAN, PRE_WIN)';
                temp_dipole_post = calculate_CSD_dipole(temp_CSD, HIGH_CHAN, LOW_CHAN, POST_WIN)';
                temp_dipole = calculate_CSD_dipole(temp_CSD, HIGH_CHAN, LOW_CHAN, WIN)';
                
%                 single_animal_measures(counter).CSD_pre = [single_animal_measures(counter).CSD_pre, temp_dipole_pre];
%                 single_animal_measures(counter).CSD = [single_animal_measures(counter).CSD, temp_dipole];
%                 single_animal_measures(counter).CSD_post = [single_animal_measures(counter).CSD_post, temp_dipole_post];
%                 single_animal_measures(counter).gamma = [single_animal_measures(counter).gamma, temp_gamma_pyr];
                % Concatenate temp variable to storage variable
                %gamma_ctx = [gamma_ctx temp_gamma_ctx];
                gamma_pyr = [gamma_pyr, temp_gamma_pyr];
                %gamma_slm = [gamma_slm temp_gamma_slm];
                avg_rip = [avg_rip, temp_avg_rip];
                csd = cat(3, csd, temp_CSD);
%                 single_animal_measures(counter).events = [single_animal_measures(counter).events; LTD_events];
%                 single_animal_measures(counter).dur = [single_animal_measures(counter).dur; LTD_events(:, 2) - LTD_events(:, 1)];
%                 single_animal_measures(counter).IRI = [single_animal_measures(counter).IRI; diff(LTD_events(:, 1))];
                iri = [iri; diff(LTD_events(:, 1))];
                % Save all individual ripples for basic analysis
                if group == 1
                    [rip.DB2, label.DB2] = label_ripples(rip.DB2, label.DB2, LTD_events, counter);
                elseif group == 2
                    [rip.DB4, label.DB4] = label_ripples(rip.DB4, label.DB4, LTD_events, counter);
                elseif group == 3
                    [rip.DBDB2, label.DBDB2] = label_ripples(rip.DBDB2, label.DBDB2, LTD_events, counter);
                elseif group == 4
                    [rip.DBDB4, label.DBDB4] = label_ripples(rip.DBDB4, label.DBDB4, LTD_events, counter);
                end
            else
                disp(['SWRLTDIdx is empty/missing for animal ', num2str(cur_animal) ' recording: ' num2str(k)])
            end %is empty SWRLTD
            
        end % for k SWRLTDIdx
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Per animal analysis %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(full_LFP)
            for layer = 1:3
                %%% Slowing score
                %%%%%%%%%%%%%%%%%
                % BP filter LFP for low and high frequencies
                single_channel = full_LFP(:, chans(layer, cur_animal));
                low_freq = BPfilter(single_channel, 1250, 1, 8);
                high_freq = BPfilter(single_channel, 1250, 9, 30);
                
                % Calculate ratio of signal powers
                slowing_score(group, counter, layer) = SignalPower(low_freq, 1250) ./ SignalPower(high_freq, 1250);
%                 single_animal_measures(counter).SS(layer) = slowing_score(group, counter, layer);
                % Group, Band, Animal, Layer
                
                % calculate Spectral exponent
                [PSD, frex] = pwelch(single_channel, 2*FS, FS, [], FS);
                frBand = [1, 40];
                frBins = dsearchn(frex, frBand(1)):dsearchn(frex, frBand(2));
                XX = frex(frBins);
                YY = PSD(frBins);
                robRegMeth = 'ols';
                doPlot = 01; %figure
                thisCol = [0, 0, 1];
                [intSlo, stat, Pows, Deviants, stat0, intSlo0] = fitPowerLaw3steps(XX, YY, robRegMeth, doPlot, thisCol);
                Pows_store(group, counter,layer) = Pows;
                intSlo_Store{group, counter, layer} = intSlo;
                intSlo0_Store{group, counter, layer} = intSlo0;
%                 single_animal_measures(counter).SE(layer) = intSlo0(2);
                % Coherence and PLI
                %%%%%%%%%%%%%
                switch layer
                    case 1
                        A = 1;
                        B = 2;
                        compare = 'Cortex-Pyr';
                    case 2
                        A = 1;
                        B = 3;
                        compare = 'Cortex-Slm';
                    case 3
                        A = 2;
                        B = 3;
                        compare = 'Pyr-Slm';
                end % switch layer
                disp(compare)
                % Assign single channels to run coherence on
                
                A_LFP = full_LFP(:, chans(A, cur_animal));
                B_LFP = full_LFP(:, chans(B, cur_animal));
                % create a vector of indiviudal frequencies to calculate
                % coherence
                for freq_band = 1:7
                    switch freq_band
                        case 1
                            lower = 0.1;    upper = 3;
                        case 2
                            lower = 4;      upper = 7;
                        case 3
                            lower = 8;      upper = 13;
                        case 4
                            lower = 13;     upper = 30;
                        case 5
                            lower = 30;     upper = 58;
                        case 6
                            lower = 62;     upper = 200;
                        case 7
                            lower = 0.1;    upper = 200;
                    end % switch iBand
                    
                    range = linspace(lower, upper, 20);
                    A_filt = BPfilter(A_LFP, 1250, lower, upper);
                    B_filt = BPfilter(B_LFP, 1250, lower, upper);
                    % Run coherence and average outputs for each frequency
                    % band
                    % Group, Band, recording, Animal, Layer/layer
                    % Co(1: group id, animal, frequency band, layer
                    % comparison (check compare variable))
                    % comparison)
                    Co(group, counter, freq_band, layer) = nanmean(mscohere(A_LFP, B_LFP, hamming(12500), [], range, 1250));
                    PLI(group, counter, freq_band, layer) = PLV(angle(hilbert(A_filt))', angle(hilbert(B_filt))');
                    powers(group, counter, freq_band, layer) = SignalPower(BPfilter(full_LFP(:,chans(layer, cur_animal)), FS, lower, upper),FS);
                    
                    filtered = BPfilter(full_LFP(:,chans(layer, cur_animal)), FS, lower, upper);
                    normalized = ((filtered - mean(filtered))/std(filtered));
                    envelope = abs(hilbert(filtered));
                    smoothed = smoothvect(envelope, KERNEL);
                    GR_powers(group, counter, freq_band, layer) = SignalPower(smoothed, FS);
                    
%                     single_animal_measures(counter).Co(layer, freq_band) = Co(group, counter, freq_band, layer);
%                     single_animal_measures(counter).PLI(layer, freq_band) = PLI(group, counter, freq_band, layer);
%                     single_animal_measures(counter).LFP = full_LFP;
                    
                end % frequency band
            end % layer
        else
            disp(['LFP is empty for animal ', num2str(cur_animal)])
        end % if
        cd ..
        %UCSF_single_animal_excel(single_animal_measures(COUNT), [gn, '_', num2str(counter)], user);
        
        
        %save(['C:\COM\ePhy\dbdb\Data\Outputs\Data\SingleAnimal\' gn '_' num2str(counter) ],'single_animal_measures(COUNT)')
    end % animal
    % Save data to variable outside loop
    if group == 1
        %Gamma.DB2_Ctx = gamma_ctx;
        Gamma.DB2_Pyr = gamma_pyr;
        %Gamma.DB2_SLM = gamma_slm;
        CSD.DB2 = csd;
        rip_wav.DB2 = avg_rip;
        iri_all.DB2 = iri;
    elseif group == 2
        %Gamma.DB4_Ctx = gamma_ctx;
        Gamma.DB4_Pyr = gamma_pyr;
        %Gamma.DB4_SLM = gamma_slm;
        CSD.DB4 = csd;
        rip_wav.DB4 = avg_rip;
        iri_all.DB4 = iri;
    elseif group == 3
        %Gamma.DBDB2_Ctx = gamma_ctx;
        Gamma.DBDB2_Pyr = gamma_pyr;
        %Gamma.DBDB2_SLM = gamma_slm;
        CSD.DBDB2 = csd;
        rip_wav.DBDB2 = avg_rip;
        iri_all.DBDB2 = iri;
    elseif group == 4
        %Gamma.DBDB4_Ctx = gamma_ctx;
        Gamma.DBDB4_Pyr = gamma_pyr;
        %Gamma.DBDB4_SLM = gamma_slm;
        CSD.DBDB4 = csd;
        rip_wav.DBDB4 = avg_rip;
        iri_all.DBDB4 = iri;
    end % if
end % group
% Group, Animal, freq_band, Layer

%% save processed data
switch USER
    case 'Z'
        cd('C:\Users\ipzach\Documents\MATLAB\output\Diabetes-Saved-Files')
    case 'S'
        cd('C:\COM\ePhy\dbdb\Data\Outputs\Data')
end
save('LFP Measures', 'Gamma', 'rip', 'iri_all', 'rip_wav', 'label', 'CSD', 'Co', 'PLI', 'slowing_score', 'state_changes', 'intSlo0_Store', 'intSlo_Store', 'Pows_store','powers', 'GR_powers')

%% Create excel files
ripple_measures = {'gamma','CSD','dur','IRI','CSD_pre','CSD_post'};


animal_num = 23;
% ripple measures

for i = 1:length(ripple_measures)
    storage = NaN(995,animal_num);
    disp(['i: ' num2str(i)])
    for j = 1:animal_num
        disp(['j: ' num2str(j)])
        switch i
            case 1
                storage(1:length(single_animal_measures(j).gamma),j) = single_animal_measures(j).gamma;
                disp('test')
            case 2
                storage(1:length(single_animal_measures(j).CSD),j) = single_animal_measures(j).CSD;
            case 3
                storage(1:length(single_animal_measures(j).dur),j) = single_animal_measures(j).dur;
            case 4
                storage(1:length(single_animal_measures(j).IRI),j) = single_animal_measures(j).IRI;
            case 5
                storage(1:length(single_animal_measures(j).CSD_pre),j) = single_animal_measures(j).CSD_pre;
            case 6
                storage(1:length(single_animal_measures(j).CSD_post),j) = single_animal_measures(j).CSD_post;
        end % switch
        
    end % animal num
    
    return_path = pwd;
    Datetime = string(datetime('now'));
    if strcmp(USER,'S')
        cd('C:\COM\ePhy\dbdb\Data\Outputs\Data\SingleAnimal');
    else
        cd('C:\Users\ipzach\Documents\MATLAB\output\excel');
    end
    Filename = sprintf([ripple_measures{i} '_data_%s.xlsx'], Datetime);
    Filename = regexprep(Filename, {' ', ':', '-'}, {'_', '_', '_'});
    xlswrite(Filename, storage);
    cd(return_path)
end % ripple measure

%% single measures
single_measures = {'SS','SE'};
for i = 1:length(single_measures)
    storage = NaN(1,animal_num);
    for layer = 1:3
        switch layer
            case 1
                layer_n = '_cortex';
            case 2
                layer_n = '_pyr';
            case 3
                layer_n = '_slm';
        end %switch
        for j = 1:animal_num
            switch i
                case 1
                    storage(j) = single_animal_measures(j).SS(layer);
                case 2
                    storage(j) = single_animal_measures(j).SE(layer);
            end % switch
            
        end % animal num
    return_path = pwd;
    Datetime = string(datetime('now'));
    if strcmp(USER,'S')
        cd('C:\COM\ePhy\dbdb\Data\Outputs\Data\SingleAnimal');
    else
        cd('C:\Users\ipzach\Documents\MATLAB\output\excel');
    end
    Filename = sprintf([single_measures{i} layer_n '_data_%s.xlsx'], Datetime);
    Filename = regexprep(Filename, {' ', ':', '-'}, {'_', '_', '_'});
    xlswrite(Filename, storage);
    cd(return_path)
    end
end % single measures
block_measures = {'Co','PLI'};
%% block measures
for i = 1:length(block_measures)
    
    for freq = 1:7
        switch freq
            case 1
                freq_n = '_delta';
            case 2
                freq_n = '_theta';
            case 3
                freq_n = '_alpha';
            case 4
                freq_n = '_beta';
            case 5
                freq_n = '_gamma';
            case 6 
                freq_n = '_hgamma';
            case 7
                freq_n = '_broadband';
        end
        for layer = 1:3
            storage = NaN(1,animal_num);
            switch layer
                case 1
                    layer_n = '_ctx_pyr';
                case 2
                    layer_n = '_ctx_slm';
                case 3
                    layer_n = '_pyr_slm';
            end
            for j = 1:animal_num
                switch i
                    case 1
                        storage(j) = single_animal_measures(j).Co(layer,freq);
                    case 2
                        storage(j) = single_animal_measures(j).PLI(layer,freq);
                end % switch
                
            end % animal num
            
            return_path = pwd;
            Datetime = string(datetime('now'));
            if strcmp(USER,'S')
                cd('C:\COM\ePhy\dbdb\Data\Outputs\Data\SingleAnimal');
            else
                cd('C:\Users\ipzach\Documents\MATLAB\output\excel');
            end
            Filename = sprintf([block_measures{i} layer_n freq_n '_data_%s.xlsx'], Datetime);
            Filename = regexprep(Filename, {' ', ':', '-'}, {'_', '_', '_'});
            xlswrite(Filename, storage);
            cd(return_path)
            
            
            
            
            
        end %layer
        
    end %freq
    
    
end % block measures

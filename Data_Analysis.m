% Load in initializing variables
voltConv = 0.000000091555527603759401; % Neurolynx saves data in a unitless value, we need this to convert it to volts
Fs = 1250;
%%%%%%
cd('C:\Users\ipzach\Documents\dbdb electrophy\Diabetes-Data-Analysis')
load('SpkInfo.mat')
load('chans.mat')
cd('C:\Users\ipzach\Documents\dbdb electrophy'); % here is the data
animal_list = dir; % create a list of every folder (each one is one animal)


rip.DB2 = [];
rip.DB4 = [];
rip.DBDB2 = [];
rip.DBDB4 = [];

label.DB2 = [];
label.DB4 = [];
label.DBDB2 = [];
label.DBDB4 = [];

Co = NaN(4,7,7,3);
% Group, Animal, freq_band, Layer_comb
slowing_score = NaN(4,7,3);
% Group,Animal,Layer_comb

% Begin parsing data
%%%%%%
for group = 1:4
    % Grab indices of animals in a particular group
    if group ==1
        grouping = 3:9; % DB+ 200D
    elseif group ==2
        grouping = 10:14; % DB+ 400D
    elseif group ==3
        grouping = [15:18 20 21]; % DBDB 200D
    elseif group ==4
        grouping = [22 24:27]; % DBDB 400D
    end
    counter = 0;
    
    spec_pyr = [];
    spec_slm = [];
    spec_ctx = [];
    csd = [];
    
    for cur_animal = grouping
        disp(['Animal: ' num2str(cur_animal)])
        cd(animal_list(cur_animal).name)
        % Get file names of specific animals relevant files
        %%%%%%
        load('SWR_Index.mat'); % load SWR HTD/LTD data
        SWR_files = dir('SWR_R_*'); % Grab the number of files that have SWR event timings (usually 2, sometimes 1 or 3)
        SWR_files = {SWR_files.name}; % throw away useless info
        
        LFP_files = dir('LFP*'); % grab number of LFP files (sometimes 1 or 3 as well)
        LFP_files = {LFP_files.name}; % throw away useless info
        
        counter = counter +1;
        full_LFP = [];
        % Check Data fidelity
        %%%%%%
        for k = 1:length(SWRLTDIdx) % run through each of the LTD periods (where SWRs occur)
            if ~isempty(SWRLTDIdx(k).R) % makes sure ripple occured during this period
                % Load in animal data
                %%%%%%
                load(char(SWR_files(k)));  % Load SWR events
                load(char(LFP_files(k))); % Load LFP events
                LFP = LFPs{1,2} .*voltConv; % load LFP
                full_LFP = [full_LFP; LFP];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Per Ripple Analysis %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Preprocess signal
                %%%%%%
                gamma_LFP = BPfilter(LFP,1250,30,60); % isolate gamma frequency band
                LTD_events = SWRevents(SWRLTDIdx(k).R,1); % grab SWRs that occur during this period
                LTD_events = LTD_events(LTD_events >= 625); % Exclude SWRs that clip beginning window
                LTD_events = LTD_events(LTD_events <= length(LFP)-1250);% Exclude SWRs that exclude final window
                
                % initialize temp storage variables
                temp_spec_ctx = zeros(50,28,length(LTD_events));
                temp_spec_pyr = zeros(50,28,length(LTD_events));
                temp_spec_slm = zeros(50,28,length(LTD_events));
                
                temp_CSD = zeros(1876,12,length(LTD_events));
                % For each event, create a spectrogram and store it in the
                % temp variable
                for r = 1:length(LTD_events)
                    % Create spectrogram of each ripple
                    [~,~,~,temp_spec_ctx(:,:,r)] = spectrogram(gamma_LFP(LTD_events(r)-625:LTD_events(r)+1250,chans(1,cur_animal)),hamming(125),[],[5:5:250],1250);
                    [~,~,~,temp_spec_pyr(:,:,r)] = spectrogram(gamma_LFP(LTD_events(r)-625:LTD_events(r)+1250,chans(2,cur_animal)),hamming(125),[],[5:5:250],1250);
                    [~,freqs,time,temp_spec_slm(:,:,r)] = spectrogram(gamma_LFP(LTD_events(r)-625:LTD_events(r)+1250,chans(3,cur_animal)),hamming(125),[],[5:5:250],1250);
                    
                    % Create CSD of each ripple
                    temp_CSD(:,:,r) = CSDlite(LFP(LTD_events(r)-625:LTD_events(r)+1250,chans(2,cur_animal)-5:chans(2,cur_animal)+6),Fs,1e-4);
                end
                % Isolate Spectral Power
                Ctx_spec_power = squeeze(mean(temp_spec_ctx(5:12,:,:),[1,2]));
                Pyr_spec_power = squeeze(mean(temp_spec_pyr(5:12,:,:),[1,2]));
                Slm_spec_power = squeeze(mean(temp_spec_slm(5:12,:,:),[1,2]));
                % Concatenate temp variable to storage variable
                spec_ctx = [spec_ctx; Ctx_spec_power];
                spec_pyr = [spec_pyr; Pyr_spec_power];
                spec_slm = [spec_slm; Slm_spec_power];
                
                csd = save_check(csd,temp_CSD);
                
                % Save all individual ripples for basic analysis
                if group ==1
                    [rip.DB2,label.DB2] = label_ripples(rip.DB2,label.DB2,SWRevents,SWRLTDIdx,k,counter);
                elseif group ==2
                    [rip.DB4,label.DB4] = label_ripples(rip.DB4,label.DB4,SWRevents,SWRLTDIdx,k,counter);
                elseif group ==3
                    [rip.DBDB2,label.DBDB2] = label_ripples(rip.DBDB2,label.DBDB2,SWRevents,SWRLTDIdx,k,counter);
                elseif group ==4
                    [rip.DBDB4,label.DBDB4] = label_ripples(rip.DBDB4,label.DBDB4,SWRevents,SWRLTDIdx,k,counter);
                end
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
                single_channel = full_LFP(:, chans(layer,cur_animal));
                low_freq = BPfilter(single_channel, 1250, 1, 8);
                high_freq = BPfilter(single_channel, 1250, 9, 30);
                
                % Calculate ratio of signal powers
                slowing_score(group,counter,layer) = SignalPower(low_freq,1250) ./ SignalPower(high_freq,1250);
                % Group, Band, Animal, Layer
                
                % Coherence
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
                
                A_LFP = full_LFP(:, chans(A,cur_animal));
                B_LFP = full_LFP(:, chans(B,cur_animal));
                % create a vector of indiviudal frequencies to calculate
                % coherence
                for freq_band = 1:7
                    switch freq_band
                        case 1
                            range = linspace(0.1,3,20);
                        case 2
                            range = linspace(4,7,20);
                        case 3
                            range = linspace(8,13,20);
                        case 4
                            range = linspace(13,30,20);
                        case 5
                            range = linspace(30,58,20);
                        case 6
                            range = linspace(62,200,20);
                        case 7
                            range = linspace(0,200,50);
                    end % switch iBand
                    % Run coherence and average outputs for each frequency
                    % band
                    % Group, Band, recording, Animal, Layer/layer
                    Co(group,counter,freq_band,layer) = nanmean(mscohere(A_LFP,B_LFP,hamming(12500),[],range,1250));
                end % frequency band
            end % layer
        end % if
        cd ..
        
    end % animal
    % Save data to variable outside loop
    if group ==1
        Spec.DB2_Ctx = spec_ctx;
        Spec.DB2_Pyr = spec_pyr;
        Spec.DB2_SLM = spec_slm;
        CSD.DB2 = csd;
    elseif group ==2
        Spec.DB4_Ctx = spec_ctx;
        Spec.DB4_Pyr = spec_pyr;
        Spec.DB4_SLM = spec_slm;
        CSD.DB4 = csd;
    elseif group ==3
        Spec.DBDB2_Ctx = spec_ctx;
        Spec.DBDB2_Pyr = spec_pyr;
        Spec.DBDB2_SLM = spec_slm;
        CSD.DBDB2 = csd;
    elseif group ==4
        Spec.DBDB4_Ctx = spec_ctx;
        Spec.DBDB4_Pyr = spec_pyr;
        Spec.DBDB4_SLM = spec_slm;
        CSD.DBDB4 = csd;
    end % if
end % group
time = time -0.5; % adjust for SWR event initiation
% Group, Animal, freq_band, Layer
%% save processed data
cd('C:\Users\ipzach\Documents\dbdb electrophy\Diabetes-Saved-Files')
save('LFP Measures','Spec','rip','label','CSD','Co','slowing_score')


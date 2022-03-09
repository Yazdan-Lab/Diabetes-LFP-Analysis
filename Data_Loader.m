




% Load in initializing variables
%%%%%%
cd('C:\Users\ipzach\Documents\dbdb electrophy\Diabetes-Data-Analysis')
load('SpkInfo.mat')
load('chans.mat')
%load('PyramChans.mat')

cd('C:\Users\ipzach\Documents\dbdb electrophy'); % here is the data
animal_list = dir; % create a list of every folder (each one is one animal)

voltConv = 0.000000091555527603759401; % Neurolynx saves data in a unitless value, we need this to convert it to volts
Fs = 1250;
rip.DB2 = [];
rip.DB4 = [];
rip.DBDB2 = [];
rip.DBDB4 = [];

label.DB2 = [];
label.DB4 = [];
label.DBDB2 = [];
label.DBDB4 = [];

% Group,Animal, recording, Layer/layer
scores = NaN(4,7,2,3);

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
    
    SpecPyr = [];
    SpecSLM = [];
    SpecCtx = [];
    Csd = [];
    
    for cur_animal = grouping
        cd(animal_list(cur_animal).name)
        % Get file names of specific animals relevant files
        %%%%%%
        load('SWR_Index.mat'); % load SWR HTD/LTD data
        SWR_files = dir('SWR_R_*'); % Grab the number of files that have SWR event timings (usually 2, sometimes 1 or 3)
        SWR_files = {SWR_files.name}; % throw away useless info
        
        LFP_files = dir('LFP*'); % grab number of LFP files (sometimes 1 or 3 as well)
        LFP_files = {LFP_files.name}; % throw away useless info
        
        counter = counter +1;
        for layer = 1:3
            single_channel = LFP(:, chans(layer,cur_animal));
            low_freq = BPfilter(single_channel, 1250, 1, 8);
            high_freq = BPfilter(single_channel, 1250, 9, 30);
            
            % Group, Band, recording, Animal, Layer/layer
            scores(group,count,k,layer) = SignalPower( low_freq,1250) ./ SignalPower(high_freq,1250);
        end % layer
        % Check Data fidelity
        %%%%%%
        for k = 1:size(SWRLTDIdx,2) % run through each of the LTD periods (where SWRs occur)
            if ~isempty(SWRLTDIdx(k).R) % makes sure ripple occured during this period
                % Load in animal data
                %%%%%%
                load(char(SWR_files(k)));  % Load SWR events
                load(char(LFP_files(k))); % Load LFP events
                LFP = LFPs{1,2} .*voltConv; % load LFP
                
                % Preprocess signal
                %%%%%%
                LFP = BPfilter(LFP,1250,30,60); % isolate gamma frequency band
                LTDevents = SWRevents(SWRLTDIdx(k).R,1); % grab SWRs that occur during this period
                LTDevents = LTDevents(LTDevents >= 625); % Exclude SWRs that clip beginning window
                LTDevents = LTDevents(LTDevents <= length(LFP)-1250);% Exclude SWRs that exclude final window
                
                % initialize temp storage variables
                tempSpecCtx = zeros(50,28,length(LTDevents));
                tempSpecPyr = zeros(50,28,length(LTDevents));
                tempSpecSlm = zeros(50,28,length(LTDevents));
                
                temp_CSD = zeros(1876,12,length(LTDevents));
                % For each event, create a spectrogram and store it in the
                % temp variable
                for r = 1:length(LTDevents)
                    [~,~,~,tempSpecCtx(:,:,r)] = spectrogram(LFP(LTDevents(r)-625:LTDevents(r)+1250,chans(1,cur_animal)),hamming(125),[],[5:5:250],1250);
                    [~,~,~,tempSpecPyr(:,:,r)] = spectrogram(LFP(LTDevents(r)-625:LTDevents(r)+1250,chans(2,cur_animal)),hamming(125),[],[5:5:250],1250);
                    [~,freqs,time,tempSpecSlm(:,:,r)] = spectrogram(LFP(LTDevents(r)-625:LTDevents(r)+1250,chans(3,cur_animal)),hamming(125),[],[5:5:250],1250);
                    
                    temp_CSD(:,:,r) = CSDlite(LFP(LTDevents(r)-625:LTDevents(r)+1250,chans(2,cur_animal)-5:chans(2,cur_animal)+6),Fs,1e-4);
                end
                % Isolate Spectral Power
                Ctx_spec_power = squeeze(mean(mean(tempSpecCtx(5:12,:,:),2),1));
                Pyr_spec_power = squeeze(mean(mean(tempSpecPyr(5:12,:,:),2),1));
                Slm_spec_power = squeeze(mean(mean(tempSpecSlm(5:12,:,:),2),1));
                % Concatenate temp variable to storage variable
                SpecCtx = [SpecCtx; Ctx_spec_power];
                SpecPyr = [SpecPyr; Pyr_spec_power];
                SpecSLM = [SpecSLM; Slm_spec_power];
                
                Csd = save_check(Csd,temp_CSD);
                
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
        cd ..
        
    end % animal
    if group ==1
        Spec.DB2Ctx = SpecCtx;
        Spec.DB2Pyr = SpecPyr;
        Spec.DB2SLM = SpecSLM;
        CSD.DB2 = Csd;
    elseif group ==2
        Spec.DB4Ctx = SpecCtx;
        Spec.DB4Pyr = SpecPyr;
        Spec.DB4SLM = SpecSLM;
        CSD.DB4 = Csd;
    elseif group ==3
        Spec.DBDB2Ctx = SpecCtx;
        Spec.DBDB2Pyr = SpecPyr;
        Spec.DBDB2SLM = SpecSLM;
        CSD.DBDB2 = Csd;
    elseif group ==4
        Spec.DBDB4Ctx = SpecCtx;
        Spec.DBDB4Pyr = SpecPyr;
        Spec.DBDB4SLM = SpecSLM;
        CSD.DBDB4 = Csd;
    end % if
end % group
time = time -0.5; % adjust for SWR event initiation
%% save processed data
cd('C:\Users\ipzach\Documents\dbdb electrophy\Diabetes-Saved-Files')
save('LFP Measures','Spec','rip','labels','CSD')


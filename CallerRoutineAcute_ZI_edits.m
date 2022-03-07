%% This file loads the LFPs and calculates the mean modulation index
clear all

filepath = [filesep 'Volumes' filesep 'Ephy Mac' filesep 'Converted']; 
load([filesep 'Volumes' filesep 'Ephy Mac' filesep 'Converted' filesep 'SpkInfo.mat'])
% filepath = ['C:' filesep 'Users' filesep 'YazdanLab' filesep 'Documents' filesep 'MATLAB' filesep 'MI' filesep 'RealData'];
% load(['C:' filesep 'Users' filesep 'YazdanLab' filesep 'Documents' filesep 'MATLAB' filesep 'MI' filesep 'SpkInfo.mat'])
cd(filepath)

Fs = 1250; % Sampling Frequency; needed for filtering and plotting

%% 1.Load the right and left side signals (This part is similar to ThetaFreq.m)

% INPUT the group number and animal number
group = 14; % 1 for Acute dMCAO, 2 for Acute CCAO (See Spk_info.m)
animal = [4 5 7];%dMCAO - 3:12; % CCAO - 5 6 8 9 11 12;
numLFP = 2;
iside = 1;

% Layers
cortex = {'Ctx', 1};
CA1 = {'CA1', 2};
SLM = {'SLM', 3};

% INPUT the layer pairings you want for this file
phaseLayer = SLM; % SLM or CA1
ampLayer = CA1; % Either cortex OR SLM/CA1

%**INPUT the range of frequencies for which you want to get the mean MI**
phaseThetaLow = 4; % Lower limit frequency (Hz) for theta
phaseThetaHigh = 7; % Upper limit frequency (Hz) for theta

phaseDeltaLow = 0; % Lower limit frequency (Hz) for delta
phaseDeltaHigh = 3; % Upper limit frequency (Hz) for delta

ampHiGLow = 60; % Lower limit frequency (Hz) for high gamma
ampHiGHigh = 150; % Upper limit frequency (Hz) for high gamma

ampLoGLow = 30; % Lower limit frequency (Hz) for low gamma
ampLoGHigh = 60; % Upper limit frequency (Hz) for low gamma

% % INPUT the band to label the .mat files
% phaseBand = 'Theta';
% ampBand = 'HighGamma';

% If the group is 2, check to see if LFP2 is going to be calculated, and
% add another LFP index
CCAOhasLFP2 = false;
oldLFPvec = numLFP; % Have a place to store the original LFP values to calculate
if (group == 2)
    for numLFPidx = numLFP
        if (numLFPidx == 2)
            CCAOhasLFP2 = true;
            newNumLFP = zeros(numel(numLFP) + 1);
            newNumLFP = [numLFP max(numLFP)+1];
            numLFP = newNumLFP;
        end
    end
end
    
for animalIdx = animal
    for numLFPidx = numLFP
        pLayer{animalIdx,numLFPidx}.R = [];
        pLayer{animalIdx,numLFPidx}.L = [];
        aLayer{animalIdx,numLFPidx}.R = [];
        aLayer{animalIdx,numLFPidx}.L = [];
    end
end

for animalIdx = animal
    
    load([filepath filesep 'AllREMFiles' filesep SpkInfo{group,1} '_' num2str(animalIdx) filesep 'REM.mat' ]);
    disp([SpkInfo{group,1} '_' num2str(animalIdx) ])

    % Store the relevant LFP data for left and right side, cortical and
    % hippocampal LFPs
    for numLFPidx = numLFP

        if (numLFPidx > max(oldLFPvec))
            load([filepath filesep SpkInfo{group,1} '_' num2str(animalIdx) filesep 'LFP_' num2str(2)]);
        else 
            load([filepath filesep SpkInfo{group,1} '_' num2str(animalIdx) filesep 'LFP_' num2str(numLFPidx)]);
        end

        for iside = 1:2
            
            % Determine which LFPcellIdx to load: LFPs{1} or LFPs{2}.
            % For dMCAO: 3-9: LFPs{1} is left, {2} is right;
            %            10-12: LFPs{1} for all
            % For CCAO: LFPs{1} is left, {2} is right
            if (group == 1) && (animalIdx>2 && animalIdx<10)
                if iside == 1
                    LFPcellIdx = 2;
                elseif iside == 2
                    LFPcellIdx = 1;
                end
            elseif (group == 1) && (animalIdx<3 || animalIdx>9)
                LFPcellIdx = 1;
            elseif group == 2
                if iside == 1
                    LFPcellIdx = 2;
                elseif iside == 2
                    LFPcellIdx = 1;
                end
            elseif group == 13
                if iside == 1
                    LFPcellIdx = 2;
                elseif iside == 2
                    LFPcellIdx = 1;    
                end
            elseif group == 14
                if iside == 1
                    LFPcellIdx = 2;
                elseif iside == 2
                    LFPcellIdx = 1;    
                end
            end
            
            % load HTD periods
            if iside == 1

                % Right side
                disp(['Loading right side LFP ' num2str(numLFPidx) ' of animal ' num2str(animalIdx)])
                if ~isempty(SpkInfo{group,2}(animalIdx).R_chn{3})
                    HTDIdx = []; % Array for relevant data points indices
                    if (group == 2 && numLFPidx == 2)

                        % Get all of the indices from start to finish
                        for iRem = 1:length(rem(numLFPidx).R.start)
                            HTDIdx = [HTDIdx round(rem(numLFPidx).R.start(iRem)*Fs):round(rem(numLFPidx).R.end(iRem)*Fs)];    
                        end

                        % Only get the indices corresponding to the first 30
                        % minutes of recording
                        HTDIdx = HTDIdx(find(HTDIdx<=(30*60*Fs)));

                    elseif (group == 2 && numLFPidx == max(oldLFPvec)+1)

                        % Get all of the indices from start to finish
                        for iRem = 1:length(rem(2).R.start)
                            HTDIdx = [HTDIdx round(rem(2).R.start(iRem)*Fs):round(rem(2).R.end(iRem)*Fs)];
                        end

                        % Only get the indices corresponding to the second 30
                        % minutes of recording
                        HTDIdx = HTDIdx(find(HTDIdx>(30*60*Fs) & HTDIdx<=(60*60*Fs)));

                    else
                        for iRem = 1:length(rem(numLFPidx).R.start)
                            HTDIdx = [HTDIdx round(rem(numLFPidx).R.start(iRem)*Fs):round(rem(numLFPidx).R.end(iRem)*Fs)];
                        end

                        % If LFP 1 or 3 (any group) get first 30 minutes, if LFP2 group 1
                        % get last 30 minutes
                        if (numLFPidx == 1 || numLFPidx == 3)
                            HTDIdx = HTDIdx(find(HTDIdx<=(30*60*Fs)));
                        elseif (group == 1 && numLFPidx == 2)
                            HTDIdx = HTDIdx(find(HTDIdx>(30*60*Fs) & HTDIdx<=(60*60*Fs)));
                        end

                    end

                    % Create arrays of the relevant LFPs
                    if (numLFPidx == max(oldLFPvec)+1)
                        % First get the phase layer LFP
                        pLayerLFP = LFPs{LFPcellIdx}(:, SpkInfo{group,2}(animalIdx).R_chn{phaseLayer{1,2}}(2));
                        %CA1 = LFPs{2}(:, SpkInfo{group,2}(animal).R_chn{3}(2));
                        pLayer{animalIdx,numLFPidx}.R = pLayerLFP(HTDIdx)';

                        % Now get the amplitude layer LFP
                        aLayerLFP = LFPs{LFPcellIdx}(:, SpkInfo{group,2}(animalIdx).R_chn{ampLayer{1,2}}(2));
                        %ctx = LFPs{2}(:, SpkInfo{group,2}(animal).R_chn{1}(2));
                        aLayer{animalIdx,numLFPidx}.R = aLayerLFP(HTDIdx)';
                    else
                        % First get the phase layer LFP
                        pLayerLFP = LFPs{LFPcellIdx}(:, SpkInfo{group,2}(animalIdx).R_chn{phaseLayer{1,2}}(numLFPidx));
                        %CA1 = LFPs{2}(:, SpkInfo{group,2}(animal).R_chn{3}(numLFPidx));
                        pLayer{animalIdx,numLFPidx}.R = pLayerLFP(HTDIdx)';

                        % Now get the amplitude layer LFP
                        aLayerLFP = LFPs{LFPcellIdx}(:, SpkInfo{group,2}(animalIdx).R_chn{ampLayer{1,2}}(numLFPidx));
                        %ctx = LFPs{2}(:, SpkInfo{group,2}(animal).R_chn{1}(numLFPidx));
                        aLayer{animalIdx,numLFPidx}.R = aLayerLFP(HTDIdx)';
                    end
                end

            elseif iside == 2

                % Left side
                disp(['Loading left side LFP ' num2str(numLFPidx) ' of animal ' num2str(animalIdx)])
                if ~isempty(SpkInfo{group,2}(animalIdx).L_chn{3})
                    HTDIdx = []; % Array for relevant data points indices
                    if (group == 2 && numLFPidx == 2)

                        % Get all of the indices from start to finish
                        for iRem = 1:length(rem(numLFPidx).L.start)
                            HTDIdx = [HTDIdx round(rem(numLFPidx).L.start(iRem)*Fs):round(rem(numLFPidx).L.end(iRem)*Fs)];
                        end

                        % Only get the indices corresponding to the first 30
                        % minutes of recording
                        HTDIdx = HTDIdx(find(HTDIdx<=(30*60*Fs)));

                    elseif (group == 2 && numLFPidx == max(oldLFPvec)+1)

                        % Get all of the indices from start to finish
                        for iRem = 1:length(rem(2).L.start)
                            HTDIdx = [HTDIdx round(rem(2).L.start(iRem)*Fs):round(rem(2).L.end(iRem)*Fs)];
                        end

                        % Only get the indices corresponding to the second 30
                        % minutes of recording                    
                        HTDIdx = HTDIdx(find(HTDIdx>(30*60*Fs) & HTDIdx<=(60*60*Fs)));

                    else
                        for iRem = 1:length(rem(numLFPidx).L.start)
                            HTDIdx = [HTDIdx round(rem(numLFPidx).L.start(iRem)*Fs):round(rem(numLFPidx).L.end(iRem)*Fs)];
                        end

                        % If LFP 1 or 3 (any group) get first 30 minutes, if LFP2 group 1
                        % get last 30 minutes
                        if (numLFPidx == 1 || numLFPidx == 3)
                            HTDIdx = HTDIdx(find(HTDIdx<=(30*60*Fs)));
                        elseif (group == 1 && numLFPidx == 2)
                            HTDIdx = HTDIdx(find(HTDIdx>(30*60*Fs) & HTDIdx<=(60*60*Fs)));
                        end

                    end

                    % Creat an array of the relevant data points
                    if (numLFPidx == max(oldLFPvec)+1)
                        pLayerLFP = LFPs{LFPcellIdx}(:, SpkInfo{group,2}(animalIdx).L_chn{phaseLayer{1,2}}(2));
                        pLayer{animalIdx,numLFPidx}.L = pLayerLFP(HTDIdx)';

                        aLayerLFP = LFPs{LFPcellIdx}(:, SpkInfo{group,2}(animalIdx).L_chn{ampLayer{1,2}}(2));
                        aLayer{animalIdx,numLFPidx}.L = aLayerLFP(HTDIdx)';
                    else
                        pLayerLFP = LFPs{LFPcellIdx}(:, SpkInfo{group,2}(animalIdx).L_chn{phaseLayer{1,2}}(numLFPidx));
                        pLayer{animalIdx,numLFPidx}.L = pLayerLFP(HTDIdx)';

                        aLayerLFP = LFPs{LFPcellIdx}(:, SpkInfo{group,2}(animalIdx).L_chn{ampLayer{1,2}}(numLFPidx));
                        aLayer{animalIdx,numLFPidx}.L = aLayerLFP(HTDIdx)';
                    end

                end
            end
        end
    end
end
save([filepath filesep SpkInfo{group,1} ' ModulationIndex' phaseLayer{1} '_' ampLayer{1}], 'pLayer', 'aLayer');


%% 2.For Simulated LFP (Right and Left)               
% nonmodulatedamplitude=2; % increase this to get less modulation; you'll see that this is reflected in the MI value
% data_length = 3538900; % Length of the LFP data vectors <- this actually varies for each LFP
% t = 1:1:data_length;
% 
% Phase_Modulating_Freq=10;
% Amp_Modulated_Freq=80;
% 
% cort = cell(1,3);
% hippo = cell(1,3);
% 
% for numLFPidx = numLFP
%     
%     for iside = 1:2
%         if iside == 1
%             
%             % Simulated Right Side
%             lfp=(0.2*(sin(2*pi*t*Phase_Modulating_Freq/Fs)+1)+nonmodulatedamplitude*0.1).*sin(2*pi*t*Amp_Modulated_Freq/Fs)+sin(2*pi*t*Phase_Modulating_Freq/Fs);
%             cort{numLFP}.R = lfp+1*randn(1,length(lfp));
%             hippo{numLFPidx}.R = lfp+1*randn(1,length(lfp));
%             
%             
%         elseif iside == 2
%             % Left simulated LFPs
%             lfp=(0.2*(sin(2*pi*t*Phase_Modulating_Freq/Fs)+1)+nonmodulatedamplitude*0.1).*sin(2*pi*t*Amp_Modulated_Freq/Fs)+sin(2*pi*t*Phase_Modulating_Freq/Fs);
%             cort{numLFPidx}.L=lfp+1*randn(1,length(lfp));
%             hippo{numLFPidx}.L = lfp+1*randn(1,length(lfp));
%         end
%     end
% end
% % save('lfp', 'lfp');
% % load('lfp');

%% 3.Plot the LFPs

% for numLFPidx = numLFP
%     
%     figure
%     if ~isempty(pLayer{numLFPidx}.R)
%         subplot(2,2,1)
%         plot((0:length(pLayer{numLFPidx}.R)-1)/Fs, pLayer{numLFPidx}.R)
%         xlim([0 50])
%         set(gca,'fontsize',14)
%         xlabel('time (ms)')
%         ylabel('mV')
%         title(['Right Side SLM LFP ' num2str(numLFPidx)])
%     end
% 
%     if ~isempty(pLayer{numLFPidx}.L)
%         subplot(2,2,2)
%         plot((0:length(pLayer{numLFPidx}.L)-1)/Fs, pLayer{numLFPidx}.L)
%         xlim([0 50])
%         set(gca,'fontsize',14)
%         xlabel('time (ms)')
%         ylabel('mV')
%         title(['Left Side SLM LFP ' num2str(numLFPidx)])
%     end
% 
%     if ~isempty(aLayer{numLFPidx}.R)
%         subplot(2,2,3)
%         plot((0:length(aLayer{numLFPidx}.R)-1)/Fs, aLayer{numLFPidx}.R)
%         xlim([0 50])
%         set(gca,'fontsize',14)
%         xlabel('time (ms)')
%         ylabel('mV')
%         title(['Right Side Cortical LFP ' num2str(numLFPidx)])
%     end
% 
%     if ~isempty(aLayer{numLFPidx}.L)
%         subplot(2,2,4)
%         plot((0:length(aLayer{numLFPidx}.L)-1)/Fs, aLayer{numLFPidx}.L)
%         xlim([0 50])
%         set(gca,'fontsize',14)
%         xlabel('time (ms)')
%         ylabel('mV')
%         title(['Left Side Cortical LFP ' num2str(numLFPidx)])  
%     end
% end
%% 4.Define the Amplitude- and Phase- Frequencies

PhaseFreqVector = 0:0.5:10;% PhaseFreqVector = 0:2:10;
AmpFreqVector = 0:2.5:140; % AmpFreqVector = 0:5:150;
PhaseFreq_BandWidth = 0.5; % PhaseFreq_BandWidth = 3;
AmpFreq_BandWidth = 2.5; % AmpFreq_BandWidth = 10;

%% 5.For comodulation calculation (only has to be calculated once)
nbin = 18;
position = zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for jj = 1:nbin 
    position(jj) = -pi+(jj-1)*winsize; 
end
%% 6.Do filtering and Hilbert transform on CPU

disp('CPU filtering')
tic

Comodulogram = cell(length(animal),length(numLFP));
AmpFreqTransformed = cell(length(animal),length(numLFP));
PhaseFreqTransformed = cell(length(animal),length(numLFP));

for animalIdx = animal
    for numLFPidx = numLFP
        Comodulogram{animalIdx,numLFPidx}.R = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
        Comodulogram{animalIdx,numLFPidx}.L = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
        AmpFreqTransformed{animalIdx,numLFPidx}.R = zeros(length(AmpFreqVector), length(pLayer{animalIdx,numLFPidx}.R));
        AmpFreqTransformed{animalIdx,numLFPidx}.L = zeros(length(AmpFreqVector), length(pLayer{animalIdx,numLFPidx}.L));
        PhaseFreqTransformed{animalIdx,numLFPidx}.R = zeros(length(PhaseFreqVector), length(aLayer{animalIdx,numLFPidx}.R));
        PhaseFreqTransformed{animalIdx,numLFPidx}.L = zeros(length(PhaseFreqVector), length(aLayer{animalIdx,numLFPidx}.L));

        for iside = 1:2
            if (iside == 1 && (~isempty(pLayer{animalIdx,numLFPidx}.R) || ~isempty(aLayer{animalIdx,numLFPidx}.R)))
                % Right side
                toc
                disp(['Filtering right side LFP ' num2str(numLFPidx) ' of animal ' num2str(animalIdx)])

                for ii = 1:length(AmpFreqVector)
                    Af1 = AmpFreqVector(ii);
                    Af2 = Af1 + AmpFreq_BandWidth;
                    AmpFreq = eegfilt(pLayer{animalIdx,numLFPidx}.R, Fs, Af1, Af2); % just filtering
                    AmpFreqTransformed{animalIdx,numLFPidx}.R(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
                end

                for jj = 1:length(PhaseFreqVector)
                    Pf1 = PhaseFreqVector(jj);
                    Pf2 = Pf1 + PhaseFreq_BandWidth;
                    PhaseFreq = eegfilt(aLayer{animalIdx,numLFPidx}.R, Fs, Pf1, Pf2); % this is just filtering 
                    PhaseFreqTransformed{animalIdx,numLFPidx}.R(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
                end

            elseif (iside == 2 && (~isempty(pLayer{animalIdx,numLFPidx}.L) || ~isempty(aLayer{animalIdx,numLFPidx}.L)))
                % Left side
                toc
                disp(['Filtering left side LFP ' num2str(numLFPidx) ' of animal ' num2str(animalIdx)])

                for ii = 1:length(AmpFreqVector)
                    Af1 = AmpFreqVector(ii);
                    Af2 = Af1 + AmpFreq_BandWidth;
                    AmpFreq = eegfilt(pLayer{animalIdx,numLFPidx}.L, Fs, Af1, Af2); % just filtering
                    AmpFreqTransformed{animalIdx,numLFPidx}.L(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
                end

                for jj = 1:length(PhaseFreqVector)
                    Pf1 = PhaseFreqVector(jj);
                    Pf2 = Pf1 + PhaseFreq_BandWidth;
                    PhaseFreq = eegfilt(aLayer{animalIdx,numLFPidx}.L, Fs, Pf1, Pf2); % this is just filtering 
                    PhaseFreqTransformed{animalIdx,numLFPidx}.L(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
                end
            end
        end
    end
%     save([filepath filesep SpkInfo{group,1} ' ModulationIndex' phaseLayer{1} '_' ampLayer{1}], 'pLayer', 'aLayer');
    toc
end
disp('Filtering Completed')

%% 7.Do comodulation calculation
disp('Comodulation loop')

for animalIdx = animal
    for numLFPidx = numLFP
        disp(['Comodulation of LFP ' num2str(numLFPidx) ' of animal ' num2str(animalIdx)])
        for iside = 1:2
            if iside == 1 && (~isempty(PhaseFreqTransformed{animalIdx,numLFPidx}.R) || ~isempty(AmpFreqTransformed{animalIdx,numLFPidx}.R))

                toc
                % Right Side
                disp('Right Side LFP')

                counter1 = 0;
                for ii = 1:length(PhaseFreqVector)
                    counter1 = counter1+1;

                    Pf1 = PhaseFreqVector(ii);
                    Pf2 = Pf1+PhaseFreq_BandWidth;

                    counter2=0;
                    for jj = 1:length(AmpFreqVector)
                    counter2 = counter2+1;

                        Af1 = AmpFreqVector(jj);
                        Af2 = Af1+AmpFreq_BandWidth;
                        [MI,MeanAmp] = ModIndex_v2(PhaseFreqTransformed{animalIdx,numLFPidx}.R(ii, :), AmpFreqTransformed{animalIdx,numLFPidx}.R(jj, :), position);
                        Comodulogram{animalIdx,numLFPidx}.R(counter1,counter2) = MI;
                    end
                end
            elseif (iside == 2 && (~isempty(PhaseFreqTransformed{animalIdx,numLFPidx}.L) || ~isempty(AmpFreqTransformed{animalIdx,numLFPidx}.L)))

                toc
                % Left Side
                disp('Left Side LFP')

                counter1 = 0;
                for ii = 1:length(PhaseFreqVector)
                    counter1 = counter1+1;

                    Pf1 = PhaseFreqVector(ii);
                    Pf2 = Pf1+PhaseFreq_BandWidth;

                    counter2=0;
                    for jj = 1:length(AmpFreqVector)
                    counter2 = counter2+1;

                        Af1 = AmpFreqVector(jj);
                        Af2 = Af1+AmpFreq_BandWidth;
                        [MI,MeanAmp] = ModIndex_v2(PhaseFreqTransformed{animalIdx,numLFPidx}.L(ii, :), AmpFreqTransformed{animalIdx,numLFPidx}.L(jj, :), position);
                        Comodulogram{animalIdx,numLFPidx}.L(counter1,counter2) = MI;
                    end
                end
            end
        end
    end
end
save([filepath filesep SpkInfo{group,1} ' ModulationIndex' phaseLayer{1} '_' ampLayer{1}], 'pLayer', 'aLayer', 'Comodulogram');
toc

disp('Comodulation Loops Completed')

%% 8.Graph comodulograms

disp('Graphing Comodulograms');

% Decide where in the bin you want to show MI(beginning, middle, or end)
phaseVec = PhaseFreqVector+PhaseFreq_BandWidth/2;
comodPhaseFreqVec = PhaseFreqVector; %[0.5 phaseVec(1:end-1)];

ampVec = AmpFreqVector+AmpFreq_BandWidth/2;
comodAmpFreqVec = AmpFreqVector; %[0 ampVec(1:end-1)];

% Need to determine z axis limits for uniformity
maxVal = cell(max(animal), length(numLFP));
for animalIdx = animal
    for numLFPidx = numLFP
        maxVal{animalIdx,numLFPidx} = max([maxVal{animalIdx,numLFPidx} max(Comodulogram{animalIdx,numLFPidx}.R) max(Comodulogram{animalIdx,numLFPidx}.L)]);
    end
end

disp('Calculated Z Axis Limits');

if CCAOhasLFP2
    LFPTime = {'Baseline', 'Occlusion 0-30min', 'Reperfusion', 'Occlusion 30-60min'};
else
    LFPTime = {'Baseline', 'Occlusion 0-30min', 'Reperfusion', 'Occlusion 30-60min'};
end

% Determine the positions of each subplot
sLotPos.L = {[0.085 0.5 0.1725 0.39], [0.2925 0.5 0.1725 0.39], [0.5 0.5 0.1725 0.39], [0.7075 0.5 0.2 0.39]};
sLotPos.R = {[0.085 0.08 0.1725 0.39], [0.2925 0.08 0.1725 0.39], [0.5 0.08 0.1725 0.39], [0.7075 0.08 0.2 0.39]};

for animalIdx = animal
    
    disp([SpkInfo{group,1} ' ' num2str(animalIdx) 'Comodulgrams'])
    
    figure('Name', [SpkInfo{group,1} ' ' num2str(animalIdx) ' Right Side'], 'units', 'normalized', 'outerposition', [0 0 1 1])
    for numLFPidx = numLFP

        if group == 1
            switch numLFPidx
                case 1
                    sLotIdx = 1;
                case 2
                    sLotIdx = 3;
                case 3
                    sLotIdx = 4;
            end
        elseif group == 2
            switch numLFPidx
                case 1
                    sLotIdx = 1;
                case 2
                    sLotIdx = 2;
                case 3
                    sLotIdx = 4;
                case 4
                    sLotIdx = 3;
            end
        elseif group == 13
            switch numLFPidx
                case 1
                    sLotIdx = 1;
                case 2
                    sLotIdx = 2;
            end
         elseif group == 14
            switch numLFPidx
                case 1
                    sLotIdx = 1;
                case 2
                    sLotIdx = 2;
            end
        end
        disp(['Right LFP ' num2str(numLFPidx)])
                
        if ~any(isnan(Comodulogram{animalIdx,numLFPidx}.R'))
            
            subplot('Position', sLotPos.R{1,sLotIdx})
%             subplot(2,length(numLFP),numLFPidx+3)
            % contourf(comodPhaseFreqVec,comodAmpFreqVec,Comodulogram{animalIdx,numLFPidx}.R',30,'lines','none')
            pcolor(comodPhaseFreqVec,comodAmpFreqVec,Comodulogram{animalIdx,numLFPidx}.R')
            shading interp
            colormap(cubehelix(800,2.5,0.8,1.4,0.55))
            set(gca,'fontsize',14);
            if sLotIdx == 1
                hy = ylabel('Contralateral');
                set(hy, 'fontsize', 16);
            else
                set(gca, 'ytick', []);
            end
            if sLotIdx == 4
                colorbar
            end
            caxis([0 maxVal{animalIdx,numLFPidx}]);
            rectangle('Position', [phaseThetaLow ampHiGLow phaseThetaHigh-phaseThetaLow ampHiGHigh-ampHiGLow], 'Edgecolor', 'r', 'Linewidth', 1.5)
            rectangle('Position', [phaseThetaLow ampLoGLow phaseThetaHigh-phaseThetaLow ampLoGHigh-ampLoGLow], 'Edgecolor', 'g', 'Linewidth', 1.5)
            rectangle('Position', [phaseDeltaLow ampHiGLow phaseDeltaHigh-phaseDeltaLow ampHiGHigh-ampHiGLow], 'Edgecolor', 'm', 'Linewidth', 1.5)
            rectangle('Position', [phaseDeltaLow ampLoGLow phaseDeltaHigh-phaseDeltaLow ampLoGHigh-ampLoGLow], 'Edgecolor', 'w', 'Linewidth', 1.5)
%             ylabel([ampLayer{1} ' Amplitude Frequency (Hz)'])
%             xlabel([phaseLayer{1} ' Phase Frequency (Hz)'])
%             colorbar
%             title([LFPTime{numLFPidx} ' ' phaseLayer{1} '-' ampLayer{1} ' Right Comodulogram'])
%             saveas(gcf, [SpkInfo{group,1} '_' num2str(animal) '_RComod' LFPTime{numLFPidx} '_' phaseLayer{1} '-' ampLayer{1} '.fig'])
        else
            disp(['Unable to generate right Comodulogram for LFP ' num2str(numLFPidx)])
        end

        disp(['Left LFP ' num2str(numLFPidx)])
        
        if ~any(isnan(Comodulogram{animalIdx,numLFPidx}.L'))
            subplot('Position', sLotPos.L{1,sLotIdx})
%             subplot(2,3,numLFPidx)
            contourf(comodPhaseFreqVec,comodAmpFreqVec,Comodulogram{animalIdx,numLFPidx}.L',30,'lines','none')
            set(gca,'fontsize',14)
            if sLotIdx == 1
                hy = ylabel('Ipsilateral');
                set(hy, 'fontsize', 16);
            else
                set(gca, 'ytick', []);
            end
            if sLotIdx == 4
                colorbar
            end
            caxis([0 maxVal{animalIdx,numLFPidx}]);
            rectangle('Position', [phaseThetaLow ampHiGLow phaseThetaHigh-phaseThetaLow ampHiGHigh-ampHiGLow], 'Edgecolor', 'r', 'Linewidth', 1.5)
            rectangle('Position', [phaseThetaLow ampLoGLow phaseThetaHigh-phaseThetaLow ampLoGHigh-ampLoGLow], 'Edgecolor', 'g', 'Linewidth', 1.5)
            rectangle('Position', [phaseDeltaLow ampHiGLow phaseDeltaHigh-phaseDeltaLow ampHiGHigh-ampHiGLow], 'Edgecolor', 'm', 'Linewidth', 1.5)
            rectangle('Position', [phaseDeltaLow ampLoGLow phaseDeltaHigh-phaseDeltaLow ampLoGHigh-ampLoGLow], 'Edgecolor', 'w', 'Linewidth', 1.5)
            if group == 1 && numLFPidx == 2
                title(LFPTime{4}, 'fontweight', 'normal');
            else
                title(LFPTime{numLFPidx}, 'fontweight', 'normal');
            end
        else
            disp(['Unable to generate left Comodulogram for LFP ' num2str(numLFPidx)])
        end
    end
    
    % Set the suplabels
    [tA, htA] = suplabel([SpkInfo{group,1} ' ' num2str(animalIdx)], 't');
    [xA, hxA] = suplabel([phaseLayer{1,1} ' Phase Frequency (Hz)'], 'x');
    [yA, hyA] = suplabel([ampLayer{1,1} ' Amplitude Frequency (Hz)'], 'y');
    set(htA, 'fontsize', 20, 'Position', [0.5 1.03 0.5]);
    set(hxA, 'fontsize', 18, 'Position', [0.5 0.0 0]);
    set(hyA, 'fontsize', 18, 'Position', [-0.02 0.50 0]);
    
    % Save the figure
    saveas(gcf, [filepath filesep SpkInfo{group,1} '_' num2str(animalIdx) filesep phaseLayer{1,1} '_' ampLayer{1,1}], 'fig');
    saveas(gcf, [filepath filesep SpkInfo{group,1} '_' num2str(animalIdx) filesep phaseLayer{1,1} '_' ampLayer{1,1}], 'tif');
    
end

disp('Finished Graphing');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9.Use the routine below to look at specific pairs of frequency range:

% Pf1 = 6
% Pf2 = 12
% Af1 = 60
% Af2 = 100
% 
% [MIR,MeanAmpR] = ModIndex_v1(sides.R,srate,Pf1,Pf2,Af1,Af2)
% 
% [MIL,MeanAmpL] = ModIndex_v1(sides.L,srate,Pf1,Pf2,Af1,Af2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10.Or use the routine below to make a comodulogram using ModIndex_v1; this takes longer than
%% the method outlined above using ModIndex_v2 because in this routine multiple filtering of the same
%% frequency range is employed (the Amp frequencies are filtered multiple times, one
%% for each phase frequency). This routine might be the only choice though
%% for computers with low memory, because it does not create the matrices
%% AmpFreqTransformed and PhaseFreqTransformed as the routine above
% 
% tic
% 
% PhaseFreqVector = 2:2:50;
% AmpFreqVector = 10:5:200;
% 
% PhaseFreq_BandWidth = 4;
% AmpFreq_BandWidth = 10;
% 
% ComodulogramR = zeros(length(PhaseFreqVector),length(AmpFreqVector));
% ComodulogramL = zeros(length(PhaseFreqVector),length(AmpFreqVector));
% 
% for iside = 1:2
%     if iside == 1
%         % Right Side
%         counter1 = 0;
%         for Pf1 = PhaseFreqVector
%             counter1 = counter1+1;
%             %Pf1 % just to check the progress
%             Pf2 = Pf1+PhaseFreq_BandWidth;
% 
%             counter2 = 0;
%             for Af1 = AmpFreqVector
%                 counter2 = counter2+1;
%                 Af2 = Af1+AmpFreq_BandWidth;
% 
%                 [MI,MeanAmp] = ModIndex_v1(sides.R,srate,Pf1,Pf2,Af1,Af2);
% 
%                 ComodulogramR(counter1,counter2) = MI;
%             end 
%         end
%     else
%         toc
%         % Left Side
%         counter1 = 0;
%         for Pf1 = PhaseFreqVector
%             counter1 = counter1+1;
%             %Pf1 % just to check the progress
%             Pf2 = Pf1+PhaseFreq_BandWidth;
% 
%             counter2 = 0;
%             for Af1 = AmpFreqVector
%                 counter2 = counter2+1;
%                 Af2 = Af1+AmpFreq_BandWidth;
% 
%                 [MI,MeanAmp] = ModIndex_v1(sides.L,srate,Pf1,Pf2,Af1,Af2);
% 
%                 ComodulogramL(counter1,counter2) = MI;
%             end 
%         end
%     end
% end
% toc

%% 11.Graph Comodulograms
% 
% clf
% contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,ComodulogramR',30,'lines','none')
% set(gca,'fontsize',14)
% ylabel('Amplitude Frequency (Hz)')
% xlabel('Phase Frequency (Hz)')
% colorbar
% title('Right Comodulogram')
% 
% contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,ComodulogramL',30,'lines','none')
% set(gca,'fontsize',14)
% ylabel('Amplitude Frequency (Hz)')
% xlabel('Phase Frequency (Hz)')
% colorbar
% title('Left Comodulogram')

%% 12.Get the average MI

PFreqRangeTheta = phaseThetaLow:2:phaseThetaHigh; % phase frequency range (increments of 1 just like PhaseFreqVector)
AFreqRangeHiG = ampHiGLow:5:ampHiGHigh; % amp frequency range (increments of 5 just like AmpFreqVector)
AFreqRangeLoG = ampLoGLow:5:ampLoGHigh; % amp frequency range (increments of 5 just like AmpFreqVector)
PFreqRangeDelta = phaseDeltaLow:2:phaseDeltaHigh; % phase frequency range (increments of 1 just like PhaseFreqVector)

avgMIThetaHiG = cell(length(animal),length(numLFP));
avgMIThetaLoG = cell(length(animal),length(numLFP));
avgMIDeltaHiG = cell(length(animal),length(numLFP));
avgMIDeltaLoG = cell(length(animal),length(numLFP));

% Get the mean MI values within the specified frequency ranges
for animalIdx = animal
    for numLFPidx = numLFP
        avgMIThetaHiG{animalIdx,numLFPidx}.R = CalcAvgMI(comodPhaseFreqVec, comodAmpFreqVec, Comodulogram{animalIdx,numLFPidx}.R, PFreqRangeTheta, AFreqRangeHiG);
        avgMIThetaHiG{animalIdx,numLFPidx}.L = CalcAvgMI(comodPhaseFreqVec, comodAmpFreqVec, Comodulogram{animalIdx,numLFPidx}.L, PFreqRangeTheta, AFreqRangeHiG);
        
        avgMIThetaLoG{animalIdx,numLFPidx}.R = CalcAvgMI(comodPhaseFreqVec, comodAmpFreqVec, Comodulogram{animalIdx,numLFPidx}.R, PFreqRangeTheta, AFreqRangeLoG);
        avgMIThetaLoG{animalIdx,numLFPidx}.L = CalcAvgMI(comodPhaseFreqVec, comodAmpFreqVec, Comodulogram{animalIdx,numLFPidx}.L, PFreqRangeTheta, AFreqRangeLoG);
        
        avgMIDeltaHiG{animalIdx,numLFPidx}.R = CalcAvgMI(comodPhaseFreqVec, comodAmpFreqVec, Comodulogram{animalIdx,numLFPidx}.R, PFreqRangeDelta, AFreqRangeHiG);
        avgMIDeltaHiG{animalIdx,numLFPidx}.L = CalcAvgMI(comodPhaseFreqVec, comodAmpFreqVec, Comodulogram{animalIdx,numLFPidx}.L, PFreqRangeDelta, AFreqRangeHiG);
        
        avgMIDeltaLoG{animalIdx,numLFPidx}.R = CalcAvgMI(comodPhaseFreqVec, comodAmpFreqVec, Comodulogram{animalIdx,numLFPidx}.R, PFreqRangeDelta, AFreqRangeLoG);
        avgMIDeltaLoG{animalIdx,numLFPidx}.L = CalcAvgMI(comodPhaseFreqVec, comodAmpFreqVec, Comodulogram{animalIdx,numLFPidx}.L, PFreqRangeDelta, AFreqRangeLoG);
    end
end

disp('Calculated Average MI')
%% 13.Save to a file in the current directory

save([filepath filesep SpkInfo{group,1} ' ModulationIndex' phaseLayer{1} '_' ampLayer{1}], 'avgMIThetaHiG', 'avgMIThetaLoG', 'avgMIDeltaHiG', 'avgMIDeltaLoG', 'pLayer', 'aLayer', 'Comodulogram', 'phaseThetaLow', 'phaseThetaHigh', 'phaseDeltaLow', 'phaseDeltaHigh', 'ampHiGLow', 'ampHiGHigh', 'ampLoGLow', 'ampLoGHigh');
disp('Saved Everything')
disp('Good job!')
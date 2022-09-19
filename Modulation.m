% load in spkinf and concise channel info
% clear all, close all
load('SpkInfo.mat')
load('chans.mat')
% file path for
filepath = 'C:\Users\ipzach\Documents\dbdb electrophy';
cd(filepath)
animalList = dir;
Fs = 1250; % Sampling Frequency; needed for filtering and plotting

TD = 1; %High 1, Low 2, Full 3
phase = 0; %0 Cortical layer is phase, 1 is hippocampal layer
%%%%%%%%%%%%%%%%%%%%%%%%%
switch phase
    case 0
        disp('Cortex is phase variable')
    case 1
        disp('Hippocampus is phase variable')
end
% Define the Amplitude- and Phase- Frequencies
% Lower this range to have faster computing
PhaseFreqVector = 0:0.5:12 ;
AmpFreqVector = 0:2.5:160;
PhaseFreq_BandWidth = 0.5;
AmpFreq_BandWidth = 2.5;
% For comodulation calculation (only has to be calculated once)
nbin = 0:17; % defines how many fractions to divide phase into
winsize = 2*pi/18;
position = nbin*winsize-pi;


for group = 1:4
    % 1.Load the right and left side signals
    %load H/L TD indexes
    % set files to load skipping animals with bad channels
    if group ==1
        grouping = 3:9; % DB+ 200D
    elseif group ==2
        grouping = 10:14; % DB+ 400D
    elseif group ==3
        grouping = [15:18 20 21]; % DBDB 200D
    elseif group ==4
        grouping = [22 24:27]; % DBDB 400D
    end
    
    % run through the correct animals for each group
    for animal = grouping
        disp(['Animal: ' num2str(animal)])
        cd(animalList(animal).name)
        load('REM.mat');
        disp('loading data')
        files = dir('LFP*');
        for recording =  1:length(rem(1).R.start)
            load(files(recording).name);
            LFP = LFPs{1,1};
            for layer = 1:3
                switch layer
                    case 1
                        A = 1;
                        B = 2;
                        A_name = 'Cortex';
                        B_name =  'Pyr';
                    case 2
                        A = 1;
                        B = 3;
                        A_name = 'Cortex';
                        B_name = 'Slm';
                    case 3
                        A = 2;
                        B = 3;
                        A_name = 'Pyr';
                        B_name = 'Slm';
                end % switch layer
                switch phase
                    case 0
                        amp_sig  = LFP(:,chans(A,animal));
                        phase_sig = LFP(:,chans(B,animal));
                        name = [A_name '-' B_name];
                    case 1
                        amp_sig  = LFP(:,chans(B,animal));
                        phase_sig = LFP(:,chans(A,animal));
                        name = [B_name '-' A_name];
                end % switch phase
                
                TDIdx  = [];
                switch TD
                    case 1
                        % Get indices of High Theta/Delta
                        for iRem = 1:length(rem(recording).R.start)
                            TDIdx = [TDIdx (round(rem(recording).R.start(iRem)*Fs)):(round(rem(recording).R.end(iRem)*Fs))];
                        end
                        REMtitle = 'H';
                        
                        % Get indices of Low Theta/Delta
                    case 2
                        LowRightStart = rem(recording).R.end(1:end-1);
                        LowRightEnd = rem(recording).R.start(2:end);
                        for iRem = 1:length(LowRightStart)
                            TDIdx = [TDIdx (round(LowRightStart(iRem)*Fs)):(round(LowRightEnd(iRem)*Fs))];
                        end
                        REMtitle = 'L';
                        
                        % For entire reading
                    case 3
                        TDIdx = 1 : length(LFP);
                        REMtitle = 'F';
                end
                TD_amp_sig = amp_sig(TDIdx);
                TD_phase_sig = phase_sig(TDIdx);
                %try
                    disp('Filtering')
                    Comodulogram = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
                    AmpFreqTransformed = zeros(length(AmpFreqVector), length(TD_amp_sig));
                    PhaseFreqTransformed = zeros(length(PhaseFreqVector), length(TD_phase_sig));
                    for ii = 1:length(AmpFreqVector)
                        Af1 = AmpFreqVector(ii); % Amplitude frequency start
                        Af2 = Af1 + AmpFreq_BandWidth; % amplitude frequency end
                        AmpFreq = eegfilt(TD_amp_sig', Fs, Af1, Af2); % just filtering
                        AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
                    end
                    for jj = 1:length(PhaseFreqVector)
                        Pf1 = PhaseFreqVector(jj);
                        Pf2 = Pf1 + PhaseFreq_BandWidth;
                        PhaseFreq = eegfilt(TD_phase_sig', Fs, Pf1, Pf2); % this is just filtering
                        PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
                    end
%                 catch
%                     disp('Error during filtering')
%                     continue
%                 end % try
                disp('Calculating comodulation')
                counter1= 0;
                for ii = 1:length(PhaseFreqVector)
                    counter1 = counter1+1;
                    Pf1 = PhaseFreqVector(ii);
                    Pf2 = Pf1+PhaseFreq_BandWidth;
                    counter2=0;
                    for jj = 1:length(AmpFreqVector)
                        counter2 = counter2+1;
                        Af1 = AmpFreqVector(jj);
                        Af2 = Af1+AmpFreq_BandWidth;
                        [MI,MeanAmp] = ModIndex_v2(PhaseFreqTransformed(ii, :), AmpFreqTransformed(jj, :), position);
                        Comodulogram(counter1,counter2) = MI;
                    end % for jj
                end  % for  ii
                disp('Saving')
                save(['Comodulogram_' name '_' num2str(recording)], 'Comodulogram','TD_amp_sig','TD_phase_sig');
                
                disp('Plotting')
                figure
                contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulogram',30,'lines','none')
                title([num2str(animal) ' ' name ' ' num2str(recording)]), colormap('jet')
                colorbar
                drawnow
            end % for layer
        end % for recording
        cd ..
    end % for animal
end % for group




















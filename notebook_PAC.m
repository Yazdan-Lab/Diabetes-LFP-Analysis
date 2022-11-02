% load in spkinf and concise channel info
% clear all, close all
load('SpkInfo.mat')
load('chans.mat')

% addpath('C:\Users\ipzach\Documents\MATLAB\Toolbox Zach',...
%     'C:\Users\ipzach\Documents\MATLAB\spectral-analysis-tools')

 addpath('C:\COM\ePhy\dbdb\code\utils-toolbox\utils-toolbox')
  addpath('C:\COM\ePhy\dbdb\code\spectral-analysis-tools')

% file path for
%filepath = 'C:\Users\ipzach\Documents\MATLAB\Data\dbdb electrophy';
% savepath = 'C:\Users\ipzach\Documents\MATLAB\output\Diabetes-Saved-Files\';
filepath = 'C:\COM\ePhy\dbdb\Data\dbdb electrophy';
savepath = 'C:\COM\ePhy\dbdb\Data\Outputs\';

cd(filepath)
animalList = dir;
Fs = 1250; % Sampling Frequency; needed for filtering and plotting

TD = 1; %High 1, Low 2
switch TD
    case 1
        method = 'high';
    case 2
        method = 'low';
end
phase = false; % 0 Cortical layer is phase, 1 is hippocampal layer
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
        cd([filepath '\' animalList(animal).name])
        load('REM.mat');
        disp('loading data')
        files = dir('LFP*');
        for recording =  1:length(files)
            try
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
                        amp_signal  = LFP(:,chans(A,animal));
                        phase_signal = LFP(:,chans(B,animal));
                        name = [A_name '-' B_name];
                    case 1
                        amp_signal  = LFP(:,chans(B,animal));
                        phase_signal = LFP(:,chans(A,animal));
                        name = [B_name '-' A_name];
                end % switch phase
                
                REMtitle = 'H';
                TDIdx = isolate_state_idx(rem(recording).R.start, ...
                    rem(recording).R.end, ...
                    Fs,...
                    method);
                
                TD_amp_signal = amp_signal(TDIdx);
                TD_phase_signal = phase_signal(TDIdx);
                
                disp('Filtering')
                
                
                AmpFreqTransformed = abs(...
                    bin_filter_hilbert(TD_amp_signal,...
                        AmpFreqVector,...
                        AmpFreq_BandWidth,...
                        Fs)...
                    );
                PhaseFreqTransformed = angle(...
                        bin_filter_hilbert(TD_phase_signal,...
                        PhaseFreqVector,...
                        PhaseFreq_BandWidth,...
                        Fs)...
                    );
                
                disp('Calculating comodulation')
                Comodulogram = comodulate(AmpFreqTransformed, PhaseFreqTransformed);
                
                disp(['Saving: ' savepath 'Comodulogram_' name '_' num2str(recording)])
                save([savepath 'Comodulogram_' num2str(animal) '_' name '_' num2str(recording)], 'Comodulogram','TD_amp_signal','TD_phase_signal');
                
                disp('Plotting')
                figure
                contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulogram',30,'lines','none');
                title([num2str(animal) ' ' name ' ' num2str(recording)]), colormap('jet')
                colorbar
                drawnow
                %MS
                Datetime_PAC = string(datetime('now'));
                Filename_PAC = sprintf('PAC_Figure_%s.tiff', Datetime_PAC);
                Filename_PAC = regexprep(Filename_PAC, ' ', '_');
                Filename_PAC = regexprep(Filename_PAC, ':', '_');
                saveas(gcf, Filename_PAC);                
                %ME
            end % for layer
            catch
                disp('No recording exists, continuing')
                continue
            end % try catch
        end % for recording
        cd ..
    end % for animal
end % for group
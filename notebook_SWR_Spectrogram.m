%% SWR CSD Analysis
% This script will load LFP and grab the LFP 500ms before and 1500ms after
% the onset of a ripple (max duration recorded 1237ms) then perform CSD on
% that file, save the file to a 3D matrix, and create an average CSD
clear all; close all; clc
% initalize filepaths and load reference data files
%cd('C:\Users\ipzach\Documents\dbdb electrophy\General_Scripts')
cd('C:\COM\ePhy\dbdb\code\Diabetes-LFP-Analysis')
load('chans.mat') % This file lets us know where in the brain each electrode is
%cd('C:\Users\ipzach\Documents\dbdb electrophy'); % here is the data
cd('C:\COM\ePhy\dbdb\Data\dbdb electrophy')
animals = dir; % create a list of every folder (each one is one animal)
voltConv = 0.000000091555527603759401; % Neurolynx saves data in a unitless value, we need this to convert it to volts
for i = 1:4 % number of groups
    % Here we go through the folders of animals, defining which folders
    % belong to which animals
    if i ==1
        grouping = 3:9; % DB+ 200D
    elseif i ==2
        grouping = 10:14; % DB+ 400D
    elseif i ==3
        grouping = [15:18 20 21]; % DBDB 200D
    elseif i ==4
        grouping = [22 24:27]; % DBDB 400D
    end
    % for each animal re-initialize empty variables to hold data
    SwrFullSpecPyr = [];
    SwrFullSpecSLM = [];
    SwrFullSpecCtx = [];
    
    for j = grouping % Run through each animal in a group
        cd(animals(j).name) % Enter that animals folder
        load('SWR_Index.mat'); % load SWR HTD/LTD data
        SWR_files = dir('SWR_R_*'); % Grab the number of files that have SWR event timings (usually 2, sometimes 1 or 3)
        SWR_files = {SWR_files.name}; % throw away useless info
        
        LFP_files = dir('LFP*'); % grab number of LFP files (sometimes 1 or 3 as well)
        LFP_files = {LFP_files.name}; % throw away useless info
        
        for k = 1:size(SWRLTDIdx,2) % run through each of the LTD periods (where SWRs occur)
            if ~isempty(SWRLTDIdx(k).R) % makes sure ripple occured during this period
                load(char(SWR_files(k)));  % Load SWR events 
                load(char(LFP_files(k))); % Load LFP events 
                LFP = LFPs{1,2} .*voltConv; % load LFP
                LFP = BPfilter(LFP,1250,30,60); % isolate gamma frequency band
                LTDevents = SWRevents(SWRLTDIdx(k).R,1); % grab SWRs that occur during this period
                LTDevents = LTDevents(LTDevents >= 625); % Exclude SWRs that clip beginning window
                LTDevents = LTDevents(LTDevents <= length(LFP)-1250);% Exclude SWRs that exclude final window
                
                % initialize some temp variables
                tempSpecCtx = zeros(50,28,length(LTDevents));
                tempSpecPyr = zeros(50,28,length(LTDevents));
                tempSpecSLM = zeros(50,28,length(LTDevents));
                
                % For each event, create a spectrogram and store it in the
                % temp variable
                for l = 1:length(LTDevents)
                    [~,~,~,tempSpecCtx(:,:,l)] = spectrogram(LFP(LTDevents(l)-625:LTDevents(l)+1250,chans(1,j)),hamming(125),[],[5:5:250],1250);
                    [~,~,~,tempSpecPyr(:,:,l)] = spectrogram(LFP(LTDevents(l)-625:LTDevents(l)+1250,chans(2,j)),hamming(125),[],[5:5:250],1250);
                    [~,~,~,tempSpecSLM(:,:,l)] = spectrogram(LFP(LTDevents(l)-625:LTDevents(l)+1250,chans(3,j)),hamming(125),[],[5:5:250],1250);
                end
                
                % Concatenate temp variable to storage variable
                if isempty(SwrFullSpecPyr)
                    SwrFullSpecPyr = tempSpecPyr;
                    SwrFullSpecCtx = tempSpecCtx;
                    SwrFullSpecSLM = tempSpecSLM;
                else
                    SwrFullSpecCtx = cat(3,SwrFullSpecCtx,tempSpecCtx);
                    SwrFullSpecPyr = cat(3,SwrFullSpecPyr,tempSpecPyr);
                    SwrFullSpecSLM = cat(3,SwrFullSpecSLM,tempSpecSLM);
                    
                end % if isempty FullSPec
            end %is empty SWRLTD
            
        end % for k SWRLTDIdx
        cd .. % Jump back into the folder stack so we can enter the next animal
        
    end % grouping
    % store variable ouside loop
        if i ==1
            Spec.DB2Ctx = SwrFullSpecCtx;
            Spec.DB2Pyr = SwrFullSpecPyr;
            Spec.DB2SLM = SwrFullSpecSLM;
        elseif i ==2
            Spec.DB4Ctx = SwrFullSpecCtx;
            Spec.DB4Pyr = SwrFullSpecPyr;
            Spec.DB4SLM = SwrFullSpecSLM;
        elseif i ==3
            Spec.DBDB2Ctx = SwrFullSpecCtx;
            Spec.DBDB2Pyr = SwrFullSpecPyr;
            Spec.DBDB2SLM = SwrFullSpecSLM;
        elseif i ==4
            Spec.DBDB4Ctx = SwrFullSpecCtx;
            Spec.DBDB4Pyr = SwrFullSpecPyr;
            Spec.DBDB4SLM = SwrFullSpecSLM;
        end % if
end % treatment
% grab frequency range and time for plots
[~,freqs,time,~] = spectrogram(LFP(LTDevents(l)-625:LTDevents(l)+1250,chans(1,j)),hamming(125),[],[5:5:250],1250);

% Take the average across all SWRs
Spec.DB2CtxM = mean(Spec.DB2Ctx,3);
Spec.DB2PyrM = mean(Spec.DB2Pyr,3);
Spec.DB2SLMM = mean(Spec.DB2SLM,3);

Spec.DB4CtxM = mean(Spec.DB4Ctx,3);
Spec.DB4PyrM = mean(Spec.DB4Pyr,3);
Spec.DB4SLMM = mean(Spec.DB4SLM,3);

Spec.DBDB2CtxM = mean(Spec.DBDB2Ctx,3);
Spec.DBDB2PyrM = mean(Spec.DBDB2Pyr,3);
Spec.DBDB2SLMM = mean(Spec.DBDB2SLM,3);

Spec.DBDB4CtxM = mean(Spec.DBDB4Ctx,3);
Spec.DBDB4PyrM = mean(Spec.DBDB4Pyr,3);
Spec.DBDB4SLMM = mean(Spec.DBDB4SLM,3);

%% Plot
maxVal = 2e-10;
figure
subplot(3,4,1)
h1 = heatmap(time,flipud(freqs), flipud(Spec.DB2CtxM));
h1.GridVisible = 'off';
title('Control 200')
ylabel('Cortex')
caxis([0 maxVal])

subplot(3,4,2)
h2 = heatmap(time,flipud(freqs), flipud(Spec.DB4CtxM));
h2.GridVisible = 'off';
title('Control 400')
caxis([0 maxVal])

subplot(3,4,3)
h3 = heatmap(time,flipud(freqs), flipud(Spec.DBDB2CtxM));
h3.GridVisible = 'off';
title('DBDB 200')
caxis([0 maxVal])

subplot(3,4,4)
h4 = heatmap(time,flipud(freqs), flipud(Spec.DBDB4CtxM));
h4.GridVisible = 'off';
title('DBDB 400')
caxis([0 maxVal])

subplot(3,4,5)
h5 = heatmap(time,flipud(freqs), flipud(Spec.DB2PyrM));
h5.GridVisible = 'off';
ylabel('Pyramidal')
caxis([0 maxVal])

subplot(3,4,6)
h6 = heatmap(time,flipud(freqs), flipud(Spec.DB4PyrM));
h6.GridVisible = 'off';
caxis([0 maxVal])

subplot(3,4,7)
h7 = heatmap(time,flipud(freqs), flipud(Spec.DBDB2PyrM));
h7.GridVisible = 'off';
caxis([0 maxVal])

subplot(3,4,8)
h8 = heatmap(time,flipud(freqs), flipud(Spec.DBDB4PyrM));
h8.GridVisible = 'off';
caxis([0 maxVal])

subplot(3,4,9)
h9 = heatmap(time,flipud(freqs), flipud(Spec.DB2SLMM));
h9.GridVisible = 'off';
ylabel('SLM')
caxis([0 maxVal])

subplot(3,4,10)
h10 = heatmap(time,flipud(freqs), flipud(Spec.DB4SLMM));
h10.GridVisible = 'off';
caxis([0 maxVal])

subplot(3,4,11)
h11 = heatmap(time,flipud(freqs), flipud(Spec.DBDB2SLMM));
h11.GridVisible = 'off';
caxis([0 maxVal])

subplot(3,4,12)
h12 = heatmap(time,flipud(freqs), flipud(Spec.DBDB4SLMM));
h12.GridVisible = 'off';
caxis([0 maxVal])



ageLabs = [label.DB2age; label.DBDB2age; label.DB4age; label.DBDB4age];
treatLabs = [label.DB2treat; label.DBDB2treat; label.DB4treat; label.DBDB4treat];


% NOW we can run the stats!
[ctxP,ctxT,ctxStats] = anovan(ctxVals,{ treatLabs ageLabs},'model','interaction');
[ctxC,ctxM,~,ctxN] = multcompare(ctxStats,'Dimension',[1 2],'CType','bonferroni');

[pyrP,pyrT,pyrStats] = anovan(pyrVals,{ treatLabs ageLabs},'model','interaction');
[pyrC,pyrM,~,pyrN] = multcompare(pyrStats,'Dimension',[1 2],'CType','bonferroni');

[slmP,slmT,slmStats] = anovan(slmVals,{ treatLabs ageLabs},'model','interaction');
[slmC,slmM,~,slmN] = multcompare(slmStats,'Dimension',[1 2],'CType','bonferroni');

%% Final Figures
% Initialize figure settings
figure
set(gcf,'Color','w','Position',[100 100 1200 500])
x1 = 0.07;
x2 = 0.39;
x3 = 0.71;
y = 0.15;
h = 0.8;
w = 0.24;


subplot('Position',[x1 y w h ])
[ctxbar] = UCSF_graph([ctxM(1:2,2),ctxM(3:4,2)]',[ctxM(1:2,1),ctxM(3:4,1)]',ctxC);
ylabel('Gamma Power (v)','Fontsize',18)
legend('db/+','db/db');
legend('boxoff')
legend('Location','northoutside')
title('Cortex');
ctx_age_sig = sig_check(ctxP(2));
ctx_db_sig = sig_check(ctxP(1));

if ctxP(2) <= 0.05
A = suplabel(['age effect: ' ctx_age_sig],'t',[0.3 0.08 0.01 0.83]);
set(A,'FontSize',12)
end
if ctxP(1) <= 0.05
B = suplabel(['db effect: ' ctx_db_sig],'t',[0.3 0.08 0.01 0.83]);
set(B,'FontSize',12)
end


subplot('Position',[x2 y w h ])
[pyrbar] = UCSF_graph([pyrM(1:2,2),pyrM(3:4,2)]',[pyrM(1:2,1),pyrM(3:4,1)]',pyrC);

title('Pyramidal');
legend('db/+','db/db');
legend('boxoff')
legend('Location','northoutside')
pyr_age_sig = sig_check(pyrP(2));
pyr_db_sig = sig_check(pyrP(1));
if pyrP(2) <= 0.05
C = suplabel(['age effect: ' pyr_age_sig],'t',[0.6 0.08 0.01 0.83]);
set(C,'FontSize',12)
end
if pyrP(1) <= 0.05
D = suplabel(['db effect: ' pyr_db_sig],'t',[0.6 0.08 0.01 0.78]);
set(D,'FontSize',12)
end

subplot('Position',[x3 y w h ])
[pyrbar] = UCSF_graph([slmM(1:2,2),slmM(3:4,2)]',[slmM(1:2,1),slmM(3:4,1)]',slmC);

title('SLM');
legend('db/+','db/db');
legend('boxoff')
legend('Location','northoutside')
slm_db_sig = sig_check(slmP(1));
slm_age_sig = sig_check(slmP(2));
if slmP(2) <= 0.05
E = suplabel(['age effect: ' slm_age_sig],'t',[0.91 0.08 0.01 0.83]);
set(E,'FontSize',12)
end
if slmP(1) <= 0.05
F = suplabel(['db effect: ' slm_db_sig],'t',[0.91 0.08 0.01 0.78]);
set(F,'FontSize',12)
end

cd('C:\Users\ipzach\Documents\dbdb electrophy\General_Scripts')

function sig = sig_check(P)
if P > 0.05
    sig = 'n.s.';
elseif P <= 0.05 && P >0.01
    sig = '*';
elseif P <= 0.01 && P >0.001
    sig = '**';
elseif P <= 0.001
    sig = '***';
end
end
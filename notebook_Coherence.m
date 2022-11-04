% load in spkinf and concise channel info
 clear all; close all;clc
load('SpkInfo.mat')
load('chans.mat')
% file path for 
%filepath = 'C:\Users\ipzach\Documents\dbdb electrophy';
filepath ='C:\COM\ePhy\dbdb\Data\dbdb electrophy';

cd(filepath)
animalList = dir;
Fs = 1250; % Sampling Frequency; needed for filtering and plotting
% Theta/Delta state to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%
% Group, Band, recording, Animal, Layer/layer
Co = NaN(4,7,2,25,3,3);
%% Stroke group to analyze
% 1:4 6 ref Spk Info
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
        for i = 1:2
            try
                load(['LFP_' num2str(i) '.mat'])
                disp(['loaded ' pwd '\LFP_' num2str(i) '.mat'])
            catch
                disp(['could not load ' pwd '\LFP_' num2str(i) '.mat'])
                continue
            end
            LFP = LFPs{1,2};
            for layer = 1:3
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
                
                A_LFP = LFP(:, chans(A,animal));
                B_LFP = LFP(:, chans(B,animal));
                % create a vector of indiviudal frequencies to calculate
                % coherence
                for iBand = 1:7
                    switch iBand
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
                    Co(group,iBand,i,animal-2,A,B) = nanmean(mscohere(A_LFP,B_LFP,hamming(12500),[],range,1250));
                end % frequency band
            end % layer
        end  % LFP
        cd ..
    end % animal
    
end % group
%% Clean Data
% Group, Band, recording, Animal, Layer, layer
% (4,4,2,25,3,3);
count = 0;
lab1 = {'Ctx','Pyr'};
lab2 = {'Pyr','SLM'};
for i = 1:4
    switch i
        case 1
            treatName = 'DB+ 200';
        case 2
            treatName = 'DB+ 400';
        case 3
            treatName = 'DBDB 200';
        case 4
            treatName = 'DBDB 400';
    end
    
    for j = 1:7
        switch j
            case 1
                groupName = 'Delta';
            case 2
                groupName = 'Theta';
            case 3
                groupName = 'Alpha';
            case 4
                groupName = 'Beta';
            case 5
                groupName = 'Gamma';
            case 6
                groupName = 'High Gamma';
            case 7
                groupName = 'Full';
        end
        count = count +1;
        %average the coherence value for each recording session per animal
        % Group, Band, recording, Animal, Layer/layer
        combineRec = squeeze(nanmean(Co,3));
        
        % average coherence values of all animals to make a single plot
        % Group, Band, Animal, Layer, layer
        combineAnimal = squeeze(nanmean(combineRec,3));
        % group, band, layer, layer
        %subplot(4,7,count)
        %h = heatmap(squeeze(combineAnimal(i,j,1:2,2:3)));
        %h.ColorbarVisible = 'off';
        %set(gca,'xData',lab2,'yData',lab1);
        %caxis([ 0 1]);
        %if i ==1
         %   title(groupName)
        %end
        %if j ==1
        %    ylabel(treatName)
        %end
        
    end
end


%% Stats
% Group, Band, recording, Animal, Layer, layer
% We will compare across group, combining recordings, and considering all
% animals, so we need two for loops, one for frequency ban, and one for
% each (3) layer combination
% We already combined the recordings in combineRec
% combineRec(Group, Band, Animal, Layer, layer)
figure
set(gcf, 'Color','w','Position',[100 100 1300 370])
counter = 0;
for layComb = 1:2
    for band = 1:7
        switch band
            case 1
                groupName = 'Delta';
            case 2
                groupName = 'Theta';
            case 3
                groupName = 'Alpha';
            case 4
                groupName = 'Beta';
            case 5
                groupName = 'Gamma';
            case 6
                groupName = 'High Gamma';
            case 7
                groupName = 'Full';
        end
        
        % First we want to grab individual values, create 2-way labels for
        % them, then concatenate everything together
        switch layComb
            case 1
                Ct200 = squeeze(combineRec(1, band, :, 1, 2));
                Ct400 = squeeze(combineRec(2, band, :, 1, 2));
                Db200 = squeeze(combineRec(3, band, :, 1, 2));
                Db400 = squeeze(combineRec(4, band, :, 1, 2));
            case 2
                Ct200 = squeeze(combineRec(1, band, :, 1, 3));
                Ct400 = squeeze(combineRec(2, band, :, 1, 3));
                Db200 = squeeze(combineRec(3, band, :, 1, 3));
                Db400 = squeeze(combineRec(4, band, :, 1, 3));
            case 3
                Ct200 = squeeze(combineRec(1, band, :, 2, 3));
                Ct400 = squeeze(combineRec(2, band, :, 2, 3));
                Db200 = squeeze(combineRec(3, band, :, 2, 3));
                Db400 = squeeze(combineRec(4, band, :, 2, 3));
        end % switch layComb
        
        Ct200(isnan(Ct200)) = [];
        Ct400(isnan(Ct400)) = [];   
        Db200(isnan(Db200)) = [];
        Db400(isnan(Db400)) = [];
        
        Ct2AgeLab = cell(length(Ct200),1);
        Ct2AgeLab(:) = {'200'};
        Ct2DbLab = cell(length(Ct200),1);
        Ct2DbLab(:) = {'Ctrl'};
        
        Ct4AgeLab = cell(length(Ct400),1);
        Ct4AgeLab(:) = {'400'};
        Ct4DbLab = cell(length(Ct400),1);
        Ct4DbLab(:) = {'Ctrl'};
        
        Db2AgeLab = cell(length(Db200),1);
        Db2AgeLab(:) = {'200'};
        Db2DbLab = cell(length(Db200),1);
        Db2DbLab(:) = {'DBDB'};
        
        Db4AgeLab = cell(length(Db400),1);
        Db4AgeLab(:) = {'400'};
        Db4DbLab = cell(length(Db400),1);
        Db4DbLab(:) = {'DBDB'};
        
        vals = [Ct200; Db200; Ct400; Db400];
        ageLabs = [Ct2AgeLab; Db2AgeLab;Ct4AgeLab; Db4AgeLab];
        dbLabs = [Ct2DbLab; Db2DbLab; Ct4DbLab; Db4DbLab];
        
        [cohP,cohT,cohStats] = anovan(vals,{dbLabs ageLabs},'model','interaction','display','off');
        [cohC,cohM,~,cohN] = multcompare(cohStats,'Dimension',[1 2],'CType','bonferroni','display','off');
        %          if layComb == 1
        
        counter = counter +1;
        subplot(2,7,counter)
        UCSF_graph([cohM(1:2,2),cohM(3:4,2)]',[cohM(1:2,1),cohM(3:4,1)]',cohC);
        %MS
        T_Chrnc = cohM'; 
        Datetime_Chrnc = string(datetime('now'));
        cd ('C:\COM\ePhy\dbdb\Data\Outputs\Data\Coherence_Notebook')
        Filename_Chrnc = sprintf('Coherence_data_animal_%d_%s_%s.xlsx', animal, compare, Datetime_Chrnc);
        Filename_Chrnc = regexprep(Filename_Chrnc, ' ', '_');
        Filename_Chrnc = regexprep(Filename_Chrnc, ':', '_');
        xlswrite(Filename_Chrnc,T_Chrnc);
        %ME
        
        set(gca, 'ytick', [0 0.5 1])
        ylim([0 1.1])
        %             if counter == 1
        %                 ylabel('Coherence')
        %             end
        %             if counter == 5
        %                 legend('db/+','db/db');
        %                 legend('boxoff')
        %                 legend('Location',[-0.05 0.1 0.2 0.2])
        %             end
        
        age_sig = sig_check(cohP(1));
        db_sig = sig_check(cohP(2));
        %title(['age effect: ' age_sig ' db effect: ' db_sig])
        disp([num2str(layComb) ' ' groupName ' age:' num2str(age_sig) ' db:' num2str(db_sig)])
        
        
        if band == 1 
            if layComb == 1
            ylabel({'Cortical-Pyramidal','Coherence'})
            elseif layComb == 2
                ylabel({'Cortical-SLM','Coherence'})
            end
        end
        if layComb == 1
            title(groupName)
%             B = suplabel(['age effect: ' age_sig],'t',[0.065+(0.05*counter) 0.08 0.01 0.78]);
%             C = suplabel(['db effect: ' db_sig],'t',[0.065+(0.05*counter) 0.08 0.01 0.73]);
%             set(B,'FontSize',12)
%             set(C,'FontSize',12)
%         elseif  layComb == 2
%             B = suplabel(['age effect: ' age_sig],'t',[0.065+(0.05*(counter-7)) 0.08 0.01 0.4]);
%             C = suplabel(['db effect: ' db_sig],'t',[0.065+(0.05*(counter-7)) 0.08 0.01 0.35]);
%             set(B,'FontSize',12)
%             set(C,'FontSize',12)
            
        end
        %end
        %end
    end % for layComb
end% for band

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



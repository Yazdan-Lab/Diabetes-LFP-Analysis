disp('Initializing Files')
load('SpkInfo.mat')
load('chans.mat')
filepath = 'C:\Users\ipzach\Documents\dbdb electrophy';
cd(filepath)
animalList = dir;
Fs = 1250; % Sampling Frequency; needed for filtering and plotting
% Theta/Delta state to analyze
TD = 3; %High 1, Low 2, Full 3
%%%%%%%%%%%%%%%%%%%%%%%%%
% Group, Band, recording, Animal, Layer/layer
Co = zeros(4,4,2,25,3,3);
%% Stroke group to analyze
% 1:4 6 ref Spk Info
for group = 1:4
    % 1.Load the right and left side signals
    %load H/L TD indexes
    
    if group ==1
        grouping = 3:9; % DB+ 200D
    elseif group ==2
        grouping = 10:14; % DB+ 400D
    elseif group ==3
        grouping = [15:18 20 21]; % DBDB 200D
    elseif group ==4
        grouping = [22 24:27]; % DBDB 400D
    end
    
    
    for animal = grouping
        disp(animal)
        cd(animalList(animal).name)
%         load('REM.mat')
%         try
%             test = rem(1).R.start(1);
%             test2 = rem(2).R.start(2);
%         catch
%             disp('REM (1) or (2) start does not exist')
%             cd ..
%             continue
%         end
        for i = 1:2
            load(['LFP_' num2str(i) '.mat'])
            disp(['loaded ' pwd '\LFP_' num2str(i) '.mat'])
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
                
                for iBand = 1:5
                    switch iBand
                        case 1
                            range = linspace(0.1,3,20);
                        case 2
                            range = linspace(4,7,20);
                        case 3
                            range = linspace(30,58,20);
                        case 4
                            range = linspace(62,200,20);
                        case 5
                            range = linspace(0,200,50);
                    end % switch iBand
                    % Group, Band, recording, Animal, Layer/layer
                    Co(group,iBand,i,animal-2,A,B) = mean(mscohere(A_LFP,B_LFP,hamming(12500),[],range,1250));
                end % frequency band
            end % layer
        end  % LFP
        cd ..
    end % animal
    
end % group
%% Plot summary figures
% Group, Band, recording, Animal, Layer, layer
% (4,4,2,25,3,3);
Co(Co == 0) = NaN;
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
    
    for j = 1:5
        switch j
            case 1
                groupName = 'Delta';
            case 2
                groupName = 'Theta';
            case 3
                groupName = 'Gamma';
            case 4
                groupName = 'HGamma';
            case 5
                groupName = 'Full';
        end
        count = count +1;
        combineRec = squeeze(nanmean(Co,3));
        % Group, Band, Animal, Layer, layer
        combineAnimal = squeeze(nanmean(combineRec,3));
        % group, band, layer, layer
        subplot(4,5,count)
        h = heatmap(squeeze(combineAnimal(i,j,1:2,2:3)));
        h.ColorbarVisible = 'off';
        set(gca,'xData',lab2,'yData',lab1);
        caxis([ 0 1]);
        if i ==1
            title(groupName)
        end
        if j ==1
            ylabel(treatName)
        end
        
    end
end


%% Stats
% Group, Band, recording, Animal, Layer, layer
% We will compare across group, combining recordings, and considering all
% animals, so we need two for loops, one for frequency ban, and one for
% each (3) layer combination
% We already combined the recordings in combineRec
% combineRec(Group, Band, Animal, Layer, layer)
for band = 1:5
    for layComb = 1:3
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
        vals = [Ct200; Ct400; Db200; Db400];
        ageLabs = [Ct2AgeLab; Ct4AgeLab; Db2AgeLab; Db4AgeLab];
        dbLabs = [Ct2DbLab; Ct4DbLab; Db2DbLab; Db4DbLab];
        
        [~,cohT,cohStats] = anovan(vals,{ageLabs,dbLabs},'model','interaction');
        [cohC,cohM,~,cohN] = multcompare(cohStats,'Dimension',[1 2],'CType','bonferroni');
        % [C,M] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
    end % for layComb
end% for band




















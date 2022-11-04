disp('Initializing Files')
clear all, close all
load('SpkInfo.mat')
load('chans.mat')
%filepath = 'C:\Users\ipzach\Documents\dbdb electrophy';
filepath = 'C:\COM\ePhy\dbdb\Data\dbdb electrophy';
cd(filepath)
animalList = dir;
% THIS CODE ASSIGNS GROUPS WRONG
Fs = 1250; % Sampling Frequency; needed for filtering and plotting
% Theta/Delta state to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%
% Group, Band, recording, Animal, Layer/layer
scores = NaN(4,7,2,3);
%% Stroke group to analyze
% 1:4 6 ref Spk Info
for group = 1:4
    % 1.Load the right and left side signals
    %load H/L TD indexes
    % Commented out to break code until fixed
    if group ==1
        grouping = 3:9; % DB+ 200D
    elseif group ==2
        grouping = [15:18 20 21]; % DBDB 200D
    elseif group ==3
        grouping = 10:14; % DB+ 400D
    elseif group ==4
        grouping = [22 24:27]; % DBDB 400D
    end
    
    count = 0;
    for animal = grouping
        count = count +1;
        disp(animal)
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
                single_channel = LFP(:, chans(layer,animal));
                low_freq = BPfilter(single_channel, 1250, 1, 8);
                high_freq = BPfilter(single_channel, 1250, 9, 30);
                scores(group,count,i,layer) = SignalPower( low_freq,1250) ./ SignalPower(high_freq,1250);
            end % layer
        end % recording
        
        cd ..
    end % animal
end  % group

scores_combine = squeeze(nanmean(scores,3));
%% Plot summary figures
count = 0;

%boxplot_grouping =repmat(boxplot_grouping,7,1);
score_labels = {'Ct200','DB200','Ct400','DB400'};
for i = 1:3
    figure
    temp = scores_combine(:,:,i)';
    boxplot(temp,score_labels')
    %MS
    Datetime_SlwScr = string(datetime('now'));
    Filename_SlwScr = sprintf('SlowingScore_Figure_%s.tiff', Datetime_SlwScr);
    Filename_SlwScr = regexprep(Filename_SlwScr, ' ', '_');
    Filename_SlwScr = regexprep(Filename_SlwScr, ':', '_');
    saveas(gcf, Filename_SlwScr);
    %ME
    pause(1)
end


%% Stats
% Group, Band, recording, Animal, Layer, layer
% We will compare across group, combining recordings, and considering all
% animals, so we need two for loops, one for frequency ban, and one for
% each (3) layer combination
% We already combined the recordings in combineRec
% combineRec(Group, Band, Animal, Layer, layer)
figure
set(gcf, 'color', 'w','Position',[100 100 800 420])
for layComb = 2:3
    % First we want to grab individual values, create 2-way labels for
    % them, then concatenate everything together
    temp = scores_combine(:,:,layComb)';
    
    Ct200 = temp(:,1);
    Db200 = temp(1:6,2);
    Ct400 = temp(1:5,3);
    Db400 = temp(1:5,4);
    
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
    %MS
    Group_Ns = [length(Ct200); length(Db200); length(Ct400); length(Db400)];
    %ME
    ageLabs = [Ct2AgeLab; Db2AgeLab;Ct4AgeLab; Db4AgeLab];
    dbLabs = [Ct2DbLab; Db2DbLab; Ct4DbLab; Db4DbLab];
    
    [ssP,ssT,ssStats] = anovan(vals,{dbLabs ageLabs},'model','interaction','display','off');
    [ssC,ssM,~,ssN] = multcompare(ssStats,'Dimension',[1 2],'CType','bonferroni','display','off');
        
    subplot(1,2,layComb-1)
    UCSF_graph([ssM(1:2,2),ssM(3:4,2)]',[ssM(1:2,1),ssM(3:4,1)]',ssC);
    %MS
    T_SlwScr = ssM';
    T_SlwScr = [T_SlwScr;Group_Ns'];
    Datetime_SlwScr = string(datetime('now'));
    cd('C:\COM\ePhy\dbdb\Data\Outputs\Data\SlowingScore')
    Filename_SlwScr = sprintf('SlowingScore_data_%d_%s.xlsx', layComb, Datetime_SlwScr);
    Filename_SlwScr = regexprep(Filename_SlwScr, ' ', '_');
    Filename_SlwScr = regexprep(Filename_SlwScr, ':', '_');
    xlswrite(Filename_SlwScr,T_SlwScr);
    %ME
    set(gca,'Ytick',[0 5 10 15],'fontweight','bold')
    ylim([0 15])
    switch layComb
        case 1
            title('Cortex')
        case 2
            title('Pyramidal','fontweight','bold')
            
            ylabel('Slowing score','fontweight','bold')
        case 3
            title('SLM','fontweight','bold')
    end
    
    
    age_sig = sig_check(ssP(2));
    db_sig = sig_check(ssP(1));
    if ssP(2) <= 0.05
        B = suplabel(['age effect: ' age_sig],'t',[0.43*(layComb-1) 0.08 0.01 0.8]);
        set(B,'FontSize',12,'FontWeight','bold')
    end
    if ssP(1) <= 0.05
        C = suplabel(['db effect: ' db_sig],'t',[0.43*(layComb-1) 0.08 0.01 0.75]);
        set(C,'FontSize',12,'FontWeight','bold')
    end
    %MS
    Datetime_SlwScr = string(datetime('now'));
    Filename_SlwScr = sprintf('SlowingScore_Figure_%s.tiff', Datetime_SlwScr);
    Filename_SlwScr = regexprep(Filename_SlwScr, ' ', '_');
    Filename_SlwScr = regexprep(Filename_SlwScr, ':', '_');
    saveas(gcf, Filename_SlwScr);
    %ME
end

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



user = 'Z'; %'Z' for Zach or 'S' for Shahram for path stuff

switch user
    case 'Z'
        cd('C:\Users\ipzach\Documents\MATLAB\output\Diabetes-Saved-Files')
    case 'S'
        cd('C:\COM\ePhy\dbdb\Data\Outputs\Data')
end

load('LFP Measures')
hotcold = redblue();

%% Make labels for stats
% We need two labels for a 2-way anova test, one for each variable we are
% testing
label.DB2_age = cell(size(Gamma.DB2_Pyr, 2), 1);
label.DB2_treat = cell(size(Gamma.DB2_Pyr, 2), 1);
label.DB2_age(:) = {'200'};
label.DB2_treat(:) = {'Control'};

label.DB4_age = cell(size(Gamma.DB4_Pyr, 2), 1);
label.DB4_treat = cell(size(Gamma.DB4_Pyr, 2), 1);
label.DB4_age(:) = {'400'};
label.DB4_treat(:) = {'Control'};

label.DBDB2_age = cell(size(Gamma.DBDB2_Pyr, 2), 1);
label.DBDB2_treat = cell(size(Gamma.DBDB2_Pyr, 2), 1);
label.DBDB2_age(:) = {'200'};
label.DBDB2_treat(:) = {'DBDB'};

label.DBDB4_age = cell(size(Gamma.DBDB4_Pyr, 2), 1);
label.DBDB4_treat = cell(size(Gamma.DBDB4_Pyr, 2), 1);
label.DBDB4_age(:) = {'400'};
label.DBDB4_treat(:) = {'DBDB'};

% extra set of labels for groups without removed values for windowing
label.r_DB2_age = cell(length(rip.DB2), 1);
label.r_DB2_treat = cell(length(rip.DB2), 1);
label.r_DB2_age(:) = {'200'};
label.r_DB2_treat(:) = {'Control'};

label.r_DB4_age = cell(length(rip.DB4), 1);
label.r_DB4_treat = cell(length(rip.DB4), 1);
label.r_DB4_age(:) = {'400'};
label.r_DB4_treat(:) = {'Control'};

label.r_DBDB2_age = cell(length(rip.DBDB2), 1);
label.r_DBDB2_treat = cell(length(rip.DBDB2), 1);
label.r_DBDB2_age(:) = {'200'};
label.r_DBDB2_treat(:) = {'DBDB'};

label.r_DBDB4_age = cell(length(rip.DBDB4), 1);
label.r_DBDB4_treat = cell(length(rip.DBDB4), 1);
label.r_DBDB4_age(:) = {'400'};
label.r_DBDB4_treat(:) = {'DBDB'};

age_Labs = [label.DB2_age; label.DBDB2_age; label.DB4_age; label.DBDB4_age];
treat_Labs = [label.DB2_treat; label.DBDB2_treat; label.DB4_treat; label.DBDB4_treat];

r_age_Labs = [label.r_DB2_age; label.r_DBDB2_age; label.r_DB4_age; label.r_DBDB4_age];
r_treat_Labs = [label.r_DB2_treat; label.r_DBDB2_treat; label.r_DB4_treat; label.r_DBDB4_treat];

slowing_score_Ct2AgeLab = cell(sum(~isnan(slowing_score(1, :, 1))), 1);
slowing_score_Ct2AgeLab(:) = {'200'};
slowing_score_Ct2DbLab = cell(sum(~isnan(slowing_score(1, :, 1))), 1);
slowing_score_Ct2DbLab(:) = {'Ctrl'};

slowing_score_Ct4AgeLab = cell(sum(~isnan(slowing_score(2, :, 1))), 1);
slowing_score_Ct4AgeLab(:) = {'400'};
slowing_score_Ct4DbLab = cell(sum(~isnan(slowing_score(2, :, 1))), 1);
slowing_score_Ct4DbLab(:) = {'Ctrl'};

slowing_score_Db2AgeLab = cell(sum(~isnan(slowing_score(3, :, 1))), 1);
slowing_score_Db2AgeLab(:) = {'200'};
slowing_score_Db2DbLab = cell(sum(~isnan(slowing_score(3, :, 1))), 1);
slowing_score_Db2DbLab(:) = {'DBDB'};

slowing_score_Db4AgeLab = cell(sum(~isnan(slowing_score(4, :, 1))), 1);
slowing_score_Db4AgeLab(:) = {'400'};
slowing_score_Db4DbLab = cell(sum(~isnan(slowing_score(4, :, 1))), 1);
slowing_score_Db4DbLab(:) = {'DBDB'};


slowing_score_age_Labs = [slowing_score_Ct2AgeLab; slowing_score_Db2AgeLab; slowing_score_Ct4AgeLab; slowing_score_Db4AgeLab];
slowing_score_db_Labs = [slowing_score_Ct2DbLab; slowing_score_Db2DbLab; slowing_score_Ct4DbLab; slowing_score_Db4DbLab];


state_changes_Ct2AgeLab = cell(sum(~isnan(state_changes(1, :))), 1);
state_changes_Ct2AgeLab(:) = {'200'};
state_changes_Ct2DbLab = cell(sum(~isnan(state_changes(1, :))), 1);
state_changes_Ct2DbLab(:) = {'Ctrl'};

state_changes_Ct4AgeLab = cell(sum(~isnan(state_changes(2, :))), 1);
state_changes_Ct4AgeLab(:) = {'400'};
state_changes_Ct4DbLab = cell(sum(~isnan(state_changes(2, :))), 1);
state_changes_Ct4DbLab(:) = {'Ctrl'};

state_changes_Db2AgeLab = cell(sum(~isnan(state_changes(3, :))), 1);
state_changes_Db2AgeLab(:) = {'200'};
state_changes_Db2DbLab = cell(sum(~isnan(state_changes(3, :))), 1);
state_changes_Db2DbLab(:) = {'DBDB'};

state_changes_Db4AgeLab = cell(sum(~isnan(state_changes(4, :))), 1);
state_changes_Db4AgeLab(:) = {'400'};
state_changes_Db4DbLab = cell(sum(~isnan(state_changes(4, :))), 1);
state_changes_Db4DbLab(:) = {'DBDB'};


state_changes_age_Labs = [state_changes_Ct2AgeLab; state_changes_Db2AgeLab; state_changes_Ct4AgeLab; state_changes_Db4AgeLab];
state_changes_db_Labs = [state_changes_Ct2DbLab; state_changes_Db2DbLab; state_changes_Ct4DbLab; state_changes_Db4DbLab];

clear per_animal_Ct2AgeLab per_animal_Ct2DbLab per_animal_Ct4AgeLab per_animal_Ct4DbLab per_animal_Db2AgeLab per_animal_Db2DbLab per_animal_Db4AgeLab per_animal_Db4DbLab
% Process CSD into single value
CSD.DB2_max = squeeze(max(CSD.DB2, [], [1, 2])); %squeeze(CSD.DB2(625:1500,4,:));
CSD.DB4_max = squeeze(max(CSD.DB4, [], [1, 2])); %squeeze(CSD.DB4(625:1500,4,:));
CSD.DBDB2_max = squeeze(max(CSD.DBDB2, [], [1, 2])); %squeeze(CSD.DBDB2(625:1500,3,:));
CSD.DBDB4_max = squeeze(max(CSD.DBDB4, [], [1, 2])); %squeeze(CSD.DBDB4(625:1500,4,:));

CSD.DB2_min = squeeze(min(min(CSD.DB2, [], 1), [], 2)); %squeeze(CSD.DB2(625:1500,3,:));
CSD.DB4_min = squeeze(min(min(CSD.DB4, [], 1), [], 2)); %squeeze(CSD.DB4(625:1500,3,:));
CSD.DBDB2_min = squeeze(min(min(CSD.DBDB2, [], 1), [], 2)); %squeeze(CSD.DBDB2(625:1500,2,:));
CSD.DBDB4_min = squeeze(min(min(CSD.DBDB4, [], 1), [], 2)); %squeeze(CSD.DBDB4(625:1500,3,:));

CSD.DB2_full_amp = CSD.DB2_max - CSD.DB2_min;
CSD.DB4_full_amp = CSD.DB4_max - CSD.DB4_min;
CSD.DBDB2_full_amp = CSD.DBDB2_max - CSD.DBDB2_min;
CSD.DBDB4_full_amp = CSD.DBDB4_max - CSD.DBDB4_min;

% Try a more specific window
CSD.DB2_rip = mean(squeeze(CSD.DB2(650:1400, 4, :)), 1)';
CSD.DB4_rip = mean(squeeze(CSD.DB4(650:1400, 4, :)), 1)';
CSD.DBDB2_rip = mean(squeeze(CSD.DBDB2(650:1400, 4, :)), 1)';
CSD.DBDB4_rip = mean(squeeze(CSD.DBDB4(650:1400, 4, :)), 1)';

CSD.DB2_wav = mean(squeeze(CSD.DB2(650:1400, 2, :)), 1)';
CSD.DB4_wav = mean(squeeze(CSD.DB4(650:1400, 2, :)), 1)';
CSD.DBDB2_wav = mean(squeeze(CSD.DBDB2(650:1400, 2, :)), 1)';
CSD.DBDB4_wav = mean(squeeze(CSD.DBDB4(650:1400, 2, :)), 1)';

CSD.DB2_amp = CSD.DB2_rip - CSD.DB2_wav;
CSD.DB4_amp = CSD.DB4_rip - CSD.DB4_wav;
CSD.DBDB2_amp = CSD.DBDB2_rip - CSD.DBDB2_wav;
CSD.DBDB4_amp = CSD.DBDB4_rip - CSD.DBDB4_wav;

%%
IRI_vals = [];
IRI_age = {};
IRI_treat = {};
IRI_big = NaN(4, 2000);
for l = [1, 3, 2, 4]
    switch l
        case 1
            group = rip.DB2(:, 1);
            age = '200';
            treat = 'Control';
        case 2
            group = rip.DB4(:, 1);
            age = '400';
            treat = 'Control';
        case 3
            group = rip.DBDB2(:, 1);
            age = '200';
            treat = 'DBDB';
        case 4
            group = rip.DBDB4(:, 1);
            age = '400';
            treat = 'DBDB';
    end
    for m = 1:length(group) - 1
        if group(m+1) > group(m)
            IRI_vals = [IRI_vals; (group(m+1) - group(m))];
            IRI_age = [IRI_age; {age}];
            IRI_treat = [IRI_treat; {treat}];
            IRI_big(l, m) = (group(m+1) - group(m));
            
        end
    end
end
IRIdb2 = rmoutliers(IRI_big(1, :));
IRIdb4 = rmoutliers(IRI_big(2, :));
IRIdbdb2 = rmoutliers(IRI_big(3, :));
IRIdbdb4 = rmoutliers(IRI_big(4, :));

IRIdb2 = IRIdb2(~isnan(IRIdb2)) ./ 1250;
IRIdb4 = IRIdb4(~isnan(IRIdb4)) ./ 1250;
IRIdbdb2 = IRIdbdb2(~isnan(IRIdbdb2)) ./ 1250;
IRIdbdb4 = IRIdbdb4(~isnan(IRIdbdb4)) ./ 1250;

[cleanIRI, TF] = rmoutliers(IRI_vals, 'quartiles');
IRI_age(TF == 1) = [];
IRI_treat(TF == 1) = [];

%% Slowing score
disp('Slowing score')

figure
set(gcf, 'Position', [100, 100, 800, 400])
for lay_comb = 2:3
    % First we want to grab individual values, create 2-way labels for
    % them, then concatenate everything together
    SS_Ct200_w_nan = slowing_score(1, :, lay_comb)';
    SS_Ct200 = SS_Ct200_w_nan(~isnan(SS_Ct200_w_nan));
    
    SS_DB200_w_nan = slowing_score(2, :, lay_comb)';
    SS_DB200 = SS_DB200_w_nan(~isnan(SS_DB200_w_nan));
    
    SS_Ct400_w_nan = slowing_score(3, :, lay_comb)';
    SS_Ct400 = SS_Ct400_w_nan(~isnan(SS_Ct400_w_nan));
    
    SS_DB400_w_nan = slowing_score(4, :, lay_comb)';
    SS_DB400 = SS_DB400_w_nan(~isnan(SS_DB400_w_nan));
    
    slow_score_vals = [SS_Ct200; SS_DB200; SS_Ct400; SS_DB400];
    
    disp(num2str(lay_comb))
    [ssP, ssT, ssStats] = anovan(slow_score_vals, {slowing_score_db_Labs, slowing_score_age_Labs}, 'model', 'interaction', 'display', 'off');
    [ssC, ssM, ~, ssN] = multcompare(ssStats, 'Dimension', [1, 2], 'CType', 'bonferroni', 'display', 'off');
    subplot(1, 2, lay_comb-1)
    create_bar_figure(ssM(:, 2), ssM(:, 1), ssC);
    %sig_values(ssP(2), ssP(1));
    %ylabel('Slowing Score')
    set(gcf, 'Color', 'w');
    set(gca, 'ytick', [0, 6, 12])
    ylim([0, 15])
    xtickangle(60)
    %     switch lay_comb
    %         case 1
    %             title('Cortex')
    %         case 2
    %             title('Pyramidal')
    %         case 3
    %             title('SLM')
    %     end
end

%% Spectral exponent

% First we want to grab individual values, create 2-way labels for
% them, then concatenate everything together
idx = ~cellfun('isempty', intSlo_Store);
SE = NaN(size(intSlo_Store));
SE(idx) = cellfun(@(v)v(2), intSlo0_Store(idx));

SE_Ct200_w_nan = SE(1, :)';
SE_Ct200 = SE_Ct200_w_nan(~isnan(SE_Ct200_w_nan));

SE_DB200_w_nan = SE(2, :)';
SE_DB200 = SE_DB200_w_nan(~isnan(SE_DB200_w_nan));

SE_Ct400_w_nan = SE(3, :)';
SE_Ct400 = SE_Ct400_w_nan(~isnan(SE_Ct400_w_nan));

SE_DB400_w_nan = SE(4, :)';
SE_DB400 = SE_DB400_w_nan(~isnan(SE_DB400_w_nan));

SE_vals = [SE_Ct200; SE_DB200; SE_Ct400; SE_DB400];

%MS
Group_SP_Ns = [length(SE_Ct200); length(SE_DB200); length(SE_Ct400); length(SE_DB400)];
%ME

[seP, seT, seStats] = anovan(SE_vals, {slowing_score_db_Labs, slowing_score_age_Labs}, 'model', 'interaction', 'display', 'off');
[seC, seM, ~, seN] = multcompare(seStats, 'Dimension', [1, 2], 'CType', 'bonferroni', 'display', 'off');

%set(gca,'ytick',[0 6 12])
%ylim([0 15])

%% Spectral exponent figure
XXi = Pows_store(1, 1).frex;
YYi = NaN(7, 516);
int_Slo = NaN(7, 2);
figure
set(gcf, 'Position', [100, 100, 800, 400])
subplot(1, 2, 1)
for i = [4, 3, 2, 1]
    switch i
        case 1
            %color = [194, 194, 194] ./255;
            color = [0, 0, 204] ./ 255;
        case 2
            %color = [128, 128, 128] ./255;
            color = [204, 0, 0] ./ 255;
        case 3
            %color = [247,149,114] ./255;
            color = [0, 0, 204] ./ 255;
        case 4
            %color = [212, 120, 86] ./255;
            color = [204, 0, 0] ./ 255;
            
    end
    for j = 1:7
        try
            YYi(j, :) = Pows_store(i, j).obs;
            int_Slo(j, :) = intSlo0_Store{i, j};
        catch
            continue
        end
    end
    YYi = rmmissing(YYi);
    int = nanmean(int_Slo(:, 1));
    slo = nanmean(int_Slo(:, 2));
    Xi = log10(XXi);
    
    stdshade(log(YYi), 0.1, color, log(XXi)); hold on,
    YYpred0 = 10.^(int + slo * (Xi))';
    plot(log(XXi([1, end])), log(YYpred0([1, end])), 'LineWidth', 1.5, 'color', color);
    set(gca, 'FontSize', 14, 'TickDir', 'out');
    box off
    xlim([0, 3.6])
end
subplot(1, 2, 2)
create_bar_figure(seM(:, 2), seM(:, 1), seC);
%MS
if strcmp(user, 'S')
    T_SpecExpoB = seM';
    T_SpecExpoB = [T_SpecExpoB; Group_SP_Ns'];
    Datetime_SpecExpoB = string(datetime('now'));
    
    cd('C:\COM\ePhy\dbdb\Data\Outputs\Data\SpecExpoB')
    Filename_SpecExpoB = sprintf('Spectral_Exponent_data_%s.xlsx', Datetime_SpecExpoB);
    Filename_SpecExpoB = regexprep(Filename_SpecExpoB, ' ', '_');
    Filename_SpecExpoB = regexprep(Filename_SpecExpoB, ':', '_');
    xlswrite(Filename_SpecExpoB, T_SpecExpoB);
end
%ME
xtickangle(60)

%sig_values(seP(2), seP(1));
%ylabel('Spectral Exponent')
set(gcf, 'Color', 'w');
%MS
if strcmp(user, 'S')
    Datetime_SpecExpoB = string(datetime('now'));
    cd('C:\COM\ePhy\dbdb\Data\Outputs\Data\SpecExpoB')
    Filename_SpecExpoB = sprintf('Spectral_Exponent_Figure_%s.tiff', Datetime_SpecExpoB);
    Filename_SpecExpoB = regexprep(Filename_SpecExpoB, ' ', '_');
    Filename_SpecExpoB = regexprep(Filename_SpecExpoB, ':', '_');
    saveas(gcf, Filename_SpecExpoB);
end
%ME

%% PLI
% group,animal, band, layer
figure
set(gcf, 'Position', [100, 100, 1300, 350])
for lay_comb = 1
    switch lay_comb
        case 1
            comb_name = 'Ctx-Pyr';
        case 2
            comb_name = 'Ctx-Slm';
        case 3
            comb_name = 'Pyr-Slm';
    end % switch layComb
    %im bored
    for band = 1:5
        switch band
            case 1
                group_name = 'Delta ';
            case 2
                group_name = 'Theta ';
            case 3
                group_name = 'Alpha ';
            case 4
                group_name = 'Beta ';
            case 5
                group_name = 'Gamma ';
            case 6
                group_name = 'High Gamma ';
            case 7
                group_name = 'Full ';
        end
        
        % First we want to grab individual values, create 2-way labels for
        % them, then concatenate everything together
        PLI_Ct200_w_nan = PLI(1, :, band, lay_comb)';
        PLI_Ct200 = PLI_Ct200_w_nan(~isnan(PLI_Ct200_w_nan));
        
        PLI_DB200_w_nan = PLI(2, :, band, lay_comb)';
        PLI_DB200 = PLI_DB200_w_nan(~isnan(PLI_DB200_w_nan));
        
        PLI_Ct400_w_nan = PLI(3, :, band, lay_comb)';
        PLI_Ct400 = PLI_Ct400_w_nan(~isnan(PLI_Ct400_w_nan));
        
        PLI_DB400_w_nan = PLI(4, :, band, lay_comb)';
        PLI_DB400 = PLI_DB400_w_nan(~isnan(PLI_DB400_w_nan));
        
        PLI_vals = [PLI_Ct200; PLI_DB200; PLI_Ct400; PLI_DB400];
        %MS
        Group_PLI_Ns = [length(PLI_Ct200); length(PLI_DB200); length(PLI_Ct400); length(PLI_DB400)];
        %ME
        [pliP, pliT, pliStats] = anovan(PLI_vals, {slowing_score_db_Labs, slowing_score_age_Labs}, 'model', 'interaction', 'display', 'off');
        [pliC, pliM, ~, pliN] = multcompare(pliStats, 'Dimension', [1, 2], 'CType', 'bonferroni', 'display', 'off');
        subplot(1,5,band)
        create_bar_figure(pliM(:, 2), pliM(:, 1), pliC);
        title(group_name)
        %MS
        if strcmp(user, 'S')
            T_PLI = pliM';
            T_PLI = [T_PLI; Group_PLI_Ns'];
            Datetime_PLI = string(datetime('now'));
            cd('C:\COM\ePhy\dbdb\Data\Outputs\Data\PLI')
            Filename_PLI = sprintf('PLI_data_%s.xlsx', [group_name Datetime_PLI]);
            Filename_PLI = regexprep(Filename_PLI, ' ', '_');
            Filename_PLI = regexprep(Filename_PLI, ':', '_');
            xlswrite(Filename_PLI, T_PLI);
        end
        %ME
        
        
        sig_values(pliP(2), pliP(1));
        %ylabel('Phase Locking Index')
        set(gcf, 'Color', 'w');
        xtickangle(60)
        %MS
        if strcmp(user, 'S')
            Datetime_PLI = string(datetime('now'));
            cd('C:\COM\ePhy\dbdb\Data\Outputs\Data\PLI')
            Filename_PLI = sprintf('PLI_Figure_%s.tiff', Datetime_PLI);
            Filename_PLI = regexprep(Filename_PLI, ' ', '_');
            Filename_PLI = regexprep(Filename_PLI, ':', '_');
            saveas(gcf, Filename_PLI);
        end
        %ME
        %set(gca,'ytick',[0 1])
        %ylim([0 1])
        %         switch lay_comb
        %             case 1
        %                 title([group_name comb_name])
        %             case 2
        %                 title([group_name comb_name])
        %             case 3
        %                 title([group_name comb_name])
        %         end
    end
end

%% State changes
SC_Ct200_w_nan = state_changes(1, :)';
SC_Ct200 = SC_Ct200_w_nan(~isnan(SC_Ct200_w_nan));

SC_DB200_w_nan = state_changes(2, :)';
SC_DB200 = SC_DB200_w_nan(~isnan(SC_DB200_w_nan));

SC_Ct400_w_nan = state_changes(3, :)';
SC_Ct400 = SC_Ct400_w_nan(~isnan(SC_Ct400_w_nan));

SC_DB400_w_nan = state_changes(4, :)';
SC_DB400 = SC_DB400_w_nan(~isnan(SC_DB400_w_nan));

state_changes_vals = [SC_Ct200; SC_DB200; SC_Ct400; SC_DB400];

[scP, scT, scStats] = anovan(state_changes_vals, {state_changes_db_Labs, state_changes_age_Labs}, 'model', 'interaction', 'display', 'off');
[scC, scM, ~, scN] = multcompare(scStats, 'Dimension', [1, 2], 'CType', 'bonferroni', 'display', 'off');
figure
create_bar_figure(scM(:, 2), scM(:, 1), scC);
sig_values(scP(2), scP(1));
%title('State Changes')

%% Coherence
disp('Coherence')
count = 0;
figure
set(gcf, 'Color', 'w', 'Position', [100, 100, 1400, 650])
for lay_comb = 1:2 % 1:3
    switch lay_comb
        case 1
            comb_name = 'Ctx-Pyr';
        case 2
            comb_name = 'Ctx-Slm';
        case 3
            comb_name = 'Pyr-Slm';
    end % switch layComb
    disp(comb_name)
    for band = 2:5 %1:7
        switch band
            case 1
                group_name = 'Delta';
            case 2
                group_name = 'Theta';
            case 3
                group_name = 'Alpha';
            case 4
                group_name = 'Beta';
            case 5
                group_name = 'Gamma';
            case 6
                group_name = 'High Gamma';
            case 7
                group_name = 'Full';
        end
        
        count = count + 1;
        
        Coh_Ct200_w_nan = Co(1, :, band, lay_comb)';
        Coh_Ct200 = Coh_Ct200_w_nan(~isnan(Coh_Ct200_w_nan));
        
        Coh_Ct400_w_nan = Co(2, :, band, lay_comb)';
        Coh_Ct400 = Coh_Ct400_w_nan(~isnan(Coh_Ct400_w_nan));
        
        Coh_DB200_w_nan = Co(3, :, band, lay_comb)';
        Coh_DB200 = Coh_DB200_w_nan(~isnan(Coh_DB200_w_nan));
        
        Coh_DB400_w_nan = Co(4, :, band, lay_comb)';
        Coh_DB400 = Coh_DB400_w_nan(~isnan(Coh_DB400_w_nan));
        
        % First we want to grab individual values, create 2-way labels for
        % them, then concatenate everything together
        
        coh_vals = [Coh_Ct200; Coh_DB200; Coh_Ct400; Coh_DB400];
        %MS
        Group_Coh_Ns = [length(Coh_Ct200); length(Coh_DB200); length(Coh_Ct400); length(Coh_DB400)];
        %ME
        disp(group_name)
        [cohP, cohT, cohStats] = anovan(coh_vals, {slowing_score_db_Labs, slowing_score_age_Labs}, 'model', 'interaction', 'display', 'off');
        [cohC, cohM, ~, cohN] = multcompare(cohStats, 'Dimension', [1, 2], 'CType', 'bonferroni', 'display', 'off');
        %
        subaxis(2, 4, count, 'SpacingHoriz', 0.01, 'SpacingVert', 0.12)
        create_bar_figure(cohM(:, 2), cohM(:, 1), cohC);
        %MS
        if strcmp(user, 'S')
            T_coherence = cohM';
            T_coherence = [T_coherence; Group_Coh_Ns'];
            Datetime_coherence = string(datetime('now'));
            cd('C:\COM\ePhy\dbdb\Data\Outputs\Data\Coherence')
            Filename_coherence = sprintf('Coherence_data_%s.xlsx', Datetime_coherence);
            Filename_coherence = regexprep(Filename_coherence, ' ', '_');
            Filename_coherence = regexprep(Filename_coherence, ':', '_');
            xlswrite(Filename_coherence, T_coherence);
            
            Datetime_coherence = string(datetime('now'));
            cd('C:\COM\ePhy\dbdb\Data\Outputs\Data\Coherence')
            Filename_coherence = sprintf('Coherence_Figure_%s.tiff', Datetime_coherence);
            Filename_coherence = regexprep(Filename_coherence, ' ', '_');
            Filename_coherence = regexprep(Filename_coherence, ':', '_');
            saveas(gcf, Filename_coherence);
            %ME
        end
        set(gca, 'ytick', [0, 1], 'fontsize', 12)
        %      sig_values(cohP(2), cohP(1));
        xtickangle(25)
        ylim([0, 1])
        if count < 5
            title(group_name)
        end
        
        if band == 2
            ylabel({comb_name, 'Coherence'})
        else
            set(gca, 'YColor', 'none')
        end
    end
end

%% Not being used in manuscript anymore
% disp('Power')
% power_vals = [rip.DB2(:, 7); rip.DB4(:, 7); rip.DBDB2(:, 7); rip.DBDB4(:, 7)];
%
% [powerP, powerT, power_stats] = anovan(power_vals, {r_treat_Labs, r_age_Labs}, 'model', 'interaction', 'display', 'off');
% [powerC, powerM, ~, powerNames] = multcompare(power_stats, 'Dimension', [1, 2], 'CType', 'bonferroni', 'display', 'off');
%
%% % Duration
disp('Duration')
dur_vals = [rip.DB2(:, 2) - rip.DB2(:, 1); rip.DBDB2(:, 2) - rip.DBDB2(:, 1); rip.DB4(:, 2) - rip.DB4(:, 1); rip.DBDB4(:, 2) - rip.DBDB4(:, 1)] ./ 1250;
[durP, durT, dur_stats] = anovan(dur_vals, {r_treat_Labs, r_age_Labs}, 'model', 'interaction', 'display', 'off');
[durC, durM, ~, durNames] = multcompare(dur_stats, 'Dimension', [1, 2], 'CType', 'bonferroni', 'display', 'off');

figure
create_bar_figure(durM(:, 2), durM(:, 1), durC);
sig_values(durP(2), durP(1));
ylabel('SWR Duration (s)')
set(gca, 'ytick', [0, 0.15, 0.3])
ylim([0, 0.4])

%% IRI
[iriP, iriT, IRI_stats] = anovan(cleanIRI, {IRI_treat, IRI_age,}, 'model', 'interaction', 'display', 'off');
[iriC, iriM, ~, iriNames] = multcompare(IRI_stats, 'Dimension', [1, 2], 'CType', 'bonferroni', 'display', 'off');
figure
create_bar_figure(iriM(:, 2), iriM(:, 1), iriC);
sig_values(iriP(2), iriP(1));
ylabel('Inter-ripple interval (s)')
%set(gca, 'ytick', [0, 3000, 6000])
%ylim([0, 8500])

%% Gamma power

%ctx_Vals = [Gamma.DB2_Ctx Gamma.DBDB2_Ctx Gamma.DB4_Ctx Gamma.DBDB4_Ctx]';
pyr_Vals = [Gamma.DB2_Pyr, Gamma.DBDB2_Pyr, Gamma.DB4_Pyr, Gamma.DBDB4_Pyr]';
%slm_Vals = [Gamma.DB2_SLM Gamma.DBDB2_SLM Gamma.DB4_SLM Gamma.DBDB4_SLM]';
% disp('Ctx Gamma')
% [ctx_P,ctx_Table,ctx_Stats] = anovan(ctx_Vals,{treat_Labs age_Labs},'model','interaction','display','off');
% [ctx_Comparisons,ctx_Means,~,ctx_Names] = multcompare(ctx_Stats,'Dimension',[1 2],'CType','bonferroni','display','off');

disp('Pyr Gamma')
[pyr_P, pyr_Table, pyr_Stats] = anovan(pyr_Vals, {treat_Labs, age_Labs}, 'model', 'interaction', 'display', 'off');
[pyr_Comparions, pyr_Means, ~, pyr_Names] = multcompare(pyr_Stats, 'Dimension', [1, 2], 'CType', 'bonferroni', 'display', 'off');
figure

% disp('Slm Gamma')
% [slmP,slm_Table,slm_Stats] = anovan(slm_Vals,{treat_Labs age_Labs},'model','interaction','display','off');
% [slmC,slmM,~,slmN] = multcompare(slm_Stats,'Dimension',[1 2],'CType','bonferroni','display','off');
%

%% Basic features of SPWRs
both = 0.1;
toph = 0.55;
h = 0.32;
lw = 0.1;
mw = 0.55;
w = 0.35;
riph = 0.174;
mod = 0.175;

control_200 = [0, 0, 204] ./ 255;
control_400 = [0, 0, 204] ./ 255;
DB_200 = [204, 0, 0] ./ 255;
DB_400 = [204, 0, 0] ./ 255;


fig_lim = 5e-6;
% [][] |    |
% [][] |____|
%|    ||    |
%|____||____|
figure

set(gcf, 'color', 'w', 'Position', [100, 100, 1200, 700])
% top left ripple
subplot('Position', [lw, toph + mod, riph, riph])
plot(mean(rip_wav.DB2, 2), 'color', control_200)
% title('DB2')
axis off
ylim([-fig_lim, fig_lim])
vline(625, 'k:')

%bottom left ripple
subplot('Position', [lw + mod, toph + mod, riph, riph])
plot(mean(rip_wav.DBDB2, 2), 'color', DB_200)
ylim([-fig_lim, fig_lim])

%title('DBDB2')
vline(625, 'k:')
axis off


% top right ripple
subplot('Position', [lw, toph, riph, riph])
plot(mean(rip_wav.DB4, 2), 'color', control_400)
ylim([-fig_lim, fig_lim])
axis off
vline(625, 'k:')


%bottom right ripple
subplot('Position', [lw + mod, toph, riph, riph])
plot(mean(rip_wav.DBDB4, 2), 'color', DB_400)
ylim([-fig_lim, fig_lim])
%title('DBDB4')
vline(625, 'k:')
axis off
% top right quad
%Gamma power
subplot('Position', [mw, toph, w, h])
create_bar_figure(pyr_Means(:, 2), pyr_Means(:, 1), pyr_Comparions);
%ylabel('SWR Gamma power')
%sig_values(pyr_P(2), pyr_P(1));
%set(gca,'ytick',[0 5e-7 1e-6])
%ylim([0 1.15e-6])

% bottom left quad
% Duration
subplot('Position', [lw, both, w, h])
dur_fig = create_bar_figure(durM(:, 2), durM(:, 1), durC);
%sig_values(durP(2), durP(1));
%ylabel('SWR Duration (s)')
%set(gca,'ytick',[0 0.15 0.3])
%ylim([0 0.4])

%bottom right quad
% IRI
subplot('Position', [mw, both, w, h])
irir_fig = create_bar_figure(iriM(:, 2), iriM(:, 1), iriC);
%sig_values(iriP(2), iriP(1));
%ylabel('Inter-ripple interval (s)')
%set(gca,'ytick',[0 3000 6000])
%ylim([0 8500])

%% CSD %MS Commented this part %ME
%% Define windows of interest
high_chan = 7; %pyramidal channel
low_chan = 11; % Radiatum channel
pre_win = 1:550; % pre indices
win = 650:750; % ripple indices
post_win = 850:1300; % post indices

CSDm.DB2 = mean(CSD.DB2,3);
CSDm.DB4 = mean(CSD.DB4,3);
CSDm.DBDB2 = mean(CSD.DBDB2,3);
CSDm.DBDB4 = mean(CSD.DBDB4,3);

dipole_DB2_pre = calculate_CSD_dipole(CSD.DB2, high_chan, low_chan, pre_win);
dipole_DB4_pre = calculate_CSD_dipole(CSD.DB4, high_chan, low_chan, pre_win);
dipole_DBDB2_pre = calculate_CSD_dipole(CSD.DBDB2, high_chan, low_chan, pre_win);
dipole_DBDB4_pre = calculate_CSD_dipole(CSD.DBDB4, high_chan, low_chan, pre_win);

dipole_vals_pre = [dipole_DB2_pre; dipole_DBDB2_pre; dipole_DB4_pre; dipole_DBDB4_pre];

dipole_DB2   = calculate_CSD_dipole(CSD.DB2,  high_chan, low_chan, win);
dipole_DB4   = calculate_CSD_dipole(CSD.DB4,  high_chan, low_chan, win);
dipole_DBDB2 = calculate_CSD_dipole(CSD.DBDB2,high_chan, low_chan, win);
dipole_DBDB4 = calculate_CSD_dipole(CSD.DBDB4,high_chan, low_chan, win);

dipole_vals = [dipole_DB2; dipole_DBDB2; dipole_DB4; dipole_DBDB4];

dipole_DB2_post   = calculate_CSD_dipole(CSD.DB2,  high_chan, low_chan, post_win);
dipole_DB4_post   = calculate_CSD_dipole(CSD.DB4,  high_chan, low_chan, post_win);
dipole_DBDB2_post = calculate_CSD_dipole(CSD.DBDB2,high_chan, low_chan, post_win);
dipole_DBDB4_post = calculate_CSD_dipole(CSD.DBDB4,high_chan, low_chan, post_win);

dipole_vals_post = [dipole_DB2_post; dipole_DBDB2_post; dipole_DB4_post; dipole_DBDB4_post];
%MS
if strcmp(user, 'S')
    Group_Ripple_Pre_Ns = [length(dipole_DB2_pre); length(dipole_DBDB2_pre); length(dipole_DB4_pre); length(dipole_DBDB4_pre)];
    
    Group_Ripple_Dur_Ns = [length(dipole_DB2); length(dipole_DBDB2); length(dipole_DB4); length(dipole_DBDB4)];
    
    Group_Ripple_Post_Ns = [length(dipole_DB2_post); length(dipole_DBDB2_post); length(dipole_DB4_post); length(dipole_DBDB4_post)];
end
%ME

%% Make labels
label.DB2age = cell(size(CSD.DB2, 3), 1);
label.DB2treat = cell(size(CSD.DB2, 3), 1);
label.DB2age(:) = {'200'};
label.DB2treat(:) = {'Control'};

label.DB4age = cell(size(CSD.DB4, 3), 1);
label.DB4treat = cell(size(CSD.DB4, 3), 1);
label.DB4age(:) = {'400'};
label.DB4treat(:) = {'Control'};

label.DBDB2age = cell(size(CSD.DBDB2, 3), 1);
label.DBDB2treat = cell(size(CSD.DBDB2, 3), 1);
label.DBDB2age(:) = {'200'};
label.DBDB2treat(:) = {'DBDB'};

label.DBDB4age = cell(size(CSD.DBDB4, 3), 1);
label.DBDB4treat = cell(size(CSD.DBDB4, 3), 1);
label.DBDB4age(:) = {'400'};
label.DBDB4treat(:) = {'DBDB'};

ageLabs = [label.DB2age; label.DBDB2age; label.DB4age; label.DBDB4age];
treatLabs = [label.DB2treat; label.DBDB2treat; label.DB4treat; label.DBDB4treat];

%% Stats

[csdP, csdTable, CSD_stats] = anovan(dipole_vals, {treatLabs, ageLabs}, 'model', 'interaction','display','off');
[csdC, csdM, ~, csdNames] = multcompare(CSD_stats, 'Dimension', [1, 2], 'CType', 'bonferroni','display','off');

[csdP_pre, csdTable_pre, CSD_stats_pre] = anovan(dipole_vals_pre, {treatLabs, ageLabs}, 'model', 'interaction','display','off');
[csdC_pre, csdM_pre, ~, csdNames_pre] = multcompare(CSD_stats_pre, 'Dimension', [1, 2], 'CType', 'bonferroni','display','off');

[csdP_post, csdTable_post, CSD_stats_post] = anovan(dipole_vals_post, {treatLabs, ageLabs}, 'model', 'interaction','display','off');
[csdC_post, csdM_post, ~, csdNames_post] = multcompare(CSD_stats_post, 'Dimension', [1, 2], 'CType', 'bonferroni','display','off');

%% Plot
figure
set(gcf, 'Color', 'w', 'Position', [100, 100, 1200, 600])
xstart = 0.07;
x2start = 0.37;
toph = 0.55;
h = 0.35;
both = 0.15;
w = 0.28;
ylims = [0, 8];
clim = 4.5;
for i = 1:4
    switch i
        case 1
            subplot('Position', [xstart, toph, w, h])
            hand = pcolor(flipud(CSDm.DB2(:, 2:end-1)'));
            title('200 d', 'FontSize', 14)
            ylabel('db/+', 'FontSize', 14, 'FontWeight', 'bold')
        case 2
            subplot('Position',[x2start toph w h])
            hand = pcolor(flipud(CSDm.DB4(:,2:end-1)'));
            title('400 d','FontSize',14)%M needs alteration
            
        case 3
            subplot('Position', [xstart, both, w, h])
            hand = pcolor(flipud(CSDm.DBDB2(:, 2:end-1)'));
            ylabel('db/db', 'FontSize', 14, 'FontWeight', 'bold')
        case 4
            subplot('Position', [x2start, both, w, h])
            hand = pcolor(flipud(CSDm.DBDB4(:, 2:end-1)'));
    end
    set(hand, 'EdgeColor', 'none'), colormap(flipud(hotcold)), shading interp
    rectangle('Position', [pre_win(1), 12 - high_chan(end), pre_win(end) - pre_win(1), 1])
    rectangle('Position', [pre_win(1), 12 - low_chan(end), pre_win(end) - pre_win(1), 1])
    
    rectangle('Position', [win(1), 12 - high_chan(end), win(end) - win(1), 1])
    rectangle('Position', [win(1), 12 - low_chan(end), win(end) - win(1), 1])
    
    rectangle('Position', [post_win(1), 12 - high_chan(end), post_win(end) - post_win(1), 1])
    rectangle('Position', [post_win(1), 12 - low_chan(end), post_win(end) - post_win(1), 1])
    
    hline(6, 'k', 'Pyramidal')
    hline(2, 'k', 'Radiatum')
    
    set(gca, 'xtick', [])
    caxis([-clim, clim])
end

subplot('Position',[0.75 both+0.6 0.2 0.18])
[csdBar] = UCSF_graph([csdM_pre(1:2,2),csdM_pre(3:4,2)]',[csdM_pre(1:2,1),csdM_pre(3:4,1)]',csdC_pre);


title('Pre-Ripple')
ylim(ylims)
set(gca, 'xtick', [])

subplot('Position',[0.75 both+0.3 0.2 0.18])
[csdBar] = UCSF_graph([csdM(1:2,2),csdM(3:4,2)]',[csdM(1:2,1),csdM(3:4,1)]',csdC);


title('Ripple')
ylim(ylims)
set(gca, 'xtick', [])

subplot('Position',[0.75 both 0.2 0.18])
[csdBar] = UCSF_graph([csdM_post(1:2,2),csdM_post(3:4,2)]',[csdM_post(1:2,1),csdM_post(3:4,1)]',csdC_post);


title('Post-Ripple')
ylim(ylims)
% l = legend('db/+','db/db');
% legend('boxoff')
% legend('Location','northoutside')
% age_sig = sig_check(csdP(2));
% db_sig = sig_check(csdP(1));

A = suplabel('Time (s)', 'x', [0.1, 0.1, 0.52, 0.5]);
set(A, 'FontSize', 12, 'FontWeight', 'bold')

%MS
if strcmp(user, 'S')
    T_Pre_Ripple = csdM_pre';
    T_Pre_Ripple = [T_Pre_Ripple;Group_Ripple_Pre_Ns'];
    Datetime_Pre_Ripple = string(datetime('now'));
    cd ('C:\COM\ePhy\dbdb\Data\Outputs\Data\CSD_Notebook')
    Filename_Pre_Ripple = sprintf('Pre_Ripple_CSD_%s.xlsx', Datetime_Pre_Ripple);
    Filename_Pre_Ripple = regexprep(Filename_Pre_Ripple, ' ', '_');
    Filename_Pre_Ripple = regexprep(Filename_Pre_Ripple, ':', '_');
    xlswrite(Filename_Pre_Ripple,T_Pre_Ripple);
    
    T_Ripple = csdM';
    T_Ripple = [T_Ripple;Group_Ripple_Dur_Ns'];
    Datetime_Ripple = string(datetime('now'));
    cd ('C:\COM\ePhy\dbdb\Data\Outputs\Data\CSD_Notebook')
    Filename_Ripple = sprintf('Ripple_CSD_%s.xlsx', Datetime_Ripple);
    Filename_Ripple = regexprep(Filename_Ripple, ' ', '_');
    Filename_Ripple = regexprep(Filename_Ripple, ':', '_');
    xlswrite(Filename_Ripple,T_Ripple);
    
    T_Post_Ripple = csdM_post';
    T_Post_Ripple = [T_Post_Ripple;Group_Ripple_Post_Ns'];
    Datetime_Post_Ripple = string(datetime('now'));
    cd ('C:\COM\ePhy\dbdb\Data\Outputs\Data\CSD_Notebook')
    Filename_Post_Ripple = sprintf('Post_Ripple_CSD_%s.xlsx', Datetime_Post_Ripple);
    Filename_Post_Ripple = regexprep(Filename_Post_Ripple, ' ', '_');
    Filename_Post_Ripple = regexprep(Filename_Post_Ripple, ':', '_');
    xlswrite(Filename_Post_Ripple,T_Post_Ripple);
    
    Datetime_NotyebookCSD = string(datetime('now'));
    Filename_NotyebookCSD = sprintf('NotyebookCSD_Figure_%s.tiff', Datetime_NotyebookCSD);
    Filename_NotyebookCSD = regexprep(Filename_NotyebookCSD, ' ', '_');
    Filename_NotyebookCSD = regexprep(Filename_NotyebookCSD, ':', '_');
    saveas(gcf, Filename_NotyebookCSD);
end
%ME

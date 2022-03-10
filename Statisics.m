cd('C:\Users\ipzach\Documents\dbdb electrophy\Diabetes-Saved-Files')
load('LFP Measures')
hotcold = redblue();
%% Make labels for stats
% We need two labels for a 2-way anova test, one for each variable we are
% testing
label.DB2_age = cell(size(Spec.DB2_Pyr,1),1);
label.DB2_treat = cell(size(Spec.DB2_Pyr,1),1);
label.DB2_age(:) = {'200'};
label.DB2_treat(:) = {'Control'};

label.DB4_age = cell(size(Spec.DB4_Pyr,1),1);
label.DB4_treat = cell(size(Spec.DB4_Pyr,1),1);
label.DB4_age(:) = {'400'};
label.DB4_treat(:) = {'Control'};

label.DBDB2_age = cell(size(Spec.DBDB2_Pyr,1),1);
label.DBDB2_treat = cell(size(Spec.DBDB2_Pyr,1),1);
label.DBDB2_age(:) = {'200'};
label.DBDB2_treat(:) = {'DBDB'};

label.DBDB4_age = cell(size(Spec.DBDB4_Pyr,1),1);
label.DBDB4_treat = cell(size(Spec.DBDB4_Pyr,1),1);
label.DBDB4_age(:) = {'400'};
label.DBDB4_treat(:) = {'DBDB'};

% extra set of labels for groups without removed values for windowing
label.r_DB2_age = cell(length(rip.DB2),1);
label.r_DB2_treat = cell(length(rip.DB2),1);
label.r_DB2_age(:) = {'200'};
label.r_DB2_treat(:) = {'Control'};

label.r_DB4_age = cell(length(rip.DB4),1);
label.r_DB4_treat = cell(length(rip.DB4),1);
label.r_DB4_age(:) = {'400'};
label.r_DB4_treat(:) = {'Control'};

label.r_DBDB2_age = cell(length(rip.DBDB2),1);
label.r_DBDB2_treat = cell(length(rip.DBDB2),1);
label.r_DBDB2_age(:) = {'200'};
label.r_DBDB2_treat(:) = {'DBDB'};

label.r_DBDB4_age = cell(length(rip.DBDB4),1);
label.r_DBDB4_treat = cell(length(rip.DBDB4),1);
label.r_DBDB4_age(:) = {'400'};
label.r_DBDB4_treat(:) = {'DBDB'};

age_Labs = [label.DB2_age; label.DBDB2_age; label.DB4_age; label.DBDB4_age];
treat_Labs = [label.DB2_treat; label.DBDB2_treat; label.DB4_treat; label.DBDB4_treat];

r_age_Labs = [label.r_DB2_age; label.r_DBDB2_age; label.r_DB4_age; label.r_DBDB4_age];
r_treat_Labs = [label.r_DB2_treat; label.r_DBDB2_treat; label.r_DB4_treat; label.r_DBDB4_treat];


% Process CSD into single value
CSD.DB2_max = squeeze(max(CSD.DB2,[],[1 2])); %squeeze(CSD.DB2(625:1500,4,:));
CSD.DB4_max = squeeze(max(CSD.DB4,[],[1 2])); %squeeze(CSD.DB4(625:1500,4,:));
CSD.DBDB2_max = squeeze(max(CSD.DBDB2,[],[1 2])); %squeeze(CSD.DBDB2(625:1500,3,:));
CSD.DBDB4_max = squeeze(max(CSD.DBDB4,[],[1 2])); %squeeze(CSD.DBDB4(625:1500,4,:));

CSD.DB2_min = squeeze(min(min(CSD.DB2,[],1),[],2)); %squeeze(CSD.DB2(625:1500,3,:));
CSD.DB4_min = squeeze(min(min(CSD.DB4,[],1),[],2)); %squeeze(CSD.DB4(625:1500,3,:));
CSD.DBDB2_min = squeeze(min(min(CSD.DBDB2,[],1),[],2)); %squeeze(CSD.DBDB2(625:1500,2,:));
CSD.DBDB4_min = squeeze(min(min(CSD.DBDB4,[],1),[],2)); %squeeze(CSD.DBDB4(625:1500,3,:));

CSD.DB2_full_amp = CSD.DB2_max - CSD.DB2_min;
CSD.DB4_full_amp = CSD.DB4_max - CSD.DB4_min;
CSD.DBDB2_full_amp = CSD.DBDB2_max - CSD.DBDB2_min;
CSD.DBDB4_full_amp = CSD.DBDB4_max - CSD.DBDB4_min;

% Try a more specific window
CSD.DB2_rip = mean(squeeze(CSD.DB2(650:1400,4,:)),1)';
CSD.DB4_rip = mean(squeeze(CSD.DB4(650:1400,4,:)),1)';
CSD.DBDB2_rip = mean(squeeze(CSD.DBDB2(650:1400,4,:)),1)';
CSD.DBDB4_rip = mean(squeeze(CSD.DBDB4(650:1400,4,:)),1)';

CSD.DB2_wav = mean(squeeze(CSD.DB2(650:1400,2,:)),1)';
CSD.DB4_wav = mean(squeeze(CSD.DB4(650:1400,2,:)),1)';
CSD.DBDB2_wav = mean(squeeze(CSD.DBDB2(650:1400,2,:)),1)';
CSD.DBDB4_wav = mean(squeeze(CSD.DBDB4(650:1400,2,:)),1)';

CSD.DB2_amp = CSD.DB2_rip - CSD.DB2_wav;
CSD.DB4_amp = CSD.DB4_rip - CSD.DB4_wav;
CSD.DBDB2_amp = CSD.DBDB2_rip - CSD.DBDB2_wav;
CSD.DBDB4_amp = CSD.DBDB4_rip - CSD.DBDB4_wav;
%% Stats
CSD_vals = [CSD.DB2_amp; CSD.DBDB2_amp; CSD.DB4_amp; CSD.DBDB4_amp];
CSD_full_vals = [CSD.DB2_full_amp; CSD.DBDB2_full_amp; CSD.DB4_full_amp; CSD.DBDB4_full_amp];

[csdP,csdTable,CSD_stats] = anovan(CSD_vals,{treat_Labs age_Labs},'model','interaction');
[csdC,csdM,~,csdNames] = multcompare(CSD_stats,'Dimension',[1 2],'CType','bonferroni');


[csdfP,csdfTable,CSDf_stats] = anovan(CSD_full_vals,{treat_Labs age_Labs},'model','interaction');
[csdfC,csdfM,~,csdfNames] = multcompare(CSD_full_stats,'Dimension',[1 2],'CType','bonferroni');
% Spectrogram
ctxVals = [Spec.DB2_Ctx; Spec.DBDB2_Ctx;Spec.DB4_Ctx; Spec.DBDB4_Ctx];
pyrVals = [Spec.DB2_Pyr; Spec.DBDB2_Pyr;Spec.DB4_Pyr; Spec.DBDB4_Pyr];
slmVals = [Spec.DB2_SLM; Spec.DBDB2_SLM;Spec.DB4_SLM; Spec.DBDB4_SLM];

[ctxP,ctxT,ctxStats] = anovan(ctxVals,{treat_Labs age_Labs},'model','interaction');
[ctxC,ctxM,~,ctxN] = multcompare(ctxStats,'Dimension',[1 2],'CType','bonferroni');

[pyrP,pyrT,pyrStats] = anovan(pyrVals,{treat_Labs age_Labs},'model','interaction');
[pyrC,pyrM,~,pyrN] = multcompare(pyrStats,'Dimension',[1 2],'CType','bonferroni');

[slmP,slmT,slmStats] = anovan(slmVals,{treat_Labs age_Labs},'model','interaction');
[slmC,slmM,~,slmN] = multcompare(slmStats,'Dimension',[1 2],'CType','bonferroni');

% Power
power_vals = [rip.DB2(:,7); rip.DB4(:,7); rip.DBDB2(:,7); rip.DBDB4(:,7)];

[powerP,powerT,power_stats] = anovan(power_vals, {r_treat_Labs, r_age_Labs},'model','interaction');% ,'display','off');
[powerC,powerM,~,powerNames] = multcompare(power_stats,'Dimension',[1 2],'CType','bonferroni');

% Duration
dur_vals = [rip.DB2(:,2)-rip.DB2(:,1);  rip.DBDB2(:,2)-rip.DBDB2(:,1); rip.DB4(:,2)-rip.DB4(:,1); rip.DBDB4(:,2)-rip.DBDB4(:,1)] ./1250;

[durP,durT,dur_stats] = anovan(dur_vals, {r_treat_Labs, r_age_Labs},'model','interaction'); % ,'display','off');
[durC,durM,~,durNames] = multcompare(dur_stats,'Dimension',[1 2],'CType','bonferroni');







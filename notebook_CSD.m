%addpath('C:\Users\ipzach\Documents\MATLAB\Toolbox Zach',...
 %   'C:\Users\ipzach\Documents\MATLAB\spectral-analysis-tools')
 addpath('C:\COM\ePhy\dbdb\code\utils-toolbox\utils-toolbox')
  addpath('C:\COM\ePhy\dbdb\code\spectral-analysis-tools')

%% SWR CSD Analysis
% This script will load LFP and grab the LFP 500ms before and 1500ms after
% the onset of a ripple (max duration recorded 1237ms) then perform CSD on
% that file, save the file to a 3D matrix, and create an average CSD

load('PyramChans.mat')
%cd('C:\Users\ipzach\Documents\MATLAB\Data\dbdb electrophy');
%cd('C:\COM\ePhy\dbdb\code\Diabetes-LFP-Analysis')
cd('C:\COM\ePhy\dbdb\Data\dbdb electrophy');

animals = dir;
hotcold = redblue();
voltConv = 0.000000091555527603759401;
Fs = 1250;
for i = 1:4
    if i ==1
        grouping = 3:9; % DB+ 200D
    elseif i ==2
        grouping = 10:14; % DB+ 400D
    elseif i ==3
        grouping = [15:18 20 21]; % DBDB 200D
    elseif i ==4
        grouping = [22 24:27]; % DBDB 400D
    end
    SwrFullCsd = [];
    for j = grouping
        cd(animals(j).name)
        load('SWR_Index.mat');
        SWR_files = dir('SWR_R_*');
        SWR_files = {SWR_files.name};
        
%         SPWR_files = dir('SPWR_R_*'); % The new analysis with different
%         results that was not considered
%         SPWR_files = {SPWR_files.name};
        
        LFP_files = dir('LFP*');
        LFP_files = {LFP_files.name};
        
        for k = 1:size(SWRLTDIdx,2)
            if ~isempty(SWRLTDIdx(k).R) % makes sure ripple occured during this period
                load(char(SWR_files(k)));
                
                load(char(LFP_files(k)));
                LFP = LFPs{1,2} .*voltConv;
                
                LTDevents = clean_events(...
                    SWRevents(:,1),...
                    length(LFP),...
                    Fs/2,...
                    Fs);
                
                SwrFullCsd = calculate_CSD(...
                    SwrFullCsd,...
                    LFP,...
                    LTDevents,...
                    chans(j),...
                    Fs/2,...
                    Fs,...
                    Fs);
            end % for k s
            
        end % j grouping
        cd ..
        
    end % j
    if i ==1
        full_csd.DB2 = SwrFullCsd;
    elseif i ==2
        full_csd.DB4 = SwrFullCsd;
    elseif i ==3
        full_csd.DBDB2 = SwrFullCsd;
    elseif i ==4
        full_csd.DBDB4 = SwrFullCsd;
    end
end

CSDm.DB2 = mean(full_csd.DB2,3);
CSDm.DB4 = mean(full_csd.DB4,3);

CSDm.DBDB2 = mean(full_csd.DBDB2,3);
CSDm.DBDB4 = mean(full_csd.DBDB4,3);

set(gca,'xtick',[])
caxis([-5 5])

%% stats
% Pick a window
% [x_DB2, y_DB2] = define_window(CSDm.DB2);
% [x_DB4, y_DB4] = define_window(CSDm.DB4);
% [x_DBDB2, y_DBDB2] = define_window(CSDm.DBDB2);
% [x_DBDB4, y_DBDB4] = define_window(CSDm.DBDB4);

%% Define windows of interest
high_chan = 7; %pyramidal channel
low_chan = 11; % Radiatum channel
pre_win = 1:550; % pre indices
win = 650:750; % ripple indices 
post_win = 850:1300; % post indices 
%% Check window
% 
% test = CSDm.DB4;
% test(625:650,11)= 20;
% 
% h = pcolor(flipud(test(:,2:end-1)'));
% set(h,'EdgeColor','none'), colormap(flipud(hotcold))
% title('200 d','FontSize',14)
% ylabel('db/+','FontSize',14,'FontWeight','bold')
% set(gca,'xtick',[])
% caxis([-5 5])
% 
% 
% h = pcolor(test(:,2:end));
% set(h,'EdgeColor','none'), colormap(flipud(hotcold))
% caxis([-5 5])
% rectangle('Position',[high_chan-1 win(i) 1 win(end)-win(1)])
% rectangle('Position',[low_chan-1 win(i) 1 win(end)-win(1)])
% 
% %%
% figure
% h = pcolor(test(:,2:end-1));
% set(h,'EdgeColor','none'), colormap(flipud(hotcold))
% caxis([-5 5])
% rectangle('Position',[high_chan(1) win(1) 2 win(end)-win(1)])
%% Calculate
dipole_DB2_pre   = calculate_CSD_dipole(full_csd.DB2,  high_chan, low_chan, pre_win);
dipole_DB4_pre   = calculate_CSD_dipole(full_csd.DB4,  high_chan, low_chan, pre_win);
dipole_DBDB2_pre = calculate_CSD_dipole(full_csd.DBDB2,high_chan, low_chan, pre_win);
dipole_DBDB4_pre = calculate_CSD_dipole(full_csd.DBDB4,high_chan, low_chan, pre_win);

dipole_vals_pre = [dipole_DB2_pre; dipole_DBDB2_pre; dipole_DB4_pre; dipole_DBDB4_pre];
%MS
Group_Ripple_Pre_Ns = [length(dipole_DB2_pre); length(dipole_DBDB2_pre); length(dipole_DB4_pre); length(dipole_DBDB4_pre)];
%ME
dipole_DB2   = calculate_CSD_dipole(full_csd.DB2,  high_chan, low_chan, win);
dipole_DB4   = calculate_CSD_dipole(full_csd.DB4,  high_chan, low_chan, win);
dipole_DBDB2 = calculate_CSD_dipole(full_csd.DBDB2,high_chan, low_chan, win);
dipole_DBDB4 = calculate_CSD_dipole(full_csd.DBDB4,high_chan, low_chan, win);

dipole_vals = [dipole_DB2; dipole_DBDB2; dipole_DB4; dipole_DBDB4];
%MS
Group_Ripple_Dur_Ns = [length(dipole_DB2); length(dipole_DBDB2); length(dipole_DB4); length(dipole_DBDB4)];
%ME
dipole_DB2_post   = calculate_CSD_dipole(full_csd.DB2,  high_chan, low_chan, post_win);
dipole_DB4_post   = calculate_CSD_dipole(full_csd.DB4,  high_chan, low_chan, post_win);
dipole_DBDB2_post = calculate_CSD_dipole(full_csd.DBDB2,high_chan, low_chan, post_win);
dipole_DBDB4_post = calculate_CSD_dipole(full_csd.DBDB4,high_chan, low_chan, post_win);

dipole_vals_post = [dipole_DB2_post; dipole_DBDB2_post; dipole_DB4_post; dipole_DBDB4_post];
%MS
Group_Ripple_Post_Ns = [length(dipole_DB2_post); length(dipole_DBDB2_post); length(dipole_DB4_post); length(dipole_DBDB4_post)];
%ME
%% Make labels
label.DB2age = cell(size(full_csd.DB2,3),1);
label.DB2treat = cell(size(full_csd.DB2,3),1);
label.DB2age(:) = {'200'};
label.DB2treat(:) = {'Control'};

label.DB4age = cell(size(full_csd.DB4,3),1);
label.DB4treat = cell(size(full_csd.DB4,3),1);
label.DB4age(:) = {'400'};
label.DB4treat(:) = {'Control'};

label.DBDB2age = cell(size(full_csd.DBDB2,3),1);
label.DBDB2treat = cell(size(full_csd.DBDB2,3),1);
label.DBDB2age(:) = {'200'};
label.DBDB2treat(:) = {'DBDB'};

label.DBDB4age = cell(size(full_csd.DBDB4,3),1);
label.DBDB4treat = cell(size(full_csd.DBDB4,3),1);
label.DBDB4age(:) = {'400'};
label.DBDB4treat(:) = {'DBDB'};

ageLabs = [label.DB2age; label.DBDB2age; label.DB4age; label.DBDB4age];
treatLabs = [label.DB2treat; label.DBDB2treat; label.DB4treat; label.DBDB4treat];
%% Stats

[csdP,csdTable,CSD_stats] = anovan(dipole_vals,{treatLabs ageLabs},'model','interaction');
[csdC,csdM,~,csdNames] = multcompare(CSD_stats,'Dimension',[1 2],'CType','bonferroni');

[csdP_pre,csdTable_pre,CSD_stats_pre] = anovan(dipole_vals_pre,{treatLabs ageLabs},'model','interaction');
[csdC_pre,csdM_pre,~,csdNames_pre] = multcompare(CSD_stats_pre,'Dimension',[1 2],'CType','bonferroni');

[csdP_post,csdTable_post,CSD_stats_post] = anovan(dipole_vals_post,{treatLabs ageLabs},'model','interaction');
[csdC_post,csdM_post,~,csdNames_post] = multcompare(CSD_stats_post,'Dimension',[1 2],'CType','bonferroni');
%% Plot
figure
set(gcf,'Color','w','Position',[100 100 1200 600])
xstart = 0.07;
x2start = 0.37;
toph = 0.55;
h = 0.35;
both = 0.15;
w = 0.28;
ylims = [0 10];
clim = 4.5;
for i = 1:4
    switch i
        case 1
            subplot('Position',[xstart toph w h])
            hand = pcolor(flipud(CSDm.DB2(:,2:end-1)'));
            title('200 d','FontSize',14)
            ylabel('db/+','FontSize',14,'FontWeight','bold')
        case 2
            subplot('Position',[x2start toph w h])
            hand = pcolor(flipud(CSDm.DB4(:,2:end-1)'));
            title('400 d','FontSize',14)%M needs alteration
        case 3
            subplot('Position',[xstart both w h])
            hand = pcolor(flipud(CSDm.DBDB2(:,2:end-1)'));
            ylabel('db/db','FontSize',14,'FontWeight','bold')
        case 4
            subplot('Position',[x2start both w h])
            hand = pcolor(flipud(CSDm.DBDB4(:,2:end-1)'));
    end
    set(hand,'EdgeColor','none'), colormap(flipud(hotcold)), shading interp
    rectangle('Position',[pre_win(1) 12-high_chan(end) pre_win(end)-pre_win(1) 1])
    rectangle('Position',[pre_win(1) 12-low_chan(end) pre_win(end)-pre_win(1) 1])
    
    rectangle('Position',[win(1) 12-high_chan(end) win(end)-win(1) 1])
    rectangle('Position',[win(1) 12-low_chan(end) win(end)-win(1) 1])
    
    rectangle('Position',[post_win(1) 12-high_chan(end) post_win(end)-post_win(1) 1])
    rectangle('Position',[post_win(1) 12-low_chan(end) post_win(end)-post_win(1) 1])
    
    hline(6,'k', 'Pyramidal')
    hline(2, 'k', 'Radiatum')
    
    set(gca,'xtick',[])
    caxis([-clim clim])
end

subplot('Position',[0.75 both+0.6 0.2 0.18])
[csdBar] = UCSF_graph([csdM_pre(1:2,2),csdM_pre(3:4,2)]',[csdM_pre(1:2,1),csdM_pre(3:4,1)]',csdC_pre);
%MS
T_Pre_Ripple = csdM_pre';
T_Pre_Ripple = [T_Pre_Ripple;Group_Ripple_Pre_Ns'];
Datetime_Pre_Ripple = string(datetime('now'));
cd ('C:\COM\ePhy\dbdb\Data\Outputs\Data\CSD_Notebook')
Filename_Pre_Ripple = sprintf('Pre_Ripple_CSD_%s.xlsx', Datetime_Pre_Ripple);
Filename_Pre_Ripple = regexprep(Filename_Pre_Ripple, ' ', '_');
Filename_Pre_Ripple = regexprep(Filename_Pre_Ripple, ':', '_');
xlswrite(Filename_Pre_Ripple,T_Pre_Ripple);
%ME
title('Pre-Ripple')
ylim(ylims)

subplot('Position',[0.75 both+0.3 0.2 0.18])
[csdBar] = UCSF_graph([csdM(1:2,2),csdM(3:4,2)]',[csdM(1:2,1),csdM(3:4,1)]',csdC);
%MS
T_Ripple = csdM'; 
T_Ripple = [T_Ripple;Group_Ripple_Dur_Ns'];
Datetime_Ripple = string(datetime('now'));
cd ('C:\COM\ePhy\dbdb\Data\Outputs\Data\CSD_Notebook')
Filename_Ripple = sprintf('Ripple_CSD_%s.xlsx', Datetime_Ripple);
Filename_Ripple = regexprep(Filename_Ripple, ' ', '_');
Filename_Ripple = regexprep(Filename_Ripple, ':', '_');
xlswrite(Filename_Ripple,T_Ripple);
%ME
title('Ripple')
ylim(ylims)

subplot('Position',[0.75 both 0.2 0.18])
[csdBar] = UCSF_graph([csdM_post(1:2,2),csdM_post(3:4,2)]',[csdM_post(1:2,1),csdM_post(3:4,1)]',csdC_post);
%MS
T_Post_Ripple = csdM_post';
T_Post_Ripple = [T_Post_Ripple;Group_Ripple_Post_Ns'];
Datetime_Post_Ripple = string(datetime('now'));
cd ('C:\COM\ePhy\dbdb\Data\Outputs\Data\CSD_Notebook')
Filename_Post_Ripple = sprintf('Post_Ripple_CSD_%s.xlsx', Datetime_Post_Ripple);
Filename_Post_Ripple = regexprep(Filename_Post_Ripple, ' ', '_');
Filename_Post_Ripple = regexprep(Filename_Post_Ripple, ':', '_');
xlswrite(Filename_Post_Ripple,T_Post_Ripple);
%ME
title('Post-Ripple')
ylim(ylims)
% l = legend('db/+','db/db');
% legend('boxoff')
% legend('Location','northoutside')
% age_sig = sig_check(csdP(2));
% db_sig = sig_check(csdP(1));

A = suplabel('Time (s)','x',[0.1 0.1 0.52 0.5]);
set(A,'FontSize',12,'FontWeight','bold')

%MS
Datetime_NotyebookCSD = string(datetime('now'));
Filename_NotyebookCSD = sprintf('NotyebookCSD_Figure_%s.tiff', Datetime_NotyebookCSD);
Filename_NotyebookCSD = regexprep(Filename_NotyebookCSD, ' ', '_');
Filename_NotyebookCSD = regexprep(Filename_NotyebookCSD, ':', '_');
saveas(gcf, Filename_NotyebookCSD);
%ME
%%
% cd('C:\Users\ipzach\Documents\dbdb electrophy\General_Scripts')
cd ('C:\COM\ePhy\dbdb\Data\dbdb electrophy')
%save('DBDB_SPWR_CSD_stats','CSD_stats','CSD_vals','csdC','csdM','csdTable')
save('DBDB_SPWR_CSD_stats','CSD_stats','csdC','csdM','csdTable'); %M CSD_vals was omitted




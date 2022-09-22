addpath('C:\Users\ipzach\Documents\MATLAB\Toolbox Zach',...
    'C:\Users\ipzach\Documents\MATLAB\spectral-analysis-tools')

%% SWR CSD Analysis
% This script will load LFP and grab the LFP 500ms before and 1500ms after
% the onset of a ripple (max duration recorded 1237ms) then perform CSD on
% that file, save the file to a 3D matrix, and create an average CSD

load('PyramChans.mat')
cd('C:\Users\ipzach\Documents\MATLAB\Data\dbdb electrophy');
animals = dir;
hotcold = redblue();
voltConv = 0.000000091555527603759401;
for i = 2% :4
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
        
        LFP_files = dir('LFP*');
        LFP_files = {LFP_files.name};
        
        for k = 1:size(SWRLTDIdx,2)
            if ~isempty(SWRLTDIdx(k).R) % makes sure ripple occured during this period
                load(char(SWR_files(k)));
                
                load(char(LFP_files(k)));
                LFP = LFPs{1,2} .*voltConv;
                
                LTDevents = clean_events(...
                    SWRevents(SWRLTDIdx(k).R,1),...
                    length(LFP),...
                    Fs/2,...
                    Fs);
                
                
            end % for k SWRLTDIdx
            
        end % j grouping
        cd ..
        
    end % j 
    if i ==1
        CSD.DB2 = SwrFullCsd;
    elseif i ==2
        CSD.DB4 = SwrFullCsd;
    elseif i ==3
        CSD.DBDB2 = SwrFullCsd;
    elseif i ==4
        CSD.DBDB4 = SwrFullCsd;
    end
end

CSD.DB2m = mean(CSD.DB2,3);
CSD.DB4m = mean(CSD.DB4,3);

CSD.DBDB2m = mean(CSD.DBDB2,3);
CSD.DBDB4m = mean(CSD.DBDB4,3);

%% Make labels
label.DB2age = cell(size(CSD.DB2,3),1);
label.DB2treat = cell(size(CSD.DB2,3),1);
label.DB2age(:) = {'200'};
label.DB2treat(:) = {'Control'};

label.DB4age = cell(size(CSD.DB4,3),1);
label.DB4treat = cell(size(CSD.DB4,3),1);
label.DB4age(:) = {'400'};
label.DB4treat(:) = {'Control'};

label.DBDB2age = cell(size(CSD.DBDB2,3),1);
label.DBDB2treat = cell(size(CSD.DBDB2,3),1);
label.DBDB2age(:) = {'200'};
label.DBDB2treat(:) = {'DBDB'};

label.DBDB4age = cell(size(CSD.DBDB4,3),1);
label.DBDB4treat = cell(size(CSD.DBDB4,3),1);
label.DBDB4age(:) = {'400'};
label.DBDB4treat(:) = {'DBDB'};

%% stats
CSD.DB2max = squeeze(max(max(CSD.DB2,[],1),[],2)); %squeeze(CSD.DB2(625:1500,4,:));
CSD.DB4max = squeeze(max(max(CSD.DB4,[],1),[],2)); %squeeze(CSD.DB4(625:1500,4,:));
CSD.DBDB2max = squeeze(max(max(CSD.DBDB2,[],1),[],2)); %squeeze(CSD.DBDB2(625:1500,3,:));
CSD.DBDB4max = squeeze(max(max(CSD.DBDB4,[],1),[],2)); %squeeze(CSD.DBDB4(625:1500,4,:));

CSD.DB2min = squeeze(min(min(CSD.DB2,[],1),[],2)); %squeeze(CSD.DB2(625:1500,3,:));
CSD.DB4min = squeeze(min(min(CSD.DB4,[],1),[],2)); %squeeze(CSD.DB4(625:1500,3,:));
CSD.DBDB2min = squeeze(min(min(CSD.DBDB2,[],1),[],2)); %squeeze(CSD.DBDB2(625:1500,2,:));
CSD.DBDB4min = squeeze(min(min(CSD.DBDB4,[],1),[],2)); %squeeze(CSD.DBDB4(625:1500,3,:));

CSD.DB2i = CSD.DB2max - CSD.DB2min;
CSD.DB4i = CSD.DB4max - CSD.DB4min;
CSD.DBDB2i = CSD.DBDB2max - CSD.DBDB2min;
CSD.DBDB4i = CSD.DBDB4max - CSD.DBDB4min;

CSD_vals = [CSD.DB2i; CSD.DBDB2i; CSD.DB4i; CSD.DBDB4i];
ageLabs = [label.DB2age; label.DBDB2age; label.DB4age; label.DBDB4age];
treatLabs = [label.DB2treat; label.DBDB2treat; label.DB4treat; label.DBDB4treat];

[csdP,csdTable,CSD_stats] = anovan(CSD_vals,{treatLabs ageLabs},'model','interaction');
[csdC,csdM,~,csdNames] = multcompare(CSD_stats,'Dimension',[1 2],'CType','bonferroni');
%% Plot
figure
set(gcf,'Color','w','Position',[100 100 1200 600])
xstart = 0.07;
x2start = 0.37;
toph = 0.55;
h = 0.35;
both = 0.15;
w = 0.28;

subplot('Position',[xstart toph w h])
h1 = pcolor(flipud(CSD.DB2m(:,2:end-1))');
set(h1,'EdgeColor','none'),shading interp, colormap(flipud(hotcold))
title('200 d','FontSize',14)
ylabel('db/+','FontSize',14,'FontWeight','bold')
set(gca,'xtick',[])
caxis([-5 5])

subplot('Position',[x2start toph w h])

h2 = pcolor(flipud(CSD.DB4m(:,2:end-1)'));
set(h2,'EdgeColor','none'),shading interp, colormap(flipud(hotcold))
title('400 d','FontSize',14)
set(gca,'xtick',[], 'ytick',[])
caxis([-5 5])

subplot('Position',[xstart both w h])
h3 = pcolor(flipud(CSD.DBDB2m(:,2:end-1)'));
set(h3,'EdgeColor','none'),shading interp, colormap(flipud(hotcold))
ylabel('db/db','FontSize',14,'FontWeight','bold')
set(gca,'xtick',[0 625 1250 1875],'xticklabels',[-0.5, 0, 0.5,1])
caxis([-5 5])


subplot('Position',[x2start both w h])
h4 = pcolor(flipud(CSD.DBDB4m(:,2:end-1)'));
shading interp, colormap(flipud(hotcold))
set(gca,'xtick',[0 625 1250 1875],'xticklabels',[-0.5, 0, 0.5,1])
caxis([-5 5])
set(h4,'EdgeColor','none');

subplot('Position',[0.75 both 0.2 0.8])
[csdBar] = UCSF_graph([csdM(1:2,2),csdM(3:4,2)]',[csdM(1:2,1),csdM(3:4,1)]',csdC);
ylabel('Average Dipole Amplitude (mV)')
l = legend('db/+','db/db');
legend('boxoff')
legend('Location','northoutside')
age_sig = sig_check(csdP(2));
db_sig = sig_check(csdP(1));
if csdP(2) <= 0.05
B = suplabel(['age effect: ' age_sig],'t',[0.72 0.08 0.01 0.85]);
set(B,'FontSize',12)
end
if csdP(1) <= 0.05
C = suplabel(['db effect: ' db_sig],'t',[0.72 0.08 0.01 0.8]);
set(C,'FontSize',12)
end
A = suplabel('Time (s)','x',[0.1 0.1 0.52 0.5]);
set(A,'FontSize',12,'FontWeight','bold')

%%
% cd('C:\Users\ipzach\Documents\dbdb electrophy\General_Scripts')
% save('DBDB_SPWR_CSD_stats','CSD_stats','CSD_vals','csdC','csdM','csdTable')


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


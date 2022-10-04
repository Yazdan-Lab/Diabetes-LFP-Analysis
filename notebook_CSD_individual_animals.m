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
        disp(num2str(j))
        SwrFullCsd = [];
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
                
                SwrFullCsd = calculate_CSD(...
                    SwrFullCsd,...
                    LFP,...
                    LTDevents,...
                    chans(j),...
                    Fs/2,...
                    Fs,...
                    Fs);
                
            end % if
            
        end % for k
        temp_csd = mean(SwrFullCsd,3);
        figure
        h = pcolor(flipud(temp_csd(:,2:end-1)'));
        set(h,'EdgeColor','none'), colormap(flipud(hotcold))
        hline(6,'k', 'Pyramidal')
        hline(3, 'k', 'Radiatum')
        title(['Group: ' num2str(i) ' Animal: ' num2str(j)],'FontSize',14)
        ylabel('db/+','FontSize',14,'FontWeight','bold')
        set(gca,'xtick',[])
        colorbar
        caxis([-15 15])
        drawnow
        
        
        cd ..
        
    end % j group
    if i ==1
        full_csd.DB2 = SwrFullCsd;
    elseif i ==2
        full_csd.DB4 = SwrFullCsd;
    elseif i ==3
        full_csd.DBDB2 = SwrFullCsd;
    elseif i ==4
        full_csd.DBDB4 = SwrFullCsd;
    end % if
end % i = group

CSDm.DB2 = mean(full_csd.DB2,3);
CSDm.DB4 = mean(full_csd.DB4,3);

CSDm.DBDB2 = mean(full_csd.DBDB2,3);
CSDm.DBDB4 = mean(full_csd.DBDB4,3);

set(gca,'xtick',[])
caxis([-5 5])

%% stats
% Pick a window
[x_DB2, y_DB2] = define_window(CSDm.DB2);
[x_DB4, y_DB4] = define_window(CSDm.DB4);
[x_DBDB2, y_DBDB2] = define_window(CSDm.DBDB2);
[x_DBDB4, y_DBDB4] = define_window(CSDm.DBDB4);

%%
high_chan = 7;
low_chan = 11;
pre_win = 1:550;
win = 650:750;
post_win = 850:1300;
%% Check window

test = CSDm.DB4;
test(625:650,11)= 20;

h = pcolor(flipud(test(:,2:end-1)'));
set(h,'EdgeColor','none'), colormap(flipud(hotcold))
title('200 d','FontSize',14)
ylabel('db/+','FontSize',14,'FontWeight','bold')
set(gca,'xtick',[])
caxis([-5 5])


h = pcolor(test(:,2:end));
set(h,'EdgeColor','none'), colormap(flipud(hotcold))
caxis([-5 5])
rectangle('Position',[high_chan-1 win(i) 1 win(end)-win(1)])
rectangle('Position',[low_chan-1 win(i) 1 win(end)-win(1)])
%%

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
            title('400 d','FontSize',14)
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
    hline(3, 'k', 'Radiatum')
    
    set(gca,'xtick',[])
    caxis([-clim clim]);
end

subplot('Position',[0.75 both+0.6 0.2 0.18])
[csdBar] = UCSF_graph([csdM_pre(1:2,2),csdM_pre(3:4,2)]',[csdM_pre(1:2,1),csdM_pre(3:4,1)]',csdC_pre);
title('Pre-Ripple')
ylim(ylims)

subplot('Position',[0.75 both+0.3 0.2 0.18])
[csdBar] = UCSF_graph([csdM(1:2,2),csdM(3:4,2)]',[csdM(1:2,1),csdM(3:4,1)]',csdC);
title('Ripple')
ylim(ylims)

subplot('Position',[0.75 both 0.2 0.18])
[csdBar] = UCSF_graph([csdM_post(1:2,2),csdM_post(3:4,2)]',[csdM_post(1:2,1),csdM_post(3:4,1)]',csdC_post);
title('Post-Ripple')
ylim(ylims)
% l = legend('db/+','db/db');
% legend('boxoff')
% legend('Location','northoutside')
% age_sig = sig_check(csdP(2));
% db_sig = sig_check(csdP(1));

A = suplabel('Time (s)','x',[0.1 0.1 0.52 0.5]);
set(A,'FontSize',12,'FontWeight','bold')

%%
% cd('C:\Users\ipzach\Documents\dbdb electrophy\General_Scripts')
% save('DBDB_SPWR_CSD_stats','CSD_stats','CSD_vals','csdC','csdM','csdTable')





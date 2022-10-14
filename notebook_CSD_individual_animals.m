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
upper = 27;
wbar = waitbar(0,'Starting...');
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
        bval = j/upper;
        waitbar(bval, wbar, ['Processing animal ' num2str(j) ': ' num2str(round(bval*100)) '%']);
        SwrFullCsd = [];
        cd(animals(j).name)
        load('SWR_Index.mat');
        SPWR_files = dir('SPWR_R_*');
        SPWR_files = {SPWR_files.name};
        
        LFP_files = dir('LFP*');
        LFP_files = {LFP_files.name};
        
        for k = 1:size(SWRLTDIdx,2)
            if ~isempty(SWRLTDIdx(k).R) % makes sure ripple occured during this period
                load(char(SPWR_files(k)));
                
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
                
            end % if
            
        end % for k
        temp_csd = mean(SwrFullCsd,3);
        figure
        h = pcolor(flipud(temp_csd(:,2:end-1)'));
        set(h,'EdgeColor','none'), colormap(flipud(hotcold))
        hline(6,'k', 'Pyramidal')
        hline(3, 'k', 'Radiatum')
        title(['Group: ' num2str(i) ' Animal: ' num2str(j) ' Num Samples: ' ... 
               num2str(size(SwrFullCsd,3))],'FontSize',14)
        ylabel('Electrode','FontSize',14,'FontWeight','bold')
        set(gca,'xtick',[])
        xlabel('Time')
        colorbar
        caxis([-5 5])
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

close(wbar)

CSDm.DB2 = mean(full_csd.DB2,3);
CSDm.DB4 = mean(full_csd.DB4,3);

CSDm.DBDB2 = mean(full_csd.DBDB2,3);
CSDm.DBDB4 = mean(full_csd.DBDB4,3);





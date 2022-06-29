%% Define path and open EEGlab
clc
clear all
close all

repoPath = '/data/brunobian/Documents/Repos/';
lm_Conf = struct();
lm_Conf.path         = [repoPath 'LMM-CBP/'];
lm_Conf.eeglabpath   = [repoPath 'eeglab14_1_1b'];
lm_Conf.ftpath       = [repoPath 'fieldtrip'];
lm_Conf.custonFc     = [repoPath 'coregistration2022-analyses/functions/'];
lm_Conf.datapathBase = [repoPath 'coregistration2022-analyses/data/'];

cd([lm_Conf.datapathBase])
addpath(genpath(lm_Conf.path))
addpath(genpath(lm_Conf.custonFc))
addpath(genpath(lm_Conf.ftpath),'-end'); 
addpath(genpath(lm_Conf.eeglabpath),'-end'); 
%% Load data from subjects and export to CSV
lm_Conf.programRun = 'R'; %'matlab'

fs = {'broadband','Theta','Alpha','Beta'};
for iFolder = 1
lm_Conf.datapath = [lm_Conf.datapathBase  fs{iFolder} '/'];
lm_Conf.csvPath  = [lm_Conf.datapathBase 'csv/' fs{iFolder} '/'];
cd([lm_Conf.datapath])
[lm_Conf, SUJ]= definePath(lm_Conf);

tic
    for iSuj = 1:length(SUJ)
        % Load EEGLab matrix with information from the experiment into DATA
        EEG = pop_loadset([lm_Conf.datapath  SUJ(iSuj).fileName]); 
        
        EEG = pop_resample(EEG,128);
        EEG = pop_select(EEG, 'time', [-.1 .6]);
        EEG = pop_rmbase(EEG, [-100 -50]);
        EEG = eeg_checkset( EEG );
        
        T = generateData(EEG, iSuj);
        if strcmpi(lm_Conf.programRun, 'R')
            lm_Conf.export = 'csv';
        elseif strcmpi(lm_Conf.programRun, 'matlab')
            lm_Conf.export = 'mat';
        end

        lm_Conf.sep = ';';
        lm_Conf = lm_exportErpCsv(T, EEG, iSuj, lm_Conf); 
    end
processData(lm_Conf)
end
toc
fprintf('Done\n')
%% [R] Run LMM with nIter and nCores permutations across both ranEf

% M0_pos
    fixEf   = 'pred:sntType + pos:sntType + freq'; 
% M1_pos
    % fixEf   = 'sntTypePrePost:pred + sntType:pos + freq'; 
% M2_pos
    % fixEf   = 'sntTypePrePost:pred + sntType:pos + freq + fixDur + saccDur'; 

ranEf   = '(1|suj_id)'; 
nIter   = 1000;
modType = 'lmm'; 
nCores  = 22;

lm_Conf.lmmOutPath            = [lm_Conf.datapathBase 'LMM/results/M0/suj_id/'];
% lm_Conf.lmmOutPath            = [lm_Conf.datapathBase 'LMM/results/M1/suj_id/'];
% lm_Conf.lmmOutPath            = [lm_Conf.datapathBase 'LMM/results/M2/suj_id/'];
% lm_Conf.lmmOutPath            = [lm_Conf.datapathBase 'LMM/results/M2_Theta/suj_id/'];
% lm_Conf.lmmOutPath            = [lm_Conf.datapathBase 'LMM/results/M2_Alpha/suj_id/'];
% lm_Conf.lmmOutPath            = [lm_Conf.datapathBase 'LMM/results/M2_Beta/suj_id/'];
lm_Conf.nohupOutPath          = [lm_Conf.datapathBase 'LMM/nohupOutputs/'];
lm_Conf.permutationMatPath    = [lm_Conf.datapathBase 'LMM/permutations/']; 
lm_Conf.csvPath               = [lm_Conf.datapathBase 'csv/broadband/'];
% lm_Conf.csvPath               = [lm_Conf.datapathBase 'csv/Theta/'];
% lm_Conf.csvPath               = [lm_Conf.datapathBase 'csv/Alpha/'];
% lm_Conf.csvPath               = [lm_Conf.datapathBase 'csv/Beta/'];

lm_Conf.customFunsPath        = lm_Conf.custonFc;     
lm_Conf.rFunctionsPath        = [lm_Conf.path '/R_functions/'];
lm_Conf.bashPath              = [lm_Conf.path '/bash_functions/'];
lm_Conf.permutationVariable   = 'suj_id';
lm_Conf.nTimes = 90;
lm_Conf.startTime = 1;
lm_Conf.cfgPath = [lm_Conf.datapathBase 'LMM/cfg.csv']; 

lm_parallelRunLMM(fixEf, ranEf, nIter, modType, nCores, lm_Conf)
%% folders
cd(lm_Conf.datapathBase)
addpath(genpath(lm_Conf.custonFc))

folders = {'M0','M1','M2','M2_Theta','M2_Alpha','M2_Beta'};

alphaTval = 1.86;
load('times')
lm_Conf.nTimes = 90;
load([lm_Conf.custonFc '/coords_LNI_128_toreplaceinEEG'])
colores = colormap;
elect = 1:128;
%% Load data from LMM

for ifolder=4%length(folders)
disp(folders{ifolder})

lm_Conf.nTimes = 90;
permToLoad = {'suj_id'};

lm_Conf.lmmOutPath         = [lm_Conf.datapath 'LMM/results/'...
                              folders{ifolder}...
                              '/'];
                          
lm_Conf.matricesLoadedPath = [lm_Conf.datapath 'LMM/matrices/'...
                              folders{ifolder}...
                              '/'];

lm_Conf.programRun = 'R'; 
values = lm_loadLmmData(permToLoad, lm_Conf);
load coords_LNI_128_toreplaceinEEG

% Custum clustering
lm_Conf.clusteralpha = 0.05;
lm_Conf.minnbtime = 1;
lm_Conf.minnbchan = 2;
lm_Conf.tail = 0;
lm_Conf.alpha = 1.86;

[perms pvals]= lm_clustering(permToLoad, lm_Conf, CHANS);
end
%% Plot custom (with perms)
lm_Conf.nTimes = 90;

for ifolder=1
disp(folders{ifolder})

lm_Conf.matricesLoadedPath = [lm_Conf.datapath 'LMM/matrices/'...
                              folders{ifolder}...
                              '/'];
                          
lm_Conf.lmmOutPath         = [lm_Conf.datapath 'LMM/results/'...
                              folders{ifolder}...
                              '/'];                        

load([lm_Conf.matricesLoadedPath 'clustersNum_suj_id.mat'])
load([lm_Conf.matricesLoadedPath 'pvalsClust_suj_id.mat'])
load([lm_Conf.matricesLoadedPath 'suj_id.mat'])
perms = clusters;

close all
colores = colormap;
fields = fieldnames(values.t);
alphaTval = 1.86;
% t = EEG.times(1:lm_Conf.nTimes);
t = (1:lm_Conf.nTimes);
load('times')
for iv = 1:length(fields) 
    
    v = fields{iv};

    tval   = squeeze(values.t.(v)(:,:,1));
    
    figure(iv);clf;
    
    subplot(1,2,1)
    set(gcf,'Color','w','Position', get(0, 'Screensize'))
    colormap redbluecmap;
        col_ax = [-3 3];
        mask = abs(tval)>alphaTval;
        elect = 1:size(tval,1);
        tmp = tval.*mask;
        imagesc(t, elect, tmp(elect,:), col_ax);
%         colorbar
        hold on
            plot([0 0], [1 elect(end)], 'k--', 'LineWidth',2)
%             set(gca, 'YTickLabel',{})
%             set(gca, 'XTickLabel',{})
            box on
        hold off
        tit=regexprep(v, '_', '-');
        title(tit)
        f=getframe(gca);
        
    p1=perms.(v).pos; u=unique(p1); s=0:length(u); p2=p1;
    for i=1:length(u); p2(p1==u(i)) = s(i);end
    cp = colores(size(colores,1)-length(u)+1:end,:); cp(1,:)=[1,1,1];

    n1=perms.(v).neg; u=unique(n1); s=0:length(u); n2=n1;
    for i=1:length(u); n2(n1==u(i)) = s(i);end
    cn = flip(colores(1:length(u),:),1);

    c = [cn(1:end-1,:);cp];    
    mask_both = p2 - n2;
       
    subplot(1,2,2)
    set(gcf,'Color','w','Position', get(0, 'Screensize'))
        tmp = mask_both;
        valores = unique(mask_both);
        imagesc(t, elect, tmp(elect,:), [-1 1]);
        colormap(gca,c)
        colormap redbluecmap;
        hold on
            plot([0 0], [1 elect(end)], 'k--', 'LineWidth',2)
            set(gcf, 'Position',  [100, 100, 1100, 400])
%             set(gca, 'YTickLabel',{})
%             set(gca, 'XTickLabel',{})
            box on
        hold off
        nperms = size(values.t.Intercept,3);
        title(['Clusters significativos (' num2str(nperms), ' perms)'])
        set(gcf,'Color','w','Position', get(0, 'Screensize'))
        
        saveas(gcf, [lm_Conf.lmmOutPath tit '.png'])
        saveas(gcf, [lm_Conf.lmmOutPath tit '.eps'],'eps2c')
end
end

%% Topoplots

timelim= [0 200 350 500 600];
nT = (length(timelim)-1);

for ifolder=2%1:length(folders)
disp(folders{ifolder})
lm_Conf.matricesLoadedPath = [lm_Conf.datapath 'LMM/matrices/'...
                              folders{ifolder}...
                              '/'];
                          
lm_Conf.lmmOutPath         = [lm_Conf.datapath 'LMM/results/'...
                              folders{ifolder}...
                              '/'];                        

load([lm_Conf.matricesLoadedPath 'clustersNum_suj_id.mat'])
load([lm_Conf.matricesLoadedPath 'pvalsClust_suj_id.mat'])
load([lm_Conf.matricesLoadedPath 'suj_id.mat'])
perms = clusters;

fields = fieldnames(values.t);


for iv = 1:length(fields) 
figure(iv+10);clf;
set(gcf, 'Position', get(0, 'Screensize'));
for it = 1:nT
    v = fields{iv};
    tval = squeeze(values.t.(v)(:,:,1));

   
    lim = [-2 2];
    electrodes = [];
    elec = {electrodes};    

    % All
    subplot(2,nT,it)
    mask = abs(tval)>alphaTval;
    tmp = tval.*mask;
    ventana = tmp(:,t>timelim(it) & t<timelim(it+1));
    meanWindow = nanmean(ventana,2);
    
    topoplot(meanWindow, ...
             CHANS.chanlocs, ...
             'emarker2', elec, ...
             'maplimits',lim , ...
             'colormap', colormap('redbluecmap'));
    set(gcf,'Color','w')
    title([num2str(timelim(it)) ' - ' num2str(timelim(it+1))]) 
    
    % Clusters
    subplot(2,nT, it+nT)
    p1=perms.(v).pos; u=unique(p1); s=0:length(u); p2=p1;
    for i=1:length(u); p2(p1==u(i)) = s(i);end
    cp = colores(size(colores,1)-length(u)+1:end,:); cp(1,:)=[1,1,1];

    n1=perms.(v).neg; u=unique(n1); s=0:length(u); n2=n1;
    for i=1:length(u); n2(n1==u(i)) = s(i);end
    cn = flip(colores(1:length(u),:),1);

    c = [cn(1:end-1,:);cp];    
    mask_both = p2 + n2;  
    tmp = tval.*mask_both;
    ventana = tmp(:,t>timelim(it) & t<timelim(it+1));
    meanWindow = nanmean(ventana,2);

    topoplot(meanWindow, ...
             CHANS.chanlocs, ...
             'emarker2', elec, ...
             'maplimits',lim , ...
             'colormap', colormap('redbluecmap'));
    set(gcf,'Color','w')
    title([num2str(timelim(it)) ' - ' num2str(timelim(it+1))])    
        
end
tit=regexprep(v, '_', '-');

% saveas(gcf, [lm_Conf.lmmOutPath tit '_topp.png'])
% saveas(gcf, [lm_Conf.lmmOutPath tit '_topp.eps'],'eps2c')
end
%     f=getframe(gca);
%     if guardar; saveas(gcf,[figName '_topo_cluster.eps'],'eps2c');end
%     if guardar; saveas(gcf,[figName '_topo_cluster.png'],'png');end
end
%% Convergence failuers analyses
ts = values.t;
fn = fieldnames(ts);

a=0;
for i=1:length(fn)
    c=fn{i};
    a = a + (ts.(c)(:,:,500)==0);
end
unique(a)
mean(a(:)/7)*100
%% Topoplot pos
load coords_LNI_128_toreplaceinEEG

camposM_rel = {'memory_rel_m3' 'memory_rel_m2' 'memory_rel_m1' 'memory_rel_0' 'memory_rel_1' 'memory_rel_2' 'memory_rel_3' 'memory_rel_4' 'memory_rel_5'};
titulos_m_rel = {'RP-3' 'RP-2' 'RP-1' 'RP' 'RP+1' 'RP+2' 'RP+3' 'RP+4' 'RP+5' };
camposM_abs = {'memory_abs_3' 'memory_abs_4' 'memory_abs_5' 'memory_abs_6' 'memory_abs_7' 'memory_abs_8' 'memory_abs_9' 'memory_abs_10' 'memory_abs_11'};
camposC = {'common_3' 'common_4' 'common_5' 'common_6' 'common_7' 'common_8' 'common_9' 'common_10' 'common_11'};
titulos_c = { '3' '4' '5' '6' '7' '8' '9' '10' '11'};

camposM_rel = {'memory_rel_m4' 'memory_rel_0' 'memory_rel_1' 'memory_rel_4' 'memory_rel_7'};
camposM_abs = {'memory_abs_2'  'memory_abs_6'  'memory_abs_7' 'memory_abs_10' 'memory_abs_13'};
camposC     = {'common_2'      'common_6'      'common_7'     'common_10'     'common_13'    };

titulos_m_rel = {'RP-2-3-2-1' 'RP' 'RP' 'RP+1+2+3' 'RP+4+5+6' 'RP+7+8+9+10' };
titulos_c = { '2-3-4-5' '6' '7-8-9' '10-11-12' '13-14-15-16'};

figure(12);clf
for i = 1:length(camposM_rel)
    electrodes = [];
    elec = {electrodes}; 
    t = times > 0 & times <150;
    cM_abs = camposM_abs{i};
    cM_rel = camposM_rel{i};
    cC = camposC{i};
    lim = [-.3 .3];
    
    % common
    subplot(3,cols,i)
    meanWindow = mean(mean(media.(cC)(:,t),2),3);
    topoplot(meanWindow, ...
             CHANS.chanlocs, ...
             'emarker2', elec, ...
             'maplimits',lim , ...
             'colormap', colormap('redbluecmap'));
    title(['common | pal ' titulos_c{i}])     
    

    % memory abs
    subplot(3,cols,i+cols)
    meanWindow = mean(mean(media.(cM_rel)(:,t),2),3);
    topoplot(meanWindow, ...
             CHANS.chanlocs, ...
             'emarker2', elec, ...
             'maplimits',lim , ...
             'colormap', colormap('redbluecmap'));
    title(['memory | ' titulos_c{i}])     

    
    % memory rel
    cols = length(camposM_rel);
    subplot(3,cols,i+2*cols)
    meanWindow = mean(mean(media.(cM_abs)(:,t),2),3);
    topoplot(meanWindow, ...
             CHANS.chanlocs, ...
             'emarker2', elec, ...
             'maplimits',lim , ...
             'colormap', colormap('redbluecmap'));
    title(['memory | ' titulos_m_rel{i}])     

    

end
set(gcf,'Color','w','Position', get(0, 'Screensize'))
% saveas(gcf, [lm_Conf.lmmOutPath '../toposPorPos.png'])



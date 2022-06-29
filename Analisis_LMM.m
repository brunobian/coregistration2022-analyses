%% Define path and open EEGlab
clc
clear all
close all
lm_Conf = struct();
lm_Conf.path = '/data/brunobian/Documents/Repos/LMM-CBP/';
lm_Conf.eeglabpath  = '/data/brunobian/Documents/Repos/eeglab14_1_1b';
lm_Conf.ftpath = '/data/brunobian/Documents/Repos/fieldtrip-20180513';
lm_Conf.ftpath = '/data/brunobian/Documents/Repos/fieldtrip';
lm_Conf.custonFc = '/media/brunobian/DATABRUNO/corregistro/Analysis_bb/2020/my_functions';
lm_Conf.datapath = '/media/brunobian/DATABRUNO/Co-registro_2018/Data/12-edited_forERPs/';

cd([lm_Conf.datapath])
addpath(genpath(lm_Conf.path))
addpath(genpath(lm_Conf.custonFc))
addpath(genpath(lm_Conf.ftpath),'-end'); 
addpath(genpath(lm_Conf.eeglabpath),'-end'); 
%% Load data from subjects and export to CSV
lm_Conf.programRun = 'R'; %'matlab'

fs = {'conICA','AlphaHigh','AlphaLow','Alpha','Beta','Theta'};
for iFolder = 5
lm_Conf.datapath = ['/media/brunobian/DATABRUNO/Co-registro_2018/Data/12-edited_forERPs/' fs{iFolder} '/'];
lm_Conf.csvPath  = [lm_Conf.datapath '../../csv_tesis_center_R_causal/' fs{iFolder} '/'];
cd([lm_Conf.datapath])
[lm_Conf, SUJ]= definePath2020(lm_Conf);

tic
    for iSuj = 1:length(SUJ)
        % Load EEGLab matrix with information from the experiment into DATA
        EEG = pop_loadset([lm_Conf.datapath  SUJ(iSuj).fileName]); 
        EEG = pop_resample(EEG,128);
        
%     indE = [1 3 6 8 10 12 14 16 18 19 21 23	25 27 29 31 ...      A
%             35 37 39 41 42 44 47 49 51 52 54 56 58 60 62 64 ...  B
%             67 69 71 72 74 77 79 81 83 85 87 90 92 94 96 ...     C
%             99 101 103 105 107 109 110 113 115 117 119 121 123 126 128]; 
%     EEG = pop_select(EEG, 'channel', indE);

%     EEG = pop_select(EEG, 'time', [-.1 .6]);
        EEG = pop_rmbase(EEG, [-100 -50]);

        EEG = eeg_checkset( EEG );
        T = generateData(EEG, iSuj);

        if strcmpi(lm_Conf.programRun, 'R')
            lm_Conf.export = 'csv';
        elseif strcmpi(lm_Conf.programRun, 'matlab')
            lm_Conf.export = 'mat';
        end

        lm_Conf.sep = ';';

        keyboard
        lm_Conf = lm_exportErpCsv(T, EEG, iSuj, lm_Conf); 

    end
processData(lm_Conf)
end
toc
fprintf('Done\n')
%% [R] Run LMM with nIter and nCores permutations across both ranEf

fixEf   = 'pred:sntType + pos:sntType + freq'; % M0_pos
% fixEf   = 'sntTypePrePost:pred + sntType:pos + freq'; % M1_pos
% fixEf   = 'sntTypePrePost:pred + sntType:pos + freq + fixDur + saccDur'; % M2_pos

ranEf   = '(1|suj_id)'; 
nIter   = 500;
modType = 'lmm'; %'lmm'
% lm_Conf.fixEf_remove = '(Intercept) sntType-0.5:pos sntType0.5:pos fixDur saccDur';
% lm_Conf.fixEf_remove = '(Intercept) sntTypePrePost0:pred sntTypePrePost1:pred sntTypePrePost2:pred fixDur saccDur';
nCores  = 22;

lm_Conf.datapath = ['/media/brunobian/DATABRUNO/Co-registro_2018/Data/'];
lm_Conf.lmmOutPath            = [lm_Conf.datapath 'LMM/results/paper/M0_resample/suj_id/'];
lm_Conf.nohupOutPath          = [lm_Conf.datapath 'LMM/nohupOutputs/'];
lm_Conf.permutationMatPath    = [lm_Conf.datapath 'LMM/permutations/']; 
lm_Conf.csvPath               = [lm_Conf.datapath 'csv_tesis_center_R_causal/conICA/'];

lm_Conf.customFunsPath        = '/media/brunobian/DATABRUNO/corregistro/Analysis_bb/2020/my_functions';     
lm_Conf.rFunctionsPath        = [lm_Conf.path '/R_functions/'];
lm_Conf.bashPath              = [lm_Conf.path '/bash_functions/'];
lm_Conf.permutationVariable   = 'suj_id';
lm_Conf.nTimes = 90;
lm_Conf.startTime = 1;
lm_Conf.cfgPath = [lm_Conf.datapath 'LMM/cfg.csv']; 


lm_parallelRunLMM(fixEf, ranEf, nIter, modType, nCores, lm_Conf)
%% [matlab] Run LMM with nIter and nCores permutations across both ranEf

lm_Conf.datapath = ['/media/brunobian/DATABRUNO/Co-registro_2018/Data/'];
cfg=[];
cfg.nIter   = 100;
cfg.perVar  = 'across';
cfg.cfgPath = [lm_Conf.datapath 'LMM/cfg-Theta']; 
cfg.nhoPath = [lm_Conf.datapath 'LMM/nohupOutputs']; 
cfg.inPath  = [lm_Conf.datapath 'csv_tesis/Theta/'];
cfg.outPath = [lm_Conf.datapath 'LMM/results/tesis/Theta/across/'];
cfg.perPath = [lm_Conf.datapath 'LMM/permutations/']; 
cfg.fixEf   = 'sntTypePrePost:pred + sntType:pos + freq'; %+ saccDur:sntType + fixDur:sntType  saccDur + pred:sntType + pos:sntType'; % 
cfg.categoricals = {'sntType', 'PostRP', 'fixDurCat','sntTypePrePost'};
cfg.ranEf   = '(1|suj_id)'; 
cfg.modType = 'lmm';
cfg.tStart  = 1;
cfg.tEnd    = 90;
cfg.nTimes  = cfg.tEnd-cfg.tStart;
cfg.mPath = [lm_Conf.path '/m_functions/'];

cfg.nCores = 8;
cfg.timerPerCore = round(cfg.nTimes/cfg.nCores); 

cfg.maxNumCompThreads=1;

lm_completRun(cfg)
lm_parallelRunLMM_m(cfg)

%% folders
lm_Conf.datapath = '/media/brunobian/DATABRUNO/Co-registro_2018/Data/';
cd(lm_Conf.datapath)
lm_Conf.custonFc = '/media/brunobian/DATABRUNO/corregistro/Analysis_bb/2020/';
addpath(genpath(lm_Conf.custonFc))

folders = {'M0_resample','M2_pos_Beta','M0_pos', 'M1_pos', 'M1_rank', 'M2_pos', 'M2_rank', 'M1_pos_Theta'...
'M1_rank_Theta', 'M2_pos_Theta', 'M2_rank_Theta', 'M1_pos_Alpha', 'M1_rank_Alpha'...
'M2_pos_Alpha', 'M2_rank_Alpha', 'M1_pos_Beta', 'M1_rank_Beta', ...
'M2_rank_Beta','M4_Alpha','M4_Theta','M4_Beta' };

alphaTval = 1.86;
load('times')
lm_Conf.nTimes = 90;
load([lm_Conf.custonFc 'my_functions/coords_LNI_128_toreplaceinEEG'])
colores = colormap;
elect = 1:128;
%% Load data from LMM

for ifolder=4%length(folders)
disp(folders{ifolder})
    
addpath(genpath(lm_Conf.ftpath),'-end'); 

lm_Conf.datapath = '/media/brunobian/DATABRUNO/Co-registro_2018/Data/';
cd(lm_Conf.datapath)
% lm_Conf.custonFc = '/home/brunobian/Documents/Repos/corregistro/Analysis_bb/2020/';
addpath(genpath(lm_Conf.custonFc))


lm_Conf.nTimes = 90;
% permToLoad = {lm_Conf.permutationVariable};
permToLoad = {'suj_id'};

lm_Conf.lmmOutPath         = [lm_Conf.datapath 'LMM/results/paper/'...
                              folders{ifolder}...
                              '/'];
                          
lm_Conf.matricesLoadedPath = [lm_Conf.datapath 'LMM/matrices/paper/'...
                              folders{ifolder}...
                              '/'];

lm_Conf.programRun = 'R'; %'matlab'
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
lm_Conf.datapath = '/media/brunobian/DATABRUNO/Co-registro_2018/Data/';
cd(lm_Conf.datapath)
lm_Conf.custonFc = '/home/brunobian/Documents/Repos/corregistro/Analysis_bb/2020/';
addpath(genpath(lm_Conf.custonFc))
lm_Conf.nTimes = 90;

for ifolder=1
disp(folders{ifolder})

lm_Conf.matricesLoadedPath = [lm_Conf.datapath 'LMM/matrices/paper/'...
                              folders{ifolder}...
                              '/'];
                          
lm_Conf.lmmOutPath         = [lm_Conf.datapath 'LMM/results/paper/'...
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

%% AIC Analisis
% Primero definir los AIC
if 0
AIC = [];
for ifolder=1:length(folders)
disp(folders{ifolder})
lm_Conf.matricesLoadedPath = [lm_Conf.datapath 'LMM/matrices/paper/'...
                              folders{ifolder}...
                              '/'];
                          
lm_Conf.lmmOutPath         = [lm_Conf.datapath 'LMM/results/paper/'...
                              folders{ifolder}...
                              '/'];                        

% load([lm_Conf.matricesLoadedPath 'clustersNum_suj_id.mat'])
% load([lm_Conf.matricesLoadedPath 'pvalsClust_suj_id.mat'])
load([lm_Conf.matricesLoadedPath 'suj_id.mat'])

AIC.(folders{ifolder}) = values.AIC.AIC(:,:,1);
end
AIC.M0_pos(:,33) = AIC.M1_pos(:,33);
end

figFolder = '/media/brunobian/DATABRUNO/Co-registro_2018/Data/LMM/results/paper/AICs/';
clim = [-20 20];

figure(1);clf
colormap redbluecmap;
set(gcf,'Color','w','Position', get(0, 'Screensize'))
imagesc(t, elect, (AIC.M1_pos - AIC.M0_pos),clim)
hold on
    plot([0 0], [1 elect(end)], 'k--', 'LineWidth',2)
    box on
hold off
title("AIC: M1\_pos - M0\_pos")
colorbar
tit = 'AIC M1 - AIC M0';
saveas(gcf, [figFolder tit '.png'])
saveas(gcf, [figFolder tit '.eps'],'eps2c')


figure(2);clf
colormap redbluecmap;
set(gcf,'Color','w','Position', get(0, 'Screensize'))
imagesc(t, elect, (AIC.M2_pos - AIC.M1_pos),clim)% - M2_pos_AIC)
hold on
    plot([0 0], [1 elect(end)], 'k--', 'LineWidth',2)
    box on
hold off
title("AIC: M2\_pos - M1\_pos")
colorbar
tit = 'AIC M2 - AIC M1';
saveas(gcf, [figFolder tit '.png'])
saveas(gcf, [figFolder tit '.eps'],'eps2c')


figure(3);clf
colormap redbluecmap;
set(gcf,'Color','w','Position', get(0, 'Screensize'))
imagesc(t, elect, (AIC.M2_rank - AIC.M2_pos),clim)% - M2_pos_AIC)
hold on
    plot([0 0], [1 elect(end)], 'k--', 'LineWidth',2)
    box on
hold off
title("AIC: M2\_rank - M2\_pos")
colorbar
tit = 'AIC M2_rank - AIC M2_pos';
saveas(gcf, [figFolder tit '.png'])
saveas(gcf, [figFolder tit '.eps'],'eps2c')
%% Topoplots
lm_Conf.datapath = '/media/brunobian/DATABRUNO/Co-registro_2018/Data/';
cd(lm_Conf.datapath)
lm_Conf.custonFc = '/home/brunobian/Documents/Repos/corregistro/Analysis_bb/2020/';
addpath(genpath(lm_Conf.custonFc))

timelim= [0 200 350 500 600];
% timelim= [-50 0 150 300 450 600];
nT = (length(timelim)-1);

for ifolder=2%1:length(folders)
disp(folders{ifolder})
lm_Conf.matricesLoadedPath = [lm_Conf.datapath 'LMM/matrices/paper/'...
                              folders{ifolder}...
                              '/'];
                          
lm_Conf.lmmOutPath         = [lm_Conf.datapath 'LMM/results/paper/'...
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

saveas(gcf, [lm_Conf.lmmOutPath tit '_topp.png'])
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

%% Load remef type-pred
lm_Conf.datapath = '/media/brunobian/DATABRUNO/Co-registro_2018/Data/';
cd(lm_Conf.datapath)

% cd('csv_tesis_center_R_causal/conICA/') % Original
% cd('LMM/results/paper/M1_remef/suj_id/') % M1
cd('LMM/results/paper/M2_remef_Beta_(keep_pred)/suj_id/') % M2

files = dir("./");
ind = cellfun(@(y) ~isempty(y), arrayfun(@(x) regexp(x.name,'[Tt]\d{1,2}.csv'), files,'UniformOutput',0));
names = {files(ind).name};

% 128 elect; 90 tiempo; 6 condiciones; 20 sujetos
% CONDICIONES:
% 1.PreRP -LoPred, 2.PreRP -HiPred
% 3.PostRP-LoPred, 4.PostRP-HiPred
% 5.common-LoPred, 6.common-HiPred
media  = nan(128,90,6,20);
times = nan(1,90);
% sntTypePrePost = [2 2 0 0 1 1];
for i = 1:length(names)
    % load
    n = names(i); n = n{1};
    time = replace(n,'T',''); time = replace(time,'t',''); time=str2double(replace(time,'.csv',''));
    time
    T = readtable(n);
    times(time) = T.time(1);
        
    % divido las 6 condiciones en base al tipo de oracion y median split
    T.cond = nan(size(T,1),1); 
    medianas = varfun(@(x) quantile(x,[0.2 0.8]), T, 'GroupingVariables',{'sntTypePrePost'},'InputVariables',@isnumeric);
    for cond = 1:6
        type = ceil(cond/2);
        ind = (T.sntTypePrePost == type-1) & ...
              ( (T.pred <= medianas.Fun_pred(type,1) &  mod(cond,2)) |...
                (T.pred >= medianas.Fun_pred(type,2) & ~mod(cond,2)));
              
%         sum(ind)
        T(ind,'cond') = table(repmat(cond, sum(ind),1));            
    end
    
    % Calculo medias por sujeto y condicion
    mediasPorSuj = varfun(@mean, T, 'GroupingVariables',{'suj_id' 'cond'},'InputVariables',@isnumeric);
    
    % Guardo todo en media
    elect={};for e=1:128; elect{end+1}=['mean_E' num2str(e)];end
    Tnew = mediasPorSuj(:,[elect, {'suj_id','cond'}]);
    for j = 1:size(Tnew)
        c = table2array(Tnew(j,'cond'));
        s = table2array(Tnew(j,'suj_id'));
        media(:,time,c,s) = table2array(Tnew(j,elect));
    end
end
%% Plot remef
type = {'PreRP', 'PostRP', 'Common'};
pred = {'LoPred', 'HiPred'};
m = mean(media,4);
se= std(media,0,4)/sqrt(size(media,4));

figure(2);clf
for i = 1:6
    subplot(3,2,i)
    hold on
        imagesc(times, 1:128, m(:,:,i), [-1 1])
        title([type{ceil(i/2)} ' - ' pred{~mod(i,2)+1}])
        plot([0 0],[128 1],'k--','LineWidth',2)
        xlim([-100 600])
        ylim([1 128])
    hold off
end

figure(4);clf
for i = 1:2:6
    subplot(3,1,i/2+.5)
    hold on
        imagesc(times, 1:128, m(:,:,i)-m(:,:,i+1), [-.5 .5])
        title([type{ceil(i/2)} ' - ' pred{~mod(i,2)+1}])
        plot([0 0],[128 1],'k--','LineWidth',2)
        xlim([-100 600])
%         ylim([128 1])
    hold off
end

figure(3);clf
e = 60:68;
YLIM = [-4 4];
    subplot(1,3,1)
        hold on
        plot(times, mean(m(e,:,1),1),'r-', 'LineWidth',2)
        plot(times, mean(m(e,:,2),1),'r--','LineWidth',2)        
        plot([0 0],[-1.5 2.8],'k--','LineWidth',2)
        plot([-100 600],[0 0],'k--','LineWidth',2)
        xlim([-100 600])
        ylim(YLIM)
        legend({'LoPred', 'HiPred'})
        title('PreRP')

    subplot(1,3,2)
        hold on
        plot(times, mean(m(e,:,3),1),'m-', 'LineWidth',2)
        plot(times, mean(m(e,:,4),1),'m--','LineWidth',2)
        plot([0 0],[-1.5 2.8],'k--','LineWidth',2)
        plot([-100 600],[0 0],'k--','LineWidth',2)
        xlim([-100 600])
        ylim(YLIM)
        legend({'LoPred', 'HiPred'})
        title('PostRP')

    subplot(1,3,3)
        hold on
        plot(times, mean(m(e,:,5),1),'b-', 'LineWidth',2)
        plot(times, mean(m(e,:,6),1),'b--','LineWidth',2)
        plot([0 0],[-1.5 2.8],'k--','LineWidth',2)
        plot([-100 600],[0 0],'k--','LineWidth',2)
        xlim([-100 600])
        ylim(YLIM)
        legend({'LoPred', 'HiPred'})
        title('Common')
%% Load remef type-pos (para topo)
lm_Conf.datapath = '/media/brunobian/DATABRUNO/Co-registro_2018/Data/';
cd(lm_Conf.datapath)

% cd('csv_tesis_center_R_causal/conICA/') % Original
% cd('LMM/results/paper/M1_remef/suj_id/') % M1
cd('LMM/results/paper/M2_remef_Beta_(keep_pos)/suj_id/') % M2

files = dir("./");
ind = cellfun(@(y) ~isempty(y), arrayfun(@(x) regexp(x.name,'[Tt]\d{1,2}.csv'), files,'UniformOutput',0));
names = {files(ind).name};

% 128 elect; 90 tiempo; 6 condiciones; 20 sujetos
% CONDICIONES:

media  = nan(128,90,6,20);
times = nan(1,90);
media = [];
for i = 1:length(names)
    % load
    n = names(i); n = n{1};
    time = replace(n,'T',''); time = replace(time,'t',''); time=str2double(replace(time,'.csv',''));
    time
    T = readtable(n);
    T.sntType = T.sntType+0.5;
    times(time) = T.time(1);
    
    
    % ERP por posisción para topos
    T.condRel = repmat({'a'},size(T,1),1); 
    T.condAbs = repmat({'a'},size(T,1),1); 
    
    pares = [-4,-3,-2,-1;0,0,0,0;1,1,2,3;4,4,5,6;7,8,9,10];
    
	posiciones = 1:length(pares); 
    for j = posiciones
        
        posM = pares(j,:);
        posC = posM+5;

        indM_rel = ismember(T.posRelRP,posM) & ~T.sntType;
        indM_abs = ismember(T.pos, posC) & ~T.sntType;
        indC     = ismember(T.pos, posC) &  T.sntType;

        condM_rel = ['memory_rel_' strrep(num2str(posM(1)),'-','m')];
        condM_abs = ['memory_abs_' num2str(posC(1)+1)];        
        condC     = ['common_' num2str(posC(1)+1)];

        T(indM_rel,'condRel') = repmat({condM_rel}, sum(indM_rel),1);
        T(indM_abs,'condAbs') = repmat({condM_abs}, sum(indM_abs),1);
        T(indC,'condAbs')     = repmat({condC}, sum(indC),1);
    end
    
    % Calculo medias por sujeto y condicion
    mediasPorSuj_rel = varfun(@mean, T, 'GroupingVariables',{'suj_id' 'condRel'},'InputVariables',@isnumeric);
    mediasPorSuj_abs = varfun(@mean, T, 'GroupingVariables',{'suj_id' 'condAbs'},'InputVariables',@isnumeric);
    
    mediasPorSuj_rel.Properties.VariableNames(2) = {'cond'};
    mediasPorSuj_abs.Properties.VariableNames(2) = {'cond'};
    mediasPorSuj = [mediasPorSuj_abs; mediasPorSuj_rel];

    % Guardo todo en media
    elect={};for e=1:128; elect{end+1}=['mean_E' num2str(e)];end
    Tnew = mediasPorSuj(:,[elect, {'suj_id','cond'}]);
    for j = 1:size(Tnew)
        c = table2array(Tnew(j,'cond'));
        s = table2array(Tnew(j,'suj_id'));
        media.(c{:})(:,time,s) = table2array(Tnew(j,elect));
    end
end
media = rmfield(media, 'a' );

f = fields(media);
for i=1:length(f)
    bl= mean(media.(f{i})(:,times<-25,:),2);
    media.(f{i}) = media.(f{i})-bl;
end



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


%% 






fs = {'AlphaHigh','AlphaLow','Alpha','Beta','Theta','conICA'};
for iFolder = 6
lm_Conf.datapath = ['/media/brunobian/DATABRUNO/Co-registro_2018/Data/12-edited_forERPs/' fs{iFolder} '/'];
lm_Conf.outpath = ['/media/brunobian/DATABRUNO/Co-registro_2018/Data/13-remef/' fs{iFolder} '/'];

cd([lm_Conf.datapath])
[lm_Conf, SUJ]= definePath2020(lm_Conf);

tic
ERP=[];
for iSuj = 1:length(SUJ)
    % Load EEGLab matrix with information from the experiment into DATA
    EEG = pop_loadset([lm_Conf.datapath  SUJ(iSuj).fileName]); 
    EEG = pop_resample(EEG,128);
   
    DATA = EEG.epoch;
    RP          = [DATA.RP];
    PostRP      = RP > 0; 
    
    pred    = [DATA.pred];
    sntType = [DATA.sntType];
    
    sntTypePrePost = sntType;
    sntTypePrePost(sntType==0 & PostRP==0) = 0;
    sntTypePrePost(sntType==0 & PostRP==1) = 1;
    sntTypePrePost(sntType==1)             = 2;
    
    beta_intercept = values.p.Intercept(:,:,1);
    beta_predType0 = values.p.pred_sntTypePrePost0(:,:,1);
    beta_predType1 = values.p.pred_sntTypePrePost1(:,:,1);
    beta_predType2 = values.p.pred_sntTypePrePost2(:,:,1);
    
    for iepoch = 1:size(EEG.data,3)
        EEG.data(:,:,iepoch) = ...
        beta_intercept + ...
        beta_predType0 * pred(iepoch) * double(sntTypePrePost(iepoch)==0) +...
        beta_predType1 * pred(iepoch) * double(sntTypePrePost(iepoch)==1) +...
        beta_predType2 * pred(iepoch) * double(sntTypePrePost(iepoch)==2);
    end
    
    EEG = eeg_checkset( EEG );
    pop_saveset(EEG, [lm_Conf.outpath SUJ(iSuj).fileName]);


    nombres_vars = {'LoProvPre',  'HiProvPre',...
                    'LoProvPost', 'HiProvPost',...
                    'LoCommon',   'HiCommon'};
    
    limites = quantile(pred, [0 0.45 0.55 1]);                
    
    limites(1) = limites(1)-1;
    limites01 = [limites(1:2); zeros(1,2)];
    limites02 = [limites(3:4); zeros(1,2)];
    limites11 = [limites(1:2); ones(1,2) ];
    limites12 = [limites(3:4); ones(1,2) ];
    limites21 = [limites(1:2); ones(1,2)*2];
    limites22 = [limites(3:4); ones(1,2)*2];
    
    limites = cat(3, limites01, limites02, ...
                     limites11, limites12, ...
                     limites21, limites22);

    for c = 1:length(nombres_vars)
        variable = nombres_vars{c};
        cond1 = pred > limites(1,1,c) & pred <= limites(1,2,c);
        cond2 = sntTypePrePost == limites(2,1,c);
        indtrial = find(cond1 & cond2);

        ERP(iSuj).pred_type.(variable).m = nanmean(EEG.data(:,:,indtrial),3);
        ERP(iSuj).pred_type.(variable).n = size(EEG.data(:,:,indtrial),3); 
        ERP(iSuj).pred_type.(variable).trials = EEG.data(:,:,indtrial); 
    end

end
end

%% TCE clustering

tfcePath = '/home/brunobian/Documents/Repos/ept_TFCE-matlab';
addpath(genpath(tfcePath),'-end'); 

data = values.t;
% f = fieldnames(data);
% for v=1:length(f)
%     data.(f{v}) = cat(3,data.(f{v}),data.(f{v}));
% end

cfg = [];
cfg.GUI            = 0;
cfg.nPerm          = size(data.Intercept,3)-1; % default number of permutations
cfg.ElecFile       = CHANS.chanlocs; 
% cfg.rSample        = 128; 
cfg.saveResultName = ['ept_Results_', date, '.mat'];
cfg.savePermName   = ['ept_Permutations_', num2str(cfg.nPerm), '_', date, '.mat'];
% cfg.type           = 'across'; % wi = permutate whitin subjects | across = permutate across subjects
cfg.plots          = 0; % change to '1' to show significance plots after calculation
cfg.flag_ft        = 0; 
cfg.flag_tfce      = 1;
cfg.flag_save      = 1;
cfg.ChN            = [];
cfg.E_H            = [1.8 2.5];

Results = ept_TFCE_LMM(data, cfg);
Results
% Plot TFCE RESULTS
close all
colores = jet;
fields = fieldnames(Results.TFCE_Obs);
alphaTval = 1;
load('times')
p = 0.10;
for iv = 1:length(fields) 

    v = fields{iv};

    tval_obs  = Results.TFCE_Obs.(v);
    tval_perm = tval_obs .* (Results.P_Values.(v) < p);
    
    figure();
    subplot(1,2,1)
    set(gcf,'Color','w')
        col_ax = [-1 1]*5000;
        mask = abs(tval_obs);
        elect = 1:size(tval_obs,1);
        tmp = tval_obs;
        imagesc(t, elect, tmp(elect,:), col_ax);
        hold on
            plot([0 0], [1 elect(end)], 'w--')
%             set(gca, 'YTickLabel',{})
%             set(gca, 'XTickLabel',{})
            box on
        hold off
        tit = regexprep(v, '_', '-');
        title([tit 'TFCE orig'])
        f=getframe(gca);
    subplot(1,2,2)
        mask = abs(tval_perm);
        elect = 1:size(tval_perm,1);
        tmp = tval_perm.*mask;
        imagesc(t, elect, tmp(elect,:), col_ax);
        hold on
            plot([0 0], [1 elect(end)], 'w--')
%             set(gca, 'YTickLabel',{})
%             set(gca, 'XTickLabel',{})
            box on
        hold off
        tit = regexprep(v, '_', '-');
        title([tit 'TFCE 100 perms'])
        f=getframe(gca);
    
    
        saveas(gcf, [lm_Conf.lmmOutPath tit '.png'])
end

%% Fig 3:  ERP clásicos vs type-pred 
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

UseChan = {'Cz'};
indch= [1];

figure(4);clf
title('Type-Predictability')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-100 600];
y_lim=[-3 3];
for ind = 1:length(indch)
%     subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        [~,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type.LoProvPre.m(indch(ind),:)),  squeeze(erp.pred_type.LoProvPre.e(indch(ind),:)), [1 0 0],.1,'-');
        n1   = ['LoProverb PreRP'];   
        [~,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type.HiProvPre.m(indch(ind),:)),  squeeze(erp.pred_type.HiProvPre.e(indch(ind),:)), [1 0 0],.1,'--');
        n2   = ['HiProverb PreRP'];
        [~,h3] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type.LoProvPost.m(indch(ind),:)),  squeeze(erp.pred_type.LoProvPost.e(indch(ind),:)), [.5 .5 0],.1,'-');
        n3   = ['LoProverb PostRP'];   
        [~,h4] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type.HiProvPost.m(indch(ind),:)),  squeeze(erp.pred_type.HiProvPost.e(indch(ind),:)), [.5 .5 0],.1,'--');
        n4   = ['HiProverb PostRP'];
        [~,h5] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type.LoCommon.m(indch(ind),:)),squeeze(erp.pred_type.LoCommon.e(indch(ind),:)), [0 0 1],.1,'-');
        n5   = ['LoCommon'];
        [~,h6] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type.HiCommon.m(indch(ind),:)),squeeze(erp.pred_type.HiCommon.e(indch(ind),:)), [0 0 1],.1,'--');
        n6   = ['HiCommon'];
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1 h2 h3 h4 h5 h6], n1,n2,n3,n4,n5,n6)
%% Plot all AIC
addpath(genpath(lm_Conf.ftpath),'-end'); 

lm_Conf.custonFc = '/home/brunobian/Documents/Repos/corregistro/Analysis_bb/2020/';
addpath(genpath(lm_Conf.custonFc))


lm_Conf.nTimes = 89;
% permToLoad = {lm_Conf.permutationVariable};
permToLoad = {'across'};

f={'0-Base',...
'1-Base_RP',...
'2-Theta',...
'3-AlphaLow',...
'4-AlphaHigh',...
'5-Beta'};

figure();
hold on
for i =1:length(f)
lm_Conf.lmmOutPath         = [lm_Conf.datapath '../LMM/results/tesis/'...
                              f{i}...
                              '/'];
lm_Conf.matricesLoadedPath = [lm_Conf.datapath '../LMM/matrices/'];
values = lm_loadLmmData(permToLoad, lm_Conf);
load coords_LNI_128_toreplaceinEEG

% Plot AICs
h = hist(values.AIC.AIC(:),100);
plot(h)
end
figure();imagesc(values.AIC.AIC)



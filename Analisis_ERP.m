%% Settings
clear all
close all
% addpath('/home/brunobian/Documentos/toolbox/my_functions/');
% addpath('/home/brunobian/Documentos/toolbox/my_functions/coords/');
eeglabpath  = '/home/brunobian/Documents/Repos/eeglab14_1_1b';
addpath(genpath(eeglabpath),'-end'); 
ftpath      = '/home/brunobian/Documents/Repos/fieldtrip-20180513';
addpath(genpath(ftpath),'-end'); 

% cstpath = '/home/brunobian/Documents/Repos/corregistro/Analysis_bb/2020/';
cstpath = '/media/brunobian/DATABRUNO/corregistro/Analysis_bb/2020/';
addpath(genpath(cstpath),'-end'); 

clear fp; fp = FP_EEGanalysis2020;

eeglab

% datapath = '/media/brunobian/ExtraDrive1/Co-registro_2018/Data/';
datapath = '/media/brunobian/DATABRUNO/Co-registro_2018/Data/';
cd(datapath)
%% Generate ERPs per subj
clear fp; fp = FP_EEGanalysis2020;

[BL_range, ERP, data_filtros, win, dist, zapato] = fp.initVars();
propCorrectas = []; 

fs = {'conICA','Theta','AlphaLow','AlphaHigh','Alpha','Beta','remef'};
for iFolder = 1
folder = ['11-epochedFirstPassCausal/' fs{iFolder} '/'];
folderOut = ['12-edited_forERPs/' fs{iFolder} '/'];

names  = fp.loadSubjects(folder, 1, 0, '.set');
%for su=1:length(names)
indSujs = [2:3 5:12 14:15 17:24];
for su = indSujs 

    fprintf('%5s ---> Calculando ERPs...\n', names{su})

    EEG = pop_loadset([datapath folder names{su}]);
    EEG = fp.addIndProv(EEG);
    
    DATA = EEG.info.DATA;
    propCorrectas = [propCorrectas mean([DATA.correcta])];

    EEG = pop_select(EEG,'channel', 1:128);
    EEG = eeg_checkset( EEG );

    EEG = pop_select(EEG, 'time', [-.2 .6]);
%     EEG = pop_rmbase(EEG, [-100 -50]);
    EEG = pop_rmbase(EEG, [-100 0]);

    %%% Filter %%%%
    if strcmpi(fs{iFolder}, 'conICA')
        [gramFilt, voltFilt, durFilt, filtData, indbadepoch, EEG] = fp.filter(EEG, su);
        save([datapath folderOut '../filtros/' names{su}(1:end-4)], 'indbadepoch','voltFilt');        
    else 
        load([datapath folderOut '../filtros/' names{su}(1:end-4)]);
        for i = 1:length(EEG.epoch)
            EEG.epoch(i).filter = indbadepoch(i);
        end
    end

    armarERP = 1;
    campos  = {'prevFixDur','fixDur','saccDur', 'pred',  'predNext',...
               'predPrev', 'pred_type', 'sntc','posAbs','pre_post', ...
               'posRel', 'all', 'fixRank'};
    campos  = {'pred_type'};
    if armarERP
        for i = 1:length(campos)
            ERP = fp.generateErp(EEG, DATA, ERP, su, indbadepoch, win, campos{i}, voltFilt);
        end
%         zapato = fp.generateZapato(EEG, 115, indbadepoch, zapato);
    end
    guardar = 0;
    if guardar
        EEG = pop_select(EEG, 'trial',find([EEG.epoch.filter]));
        pop_saveset(EEG, [datapath folderOut names{su}]);
    end
    dimigenFig = 1;
    if dimigenFig
        elects = 37:44;

        if su==indSujs(1)
            im= []; pfixDur=[]; fixDur=[]; nfixDur=[];
        end
        
        tmp = [EEG.epoch(indbadepoch).prevFixDur];
        pfixDur = [pfixDur tmp];
        
        tmp = [EEG.epoch(indbadepoch).fixDur];
        fixDur = [fixDur tmp];
        
        tmp = [EEG.epoch(indbadepoch).nextFixDur];
        nfixDur = [nfixDur tmp];
        
        tmp = squeeze(nanmean(EEG.data(elects,:,indbadepoch),1));
        im = cat(2, im, tmp);
    end
end
end

if armarERP
    ERP([1 4 13 16])=[];
end
disp('LISTO!!')
%% Calculo erp dif pre-post
for i = 1:length(ERP)
    posrel   = 0;
    pre_post = 0;
    posrel_tipo = 1;
    posabs_tipo = 1;
    
    if pre_post
        ERP(i).diff_pre_post.memory.m = ERP(i).pre_post.memory_post.m - ERP(i).pre_post.memory_pre.m;
        ERP(i).diff_pre_post.memory.n = [ERP(i).pre_post.memory_post.n, ERP(i).pre_post.memory_pre.n];
        ERP(i).diff_pre_post.common.m = ERP(i).pre_post.common_post.m - ERP(i).pre_post.common_pre.m;
        ERP(i).diff_pre_post.common.n = [ERP(i).pre_post.common_post.n,  ERP(i).pre_post.common_pre.n];
    end
    if posrel
        ERP(i).diff_Posrel.pos_m2.m = ERP(i).PosRel.pos_m2.m - ERP(i).PosRel.pos_0.m;
        ERP(i).diff_PosRel.pos_m1.m = ERP(i).PosRel.pos_m1.m - ERP(i).PosRel.pos_0.m;
        ERP(i).diff_PosRel.pos_0.m  = ERP(i).PosRel.pos_0.m  - ERP(i).PosRel.pos_0.m;
        ERP(i).diff_PosRel.pos_1.m  = ERP(i).PosRel.pos_1.m  - ERP(i).PosRel.pos_0.m;
        ERP(i).diff_PosRel.pos_2.m  = ERP(i).PosRel.pos_2.m  - ERP(i).PosRel.pos_0.m;
        
        ERP(i).diff_PosRel.pos_m2.n = ERP(i).PosRel.pos_m2.n;
        ERP(i).diff_PosRel.pos_m1.n = ERP(i).PosRel.pos_m1.n;
        ERP(i).diff_PosRel.pos_0.n  = ERP(i).PosRel.pos_0.n;
        ERP(i).diff_PosRel.pos_1.n  = ERP(i).PosRel.pos_1.n;
        ERP(i).diff_PosRel.pos_2.n  = ERP(i).PosRel.pos_2.n;
    end
    if posrel_tipo
        ERP(i).diff_PosRel_tipo.pos_m2.m = ERP(i).PosRel_tipo.pos_m2_m.m - ERP(i).PosRel_tipo.pos_m2_c.m;
        ERP(i).diff_PosRel_tipo.pos_m1.m = ERP(i).PosRel_tipo.pos_m1_m.m - ERP(i).PosRel_tipo.pos_m1_c.m;
        ERP(i).diff_PosRel_tipo.pos_0.m  = ERP(i).PosRel_tipo.pos_0_m.m  - ERP(i).PosRel_tipo.pos_0_c.m;
        ERP(i).diff_PosRel_tipo.pos_1.m  = ERP(i).PosRel_tipo.pos_1_m.m  - ERP(i).PosRel_tipo.pos_1_c.m;
        ERP(i).diff_PosRel_tipo.pos_2.m  = ERP(i).PosRel_tipo.pos_2_m.m  - ERP(i).PosRel_tipo.pos_2_c.m;
        
        ERP(i).diff_PosRel_tipo.pos_m2.n = ERP(i).PosRel_tipo.pos_m2_m.n;
        ERP(i).diff_PosRel_tipo.pos_m1.n = ERP(i).PosRel_tipo.pos_m1_m.n;
        ERP(i).diff_PosRel_tipo.pos_0.n  = ERP(i).PosRel_tipo.pos_0_m.n;
        ERP(i).diff_PosRel_tipo.pos_1.n  = ERP(i).PosRel_tipo.pos_1_m.n;
        ERP(i).diff_PosRel_tipo.pos_2.n  = ERP(i).PosRel_tipo.pos_2_m.n;
    end
    if posabs_tipo
%         ERP(i).diff_PosAbs_tipo.pos2.m  = ERP(i).PosAbs_tipo.pos2_m.m  - ERP(i).PosAbs_tipo.pos2_c.m;
%         ERP(i).diff_PosAbs_tipo.pos3.m  = ERP(i).PosAbs_tipo.pos3_m.m  - ERP(i).PosAbs_tipo.pos3_c.m;
        ERP(i).diff_PosAbs_tipo.pos2.m  = ERP(i).PosAbs_tipo.pos23_m.m  - ERP(i).PosAbs_tipo.pos23_c.m;

        ERP(i).diff_PosAbs_tipo.pos4.m  = ERP(i).PosAbs_tipo.pos4_m.m  - ERP(i).PosAbs_tipo.pos4_c.m;
        ERP(i).diff_PosAbs_tipo.pos5.m  = ERP(i).PosAbs_tipo.pos5_m.m  - ERP(i).PosAbs_tipo.pos5_c.m;
        ERP(i).diff_PosAbs_tipo.pos6.m  = ERP(i).PosAbs_tipo.pos6_m.m  - ERP(i).PosAbs_tipo.pos6_c.m;

%         ERP(i).diff_PosAbs_tipo.pos7.m  = ERP(i).PosAbs_tipo.pos7_m.m  - ERP(i).PosAbs_tipo.pos7_c.m;
%         ERP(i).diff_PosAbs_tipo.pos8.m  = ERP(i).PosAbs_tipo.pos8_m.m  - ERP(i).PosAbs_tipo.pos8_c.m;
        ERP(i).diff_PosAbs_tipo.pos7.m  = ERP(i).PosAbs_tipo.pos78_m.m  - ERP(i).PosAbs_tipo.pos78_c.m;
        
%         ERP(i).diff_PosAbs_tipo.pos9.m  = ERP(i).PosAbs_tipo.pos9_m.m  - ERP(i).PosAbs_tipo.pos9_c.m;
%         ERP(i).diff_PosAbs_tipo.pos10.m = ERP(i).PosAbs_tipo.pos10_m.m - ERP(i).PosAbs_tipo.pos10_c.m;
        ERP(i).diff_PosAbs_tipo.pos9.m  = ERP(i).PosAbs_tipo.pos910_m.m  - ERP(i).PosAbs_tipo.pos910_c.m;
        
%         ERP(i).diff_PosAbs_tipo.pos2.n = ERP(i).PosAbs_tipo.pos2_c.n;
%         ERP(i).diff_PosAbs_tipo.pos3.n = ERP(i).PosAbs_tipo.pos3_c.n;
        ERP(i).diff_PosAbs_tipo.pos2.n = ERP(i).PosAbs_tipo.pos23_c.n;
        ERP(i).diff_PosAbs_tipo.pos4.n = ERP(i).PosAbs_tipo.pos4_c.n;
        ERP(i).diff_PosAbs_tipo.pos5.n = ERP(i).PosAbs_tipo.pos5_c.n;
        ERP(i).diff_PosAbs_tipo.pos6.n = ERP(i).PosAbs_tipo.pos6_c.n;
%         ERP(i).diff_PosAbs_tipo.pos7.n = ERP(i).PosAbs_tipo.pos7_c.n;
%         ERP(i).diff_PosAbs_tipo.pos8.n = ERP(i).PosAbs_tipo.pos8_c.n;        
        ERP(i).diff_PosAbs_tipo.pos7.n = ERP(i).PosAbs_tipo.pos78_c.n;
        
%         ERP(i).diff_PosAbs_tipo.pos9.n = ERP(i).PosAbs_tipo.pos9_c.n;
%         ERP(i).diff_PosAbs_tipo.pos10.n = ERP(i).PosAbs_tipo.pos10_c.n;
         ERP(i).diff_PosAbs_tipo.pos9.n = ERP(i).PosAbs_tipo.pos910_c.n;        
    end
    
end
%% Calculo ERPs totales
% [BL_range, ~, data_filtros, win,~] = fp.initVars();
%load ERP_freq
% load DATA

load /media/brunobian/DATABRUNO/Co-registro_2018/Data/codigo/expe_corregistro_2018/my_functions/coords_LNI_128_toreplaceinEEG
erp = [];
% campos  = {'prevFixDur','fixDur','saccDur', ...
%            'pred',  'predNext','predPrev', ...
%            'pred_type', 'sntc','posAbs','pre_post', 'posRel',...
%            'all','fixRank',...
%            'diff_PosAbs_tipo','diff_PosRel_tipo','PosAbs_tipo','PosRel_tipo'};
campos  = {'pred_type'};
erp.Nsuj  = length(ERP);

erp.times.clasicos = EEG.times;
% erp.times.rel_BL   = ERP(1).abspos.tipo0(1).t;

for c = 1:length(campos)
    camp  = campos{c};
    
    niveles = fieldnames(ERP(1).(campos{c}));

    for n = 1:length(niveles)
        Nivel    = niveles{n};

        n_times = size(ERP(1).(camp).(Nivel).m,2);
        n_elect = size(ERP(1).(camp).(Nivel).m,1);
        
        temp    = nan(n_elect, n_times, erp.Nsuj);
        N = 0;
        for su = 1:length(ERP)
            temp(:,:,su) = ERP(su).(camp).(Nivel).m;
            N = N + ERP(su).(camp).(Nivel).n;
        end
        
        erp.(camp).(Nivel).m = squeeze(nanmean(temp,3));
        erp.(camp).(Nivel).e = squeeze(nanstd(temp,0,3))/sqrt(size(temp,3));
        erp.(camp).(Nivel).Nsuj = sum(~isnan(temp(1,1,:)));
        erp.(camp).(Nivel).Ntr  = N;
    end 
end
save('erp_totales','erp')
% save('erp_totales','erp')
disp('listo!!')

%% Load erps
% load('erp_Theta')
% load('erp_AlphaLow')
% load('erp_AlphaHigh')
% load('erp_Beta')
load('erp_totales')

savePath = '/media/brunobian/DATABRUNO/Co-registro_2018/Data/figs/';

%% Fig 0:  ERP clásicos (all)
errorbars = 1;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(2);clf
title('Predictability')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-100 600];
y_lim=[-3 3];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        if errorbars == 1
        % Con las barras de error
            [~,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.all.all.m(indch(ind),:)), squeeze(erp.all.all.e(indch(ind),:)), [1 0  0], .1);
        else
        % Sin las barras de error
            [~,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.all.all.m(indch(ind),:)), ceros(1,:), [1 0  0], .1);            
        end
    n1 = ['All Epochs (' num2str(erp.all.all.Ntr) ')'];
    hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1], n1)
f=getframe(gca);
%% Fig 1:  ERP clásicos vs pred 
errorbars = 1;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(2);clf
title('Predictability')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-100 600];
y_lim=[-3 3];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        if errorbars == 1
        % Con las barras de error
            [~,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred.q1.m(indch(ind),:)), squeeze(erp.pred.q1.e(indch(ind),:)), [1 0  0], .1);
            [~,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred.q2.m(indch(ind),:)), squeeze(erp.pred.q2.e(indch(ind),:)), [0 .5 0], .1 );
            [~,h3] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred.q3.m(indch(ind),:)), squeeze(erp.pred.q3.e(indch(ind),:)), [0 0  1], .1 );
        else
        % Sin las barras de error
            [~,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred.q1.m(indch(ind),:)), ceros(1,:), [1 0  0], .1);
            [~,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred.q2.m(indch(ind),:)), ceros(1,:), [0 0 1], .1 );
            [~,h3] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred.q3.m(indch(ind),:)), ceros(1,:), [0 .5 0], .1 );
        end
    n1 = ['Low Pred (' num2str(erp.pred.q1.Ntr) ')'];
    n2 = ['Mid Pred (' num2str(erp.pred.q2.Ntr) ')'];
    n3 = ['Hi Pred (' num2str(erp.pred.q3.Ntr) ')'];
    hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1 h2 h3], n1, n2, n3)
f=getframe(gca);
% saveas(gcf, [savePath 'ERP_pred.svg'])
%% Fig 2:  ERP clásicos vs type
errorbars = 1;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(3);clf
title('Type')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-100 600];
y_lim=[-.2 .2];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

            [~,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.sntc.Memory.m(indch(ind),:)), squeeze(erp.sntc.Memory.e(indch(ind),:)), [1 0  0], .1);
            n1 = ['Memory-Encoded (' num2str(erp.sntc.Memory.Ntr) ')'];
            [~,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.sntc.Common.m(indch(ind),:)), squeeze(erp.sntc.Common.e(indch(ind),:)), [0 .5 0], .1 );
            n2 = ['Common (' num2str(erp.sntc.Common.Ntr) ')'];
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1 h2], n1,n2)
f=getframe(gca);
% saveas(gcf, [savePath 'ERP_type.svg'])
%% Fig 3:  ERP clásicos vs type-pred 
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(4);clf
title('Type-Predictability')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-100 600];
y_lim=[-3 3];
for ind = 1:length(indch)
%     subplot(3,3,ind)
    figure(4*10+ind);clf
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        [~,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type.LoProv.m(indch(ind),:)),  squeeze(erp.pred_type.LoProv.e(indch(ind),:)), [1 0 0],.1,'-');
        n1   = ['LoProverb (' num2str(erp.pred_type.LoProv.Ntr) ')'];   
        [~,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type.HiProv.m(indch(ind),:)),  squeeze(erp.pred_type.HiProv.e(indch(ind),:)), [1 0 0],.1,'--');
        n2   = ['HiProverb (' num2str(erp.pred_type.HiProv.Ntr) ')'];
        [~,h3] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type.LoCommon.m(indch(ind),:)),squeeze(erp.pred_type.LoCommon.e(indch(ind),:)), [0 0 1],.1,'-');
        n3   = ['LoCommon (' num2str(erp.pred_type.LoCommon.Ntr) ')'];
        [~,h4] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type.HiCommon.m(indch(ind),:)),squeeze(erp.pred_type.HiCommon.e(indch(ind),:)), [0 0 1],.1,'--');
        n4   = ['HiCommon (' num2str(erp.pred_type.HiCommon.Ntr) ')'];
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)

    if ind == length(indch)
        legend([h1 h2 h3 h4], n1,n2,n3,n4)
        f=getframe(gca);
    end
    saveas(gcf, [savePath 'FRP_pred/ERP_pred-type' UseChan{ind} '.svg'])
end
% saveas(gcf, [savePath 'ERP_pred-type.svg'])
%% Fig 4:  ERP clásicos vs predNext
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(5);clf
title('Predictability')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-200 600];
y_lim=[-3 3];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        [errorPatch,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.predNext.q1.m(indch(ind),:)), ceros(1,:), [1 0  0], .1);
        n1   = ['LoPred N+1 (' num2str(erp.predNext.q1.Ntr) ')'];   
        [errorPatch,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.predNext.q2.m(indch(ind),:)), ceros(1,:), [0 0 1], .1 );
        n2   = ['MidPred N+1 (' num2str(erp.predNext.q2.Ntr) ')'];   
        [errorPatch,h3] = niceBars_bb( erp.times.clasicos, squeeze(erp.predNext.q3.m(indch(ind),:)), ceros(1,:), [0 .5 0], .1 );
        n3   = ['HiPred N+1 (' num2str(erp.predNext.q3.Ntr) ')'];   

        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1 h2 h3], n1,n2,n3)
%% Fig 5:  ERP clásicos vs predPrev
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(6);clf
title('Predictability')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-200 600];
y_lim=[-3 3];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        [errorPatch,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.predPrev.q1.m(indch(ind),:)), ceros(1,:), [1 0  0], .1);
        n1   = ['LoPred N-1 (' num2str(erp.predPrev.q1.Ntr) ')'];   
        [errorPatch,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.predPrev.q2.m(indch(ind),:)), ceros(1,:), [0 0 1], .1 );
        n2   = ['MidPred N-1 (' num2str(erp.predPrev.q2.Ntr) ')'];   
        [errorPatch,h3] = niceBars_bb( erp.times.clasicos, squeeze(erp.predPrev.q3.m(indch(ind),:)), ceros(1,:), [0 .5 0], .1 );
        n3   = ['HiPred N-1 (' num2str(erp.predPrev.q3.Ntr) ')'];   
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1 h2 h3], n1, n2, n3)
%% Fig 6:  ERP clásicos vs FixDur
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(7);clf
title('Predictability')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-200 600];
y_lim=[-3 4];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        [errorPatch,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.fixDur.Short.m(indch(ind),:)), ceros(1,:), [1 0 0], .1);
        n1   = ['[50 175] (' num2str(erp.fixDur.Short.Ntr) ')'];   
        [errorPatch,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.fixDur.Mid.m(indch(ind),:)), ceros(1,:), [0 .5 0], .1 );
        n2   = ['[175 250] (' num2str(erp.fixDur.Mid.Ntr) ')'];   
        [errorPatch,h3] = niceBars_bb( erp.times.clasicos, squeeze(erp.fixDur.Long.m(indch(ind),:)), ceros(1,:), [0 0 1], .1 );
        n3   = ['[250 600] (' num2str(erp.fixDur.Long.Ntr) ')'];   
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
%legend([h1 h2 ], n1, n2)
legend([h1 h2 h3], n1, n2, n3)
%% Fig 7:  ERP clásicos vs SaccDur
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(7);clf
title('Predictability')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-200 600];
y_lim=[-3 4];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        [errorPatch,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.saccDur.Short.m(indch(ind),:)), ceros(1,:), [1 0 0], .1);
        n1   = ['[50 200] (' num2str(erp.saccDur.Short.Ntr) ')'];   
        [errorPatch,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.saccDur.Mid.m(indch(ind),:)), ceros(1,:), [0 .5 0], .1 );
        n2   = ['[200 350] (' num2str(erp.saccDur.Mid.Ntr) ')'];   
        [errorPatch,h3] = niceBars_bb( erp.times.clasicos, squeeze(erp.saccDur.Long.m(indch(ind),:)), ceros(1,:), [0 0 1], .1 );
        n3   = ['[350 600] (' num2str(erp.saccDur.Long.Ntr) ')'];   
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
% legend([h1 h2 ], n1, n2)
legend([h1 h2 h3], n1, n2, n3)
%% Fig 7:  pred distribution

pred = [EEG.epoch.pred];
predProv   = pred([EEG.epoch.sntType]==0);
predCommon = pred([EEG.epoch.sntType]==1);

[h1, x1] = hist(predProv);
[h2, x2] = hist(predCommon);

figure(8)
hold on
l1 = plot(x1,h1, 'LineWidth',2);
l2 = plot(x2,h2, 'LineWidth',2);
legend([l1 l2], 'Proverb', 'Common')
xlabel('logit(pred)')
% ylim([0 320])

% limites = [-1.6128   -1.2788   -0.6090    1.5798];
% for i = 1:length(limites)
%     plot([limites(i) limites(i)], [0 600], '--k')
% end

hold off
%% Fig 8:  ERP clásicos vs type-PrePost 
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(8);clf
title('Type-PrePost')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-100 600];
y_lim=[-.2 .2];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        [~,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.pre_post.memory_pre.m(indch(ind),:)), squeeze(erp.pre_post.memory_pre.e(indch(ind),:)), [1 0  0], .1,'-');
        n1   = ['Prov Pre'];   
        [~,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.pre_post.memory_post.m(indch(ind),:)), squeeze(erp.pre_post.memory_post.e(indch(ind),:)), [1 0  0], .1,'--');
        n2   = ['Prov Post'];
        [~,h3] = niceBars_bb( erp.times.clasicos, squeeze(erp.pre_post.common_pre.m(indch(ind),:)), squeeze(erp.pre_post.common_pre.e(indch(ind),:)), [0 0  1], .1,'-');
        n3   = ['Common Pre'];
        [~,h4] = niceBars_bb( erp.times.clasicos, squeeze(erp.pre_post.common_post.m(indch(ind),:)), squeeze(erp.pre_post.common_post.e(indch(ind),:)), [0 0  1], .1,'--');
        n4   = ['Common Post'];
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1 h2 h3 h4], n1,n2,n3,n4)
% legend([h1 h2], n1,n2)
%% Fig 8:  ERP clásicos vs type-DifPrePost 
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(8);clf
title('Type-DifPrePost')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-100 600];
y_lim=[-.1 .1];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        [~,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.diff_pre_post.memory.m(indch(ind),:)), squeeze(erp.diff_pre_post.memory.e(indch(ind),:)), [1 0  0], .1,'-');
        n1   = ['Prov Post-Pre'];   
        [~,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.diff_pre_post.common.m(indch(ind),:)), squeeze(erp.diff_pre_post.common.e(indch(ind),:)), [0 0  1], .1,'-');
        n2   = ['Common Post-Pre'];
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1 h2], n1,n2)
% legend([h1 h2], n1,n2)
%% Fig 9:  ERP clásicos vs posAbs 
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(8);clf
title('Type-PrePost')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-200 600];
y_lim=[-3 3];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        [h1] = plot( erp.times.clasicos, squeeze(erp.posAbs.pos2.m(indch(ind),:)), '-', 'Color', [1 0  0], 'LineWidth',2);
        n1   = ['2'];   
        [h2] = plot( erp.times.clasicos, squeeze(erp.posAbs.pos3.m(indch(ind),:)), '-', 'Color', [.8 .2  0], 'LineWidth',2);
        n2   = ['3'];   
        [h3] = plot( erp.times.clasicos, squeeze(erp.posAbs.pos4.m(indch(ind),:)), '-', 'Color', [.6 .4  0], 'LineWidth',2);
        n3   = ['4'];   
        [h4] = plot( erp.times.clasicos, squeeze(erp.posAbs.pos5.m(indch(ind),:)), '-', 'Color', [.4 .6  0], 'LineWidth',2);
        n4   = ['5'];   
        [h5] = plot( erp.times.clasicos, squeeze(erp.posAbs.pos6.m(indch(ind),:)), '-', 'Color', [.2 .8  0], 'LineWidth',2);
        n5   = ['6'];   
        [h6] = plot( erp.times.clasicos, squeeze(erp.posAbs.pos7.m(indch(ind),:)), '-', 'Color', [0 1  0], 'LineWidth',2);
        n6   = ['7'];   
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1 h2 h3 h4 h5 h6], n1,n2,n3,n4,n5,n6)
%% Fig 10:  ERP clásicos vs posRel
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(8);clf
title('Type-PrePost')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-200 600];
y_lim=[-.2 .2];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        tmp = squeeze(erp.diff_posrel.pos_m2.m(indch(ind),:));
        [h1] = plot( erp.times.clasicos, tmp, '-', 'Color', [.8 0  0], 'LineWidth',2);
        n1   = ['-2'];   
        
        tmp = squeeze(erp.diff_posrel.pos_m1.m(indch(ind),:));
        [h2] = plot( erp.times.clasicos, tmp, '-', 'Color', [.6 0  0], 'LineWidth',2);
        n2   = ['-1'];   
        
        tmp = squeeze(erp.diff_posrel.pos_0.m(indch(ind),:));
        [h3] = plot( erp.times.clasicos, tmp, '-', 'Color', [0 0  1], 'LineWidth',2);
        n3   = ['0'];   
        
        tmp = squeeze(erp.diff_posrel.pos_1.m(indch(ind),:));
        [h4] = plot( erp.times.clasicos, tmp, '-', 'Color', [0 .6  0], 'LineWidth',2);
        n4   = ['1'];   
        
        tmp = squeeze(erp.diff_posrel.pos_2.m(indch(ind),:));
        [h5] = plot( erp.times.clasicos, tmp, '-', 'Color', [0 .8  0], 'LineWidth',2);
        n5   = ['2'];   
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1 h2 h3 h4 h5], n1,n2,n3,n4,n5)
% legend([h1 h3 h5], n1,n3,n5)
%% Fig 11:  mean power vs posRel
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];

figure(11);clf
title('Type-PrePost')
x_lim=[-3.1 3.1];
y_lim=[-.08 .12];
N = 20;
t = erp.times.clasicos > 0 & erp.times.clasicos < 600;
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))

        mean_m3_m = mean(erp.posRel_tipo.pos_m3_m.m(indch(ind),t));
        ste_m3_m  = std(erp.posRel_tipo.pos_m3_m.m(indch(ind),t))/sqrt(N);
        mean_m2_m = mean(erp.posRel_tipo.pos_m2_m.m(indch(ind),t));
        ste_m2_m  = std(erp.posRel_tipo.pos_m2_m.m(indch(ind),t))/sqrt(N);
        mean_m1_m = mean(erp.posRel_tipo.pos_m1_m.m(indch(ind),t));
        ste_m1_m  = std(erp.posRel_tipo.pos_m1_m.m(indch(ind),t))/sqrt(N);
        mean_0_m  = mean(erp.posRel_tipo.pos_0_m.m(indch(ind),t));
        ste_0_m   = std(erp.posRel_tipo.pos_0_m.m(indch(ind),t))/sqrt(N);
        mean_1_m  = mean(erp.posRel_tipo.pos_1_m.m(indch(ind),t));
        ste_1_m   = std(erp.posRel_tipo.pos_1_m.m(indch(ind),t))/sqrt(N);
        mean_2_m  = mean(erp.posRel_tipo.pos_2_m.m(indch(ind),t));
        ste_2_m   = std(erp.posRel_tipo.pos_2_m.m(indch(ind),t))/sqrt(N);
        mean_3_m  = mean(erp.posRel_tipo.pos_3_m.m(indch(ind),t));
        ste_3_m   = std(erp.posRel_tipo.pos_3_m.m(indch(ind),t))/sqrt(N);

        mean_m3_c = mean(erp.posRel_tipo.pos_m3_c.m(indch(ind),t));
        ste_m3_c  = std(erp.posRel_tipo.pos_m3_c.m(indch(ind),t))/sqrt(N);
        mean_m2_c = mean(erp.posRel_tipo.pos_m2_c.m(indch(ind),t));
        ste_m2_c  = std(erp.posRel_tipo.pos_m2_c.m(indch(ind),t))/sqrt(N);
        mean_m1_c = mean(erp.posRel_tipo.pos_m1_c.m(indch(ind),t));
        ste_m1_c  = std(erp.posRel_tipo.pos_m1_c.m(indch(ind),t))/sqrt(N);
        mean_0_c  = mean(erp.posRel_tipo.pos_0_c.m(indch(ind),t));
        ste_0_c   = std(erp.posRel_tipo.pos_0_c.m(indch(ind),t))/sqrt(N);
        mean_1_c  = mean(erp.posRel_tipo.pos_1_c.m(indch(ind),t));
        ste_1_c   = std(erp.posRel_tipo.pos_1_c.m(indch(ind),t))/sqrt(N);
        mean_2_c  = mean(erp.posRel_tipo.pos_2_c.m(indch(ind),t));
        ste_2_c   = std(erp.posRel_tipo.pos_2_c.m(indch(ind),t))/sqrt(N);
        mean_3_c  = mean(erp.posRel_tipo.pos_3_c.m(indch(ind),t));
        ste_3_c   = std(erp.posRel_tipo.pos_3_c.m(indch(ind),t))/sqrt(N);

        
        x = [-3:3];
        hold on
            means = [mean_m3_m mean_m2_m mean_m1_m mean_0_m mean_1_m mean_2_m mean_3_m];
            ste   = [ste_m3_m ste_m2_m ste_m1_m ste_0_m ste_1_m ste_2_m ste_3_m];
            errorbar(x,means,ste)
            
            means = [mean_m3_c mean_m2_c mean_m1_c mean_0_c mean_1_c mean_2_c mean_3_c];
            ste   = [ste_m3_c ste_m2_c ste_m1_c ste_0_c ste_1_c ste_2_c ste_3_c];
            errorbar(x,means,ste)

        hold off
        if ind == 9
        legend('memory', 'common')   
        end
        xlabel('Relative position to RP')
        ylabel('Mean beta Power')
        hold off
        ylim(y_lim)
        xlim(x_lim)
    box on
end
%% Fig 12 topoplot posRel_tipo
campos_rel = {'pos_m4' 'pos_m3' 'pos_m2' 'pos_m1' 'pos_0' 'pos_1' 'pos_2' 'pos_3' 'pos_4'};
campos_abs = {'pos2'   'pos3'   'pos4'   'pos5'   'pos6'  'pos7'  'pos8'  'pos9' 'pos10'};
titulos_m = {'RP-4' 'RP-3' 'RP-2' 'RP-1' 'RP' 'RP+1' 'RP+2' 'RP+3' 'RP+4'};
titulos_c = {'2' '3' '4' '5' '6' '7' '8' '9' '10'};

l=length(campos_rel);
step = 50;
winSize = 50;
ts = 150:step:200;
% ts = 0:step:50;
for j = 1:(length(ts)-1)
    it = ts(j);
    figure(12);clf
    t = erp.times.clasicos > it & erp.times.clasicos <it+winSize;
    txt = [num2str(it) ' - ' num2str(it+winSize) ];
    sgtitle(txt) 

    for i = 1:l
        electrodes = [];
        elec = {electrodes}; 
        cr = campos_rel{i};
        ca = campos_abs{i};
        lim = [-1 1]*1;

        % common abs
        subplot(4,l,i)
        meanWindow = (mean(erp.PosAbs_tipo.([ca '_c']).m(:,t),2));
        topoplot(meanWindow, ...
                 CHANS.chanlocs, ...
                 'emarker2', elec, ...
                 'maplimits',lim , ...
                 'colormap', colormap('redbluecmap'));
        title(['common | pal ' titulos_c{i}])     

        % memory abs
        subplot(4,l,i+l)
        meanWindow = (mean(erp.PosAbs_tipo.([ca '_m']).m(:,t),2));
        topoplot(meanWindow, ...
                 CHANS.chanlocs, ...
                 'emarker2', elec, ...
                 'maplimits',lim , ...
                 'colormap', colormap('redbluecmap'));
        title(['memory | pal ' titulos_c{i}])     

%         % diff abs
%         subplot(4,l,i+2*l)
%         meanWindow = (mean(erp.diff_posabs_tipo.(ca).m(:,t),2));
%         topoplot(meanWindow, ...
%                  CHANS.chanlocs, ...
%                  'emarker2', elec, ...
%                  'maplimits',lim , ...
%                  'colormap', colormap('redbluecmap'));
%         title(['Diff M-C | pal ' titulos_c{i}])    
        
        % Common rel
        subplot(4,l,i+2*l)
        meanWindow = (mean(erp.PosRel_tipo.([cr '_c']).m(:,t),2));
        topoplot(meanWindow, ...
                 CHANS.chanlocs, ...
                 'emarker2', elec, ...
                 'maplimits',lim , ...
                 'colormap', colormap('redbluecmap'));
        title(['common | ' titulos_m{i}])    
        
        % memory rel
        subplot(4,l,i+3*l)
        meanWindow = (mean(erp.PosRel_tipo.([cr '_m']).m(:,t),2));
        topoplot(meanWindow, ...
                 CHANS.chanlocs, ...
                 'emarker2', elec, ...
                 'maplimits',lim , ...
                 'colormap', colormap('redbluecmap'));
        title(['memory | ' titulos_m{i}])     
    end
    set(gcf,'Color','w','Position', get(0, 'Screensize'))
    txt = num2str(j);
    if length(txt) == 1
        txt = ['0' txt];
    end
%     saveas(gcf, ['figs/gifTopoTheta/' txt '.png'])
end
%% Fig 13 topoplot posRel dif_tipo

campos = {'pos_m2' 'pos_m1' 'pos_0' 'pos_1' 'pos_2'};
titulos = {'-2' '-1' '0' '1' '2'};
figure(13);clf
for i = 1:length(campos)
    electrodes = [];
    elec = {electrodes}; 
    c = campos{i};
    lim = [-.1 .1];
    meanWindow = (mean(erp.diff_PosRel_tipo.(c).m,2));
    subplot(2,3,i)
    topoplot(meanWindow, ...
             CHANS.chanlocs, ...
             'emarker2', elec, ...
             'maplimits',lim , ...
             'colormap', colormap('redbluecmap'));
    title(['memory - common pal ' titulos{i}])     
    
    set(gcf,'Color','w')
end
% colorbar
%% Fig 14 topoplot posAbs dif_tipo
itV = 0:100:500; %100
winSize = 100; %50
for it = itV
campos = {'pos2' 'pos4' 'pos5' 'pos6' 'pos7' 'pos9'};
titulos = {'2-3' '4' '5' '6' '7-8' '9-10'};
figure(14);clf
stdDiff = [arrayfun(@(x) std(mean(x.diff_PosAbs_tipo.pos2.m(:,t),2)), ERP)',...
arrayfun(@(x) std(mean(x.diff_PosAbs_tipo.pos4.m(:,t),2)), ERP)',...
arrayfun(@(x) std(mean(x.diff_PosAbs_tipo.pos5.m(:,t),2)), ERP)',...
arrayfun(@(x) std(mean(x.diff_PosAbs_tipo.pos6.m(:,t),2)), ERP)',...
arrayfun(@(x) std(mean(x.diff_PosAbs_tipo.pos7.m(:,t),2)), ERP)',...
arrayfun(@(x) std(mean(x.diff_PosAbs_tipo.pos9.m(:,t),2)), ERP)'];

for i = 1:length(campos)
    electrodes = [];
    elec = {electrodes}; 
    c = campos{i};
    lim = [-.1 .1]*5;
    t = erp.times.clasicos > it & erp.times.clasicos <it+winSize;
    tit = [num2str(it) ' - ' num2str(it+winSize) ];
    sgtitle(tit) 
    meanWindow = (mean(erp.diff_PosAbs_tipo.(c).m(:,t),2));
    subplot(2,6,i)
    topoplot(meanWindow, ...
             CHANS.chanlocs, ...
             'emarker2', elec, ...
             'maplimits',lim , ...
             'colormap', colormap('redbluecmap'));
    title(['diff(m,c) pal ' titulos{i}])     
       
    set(gcf,'Color','w')
end
subplot(2,1,2)
sn=sqrt(20);
errorbar([2.5,4,5,6,7.5,9.5],mean(stdDiff),std(stdDiff)/sn)
title('diferencia media')

set(gcf,'Color','w','Position', get(0, 'Screensize'))
txt = 'diffPosAbs';
saveas(gcf, ['figs/TopoFRP/' txt '.png'])

figure(140);clf
set(gcf,'Color','w','Position', get(0, 'Screensize'))
errorbar([2.5,4,5,6,7.5,9.5],mean(stdDiff),std(stdDiff)/sn,'LineWidth',2)
txt = ['diffPosAbs - ' tit];
sgtitle(txt) 
print(['figs/TopoFRP/' txt '.eps'],'-depsc2');
end

%% Fig 15 topoplot fixRank
campos = {'fix1' 'fix2' 'fix3' 'fix4' 'fix5' 'fix6' 'fix7' 'fix8' 'fix9' 'fix10'};
titulos = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'};
campos = {'fix1' 'fix3' 'fix5' 'fix7' 'fix9'};
titulos = {'1-2' '3-4' '5-6' '7-8' '9-10'};

l=length(campos);
step = 150;
winSize = 150;
ts = 0:step:150;
% ts = 0:step:50;
for j = 1:(length(ts)-1)
    it = ts(j);
    figure(12);clf
    t = erp.times.clasicos > it & erp.times.clasicos <it+winSize;
    txt = [num2str(it) ' - ' num2str(it+winSize) ];
    sgtitle(txt) 

    for i = 1:l
        electrodes = [];
        elec = {electrodes}; 
        ca = campos{i};
        lim = [-1 1]*.1;

        % common abs
        subplot(2,l,i)
        meanWindow = (mean(erp.fixRank.([ca '_c']).m(:,t),2));
        topoplot(meanWindow, ...
                 CHANS.chanlocs, ...
                 'emarker2', elec, ...
                 'maplimits',lim , ...
                 'colormap', colormap('redbluecmap'));
        N = erp.fixRank.([ca '_c']).Ntr;
        title(['com | fix ' titulos{i} ' | N=' num2str(N)])     

        % memory abs
        subplot(2,l,i+l)
        meanWindow = (mean(erp.fixRank.([ca '_m']).m(:,t),2));
        topoplot(meanWindow, ...
                 CHANS.chanlocs, ...
                 'emarker2', elec, ...
                 'maplimits',lim , ...
                 'colormap', colormap('redbluecmap'));
        N = erp.fixRank.([ca '_m']).Ntr;
        title(['mem | fix ' titulos{i} ' | N=' num2str(N)])

%         % diff abs
%         subplot(4,l,i+2*l)
%         meanWindow = (mean(erp.diff_posabs_tipo.(ca).m(:,t),2));
%         topoplot(meanWindow, ...
%                  CHANS.chanlocs, ...
%                  'emarker2', elec, ...
%                  'maplimits',lim , ...
%                  'colormap', colormap('redbluecmap'));
%         title(['Diff M-C | pal ' titulos_c{i}])    

   
    end
    set(gcf,'Color','w','Position', get(0, 'Screensize'))
    txt = num2str(j);
    if length(txt) == 1
        txt = ['0' txt];
    end
%     saveas(gcf, ['figs/gifTopoTheta/' txt '.png'])
end

%% Fig 16:  ERP clásicos vs type-predDeciles 
errorbars = 0;
clc
UseChan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
indch= [100 85 68 115 1 54 7 19 36];


figure(4);clf
title('Type-Predictability')
ceros = zeros(length(erp.times.clasicos));
x_lim=[-100 600];
y_lim=[-3 3];
for ind = 1:length(indch)
    subplot(3,3,ind)
    set(gcf,'Color','w')
    hold on
        title(UseChan(ind))
        plot([0 0],y_lim,'k-'); 
        plot(x_lim,[0 0],'k-');
        plot([x_lim(1)+1 x_lim(1)+1],y_lim,'k-'); %porque no funciona bien el box on
        plot(x_lim, [y_lim(1)+0.0001 y_lim(1)+0.0001],'k-'); %porque no funciona bien el box on

        [~,h1] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type10.Prov3.m(indch(ind),:)),squeeze(erp.pred_type10.Prov3.e(indch(ind),:)), [1 0 0],.1,'-');
        n1   = ['Mem-Enc 4/5 (' num2str(erp.pred_type10.Prov3.Ntr) ')'];   
        [~,h2] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type10.Prov4.m(indch(ind),:)),squeeze(erp.pred_type10.Prov4.e(indch(ind),:)), [1 0 0],.1,'--');
        n2   = ['Mem-Enc 5/5 (' num2str(erp.pred_type10.Prov4.Ntr) ')'];
        [~,h3] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type10.Comm3.m(indch(ind),:)),squeeze(erp.pred_type10.Comm3.e(indch(ind),:)), [0 0 1],.1,'-');
        n3   = ['Common 4/5 (' num2str(erp.pred_type10.Comm3.Ntr) ')'];
        [~,h4] = niceBars_bb( erp.times.clasicos, squeeze(erp.pred_type10.Comm4.m(indch(ind),:)),squeeze(erp.pred_type10.Comm4.e(indch(ind),:)), [0 0 1],.1,'--');
        n4   = ['Common 5/5 (' num2str(erp.pred_type10.Comm4.Ntr) ')'];
        hold off
    box on
    xlim(x_lim)
    ylim(y_lim)
end
legend([h1 h2 h3 h4], n1,n2,n3,n4)
f=getframe(gca);
saveas(gcf, [savePath 'ERP_pred-type.svg'])

%% Fig 17: trials by fixDur (1 subj)

onsetPrev  = -pfixDur;
onsetNext  = fixDur;
insetNNext = nfixDur+onsetNext;


[~,ind] = sort(fixDur);
imSorted = im(:,ind)';

imProm=[];
ventana = 50;
for i=1:(length(imSorted)-ventana)
    imProm = cat(1,imProm,mean(imSorted(i:i+ventana-1,:),1));
end
disp('Listo')

t = EEG.times;
y = 1:length(im);
lim = [-1 1]*3.5;
edges = -200:10:600;

figure(17);clf
set(gcf,'Color','w')

    subplot(3,1,1:2) 
    hold on
        imagesc(t,y,imProm,lim)
        set(gca,'YDir','normal')
        plot([0,0],[min(y) max(y)], 'k-')
        plot(fixDur(ind),y, 'k--','LineWidth',2)        
        xlim([-200 600])
        ylim([min(y) max(y)])
        ylabel('EEG segments')
        set(gca,'xticklabel',{[]})
        colorbar
        box on
    hold off
    
    subplot(3,1,3) 
    hold on
        histogram(onsetNext, edges,'DisplayStyle','stairs','Normalization','count','EdgeColor','k')
        histogram(onsetPrev, edges,'DisplayStyle','stairs','Normalization','count','EdgeColor','k')
        histogram(insetNNext,edges,'DisplayStyle','stairs','Normalization','count','EdgeColor','k')
        plot([0,0],[0 1500], 'k-')
        xlim([-200 600])
        ylim([0 1500])
        ylabel('Fixations')
        xlabel('Time after fixation onset [ms]')
    hold off

colormap(jet)
saveas(gcf, [savePath 'erpsByFixDurDimigen.svg'])

%% Fig 18:  Zapato
data  = zapato.data;
variable = [zapato.info.fixDur];
t = zapato.times;
figure(9);clf

ind0 = [zapato.info.sntType] == 0;
ind1 = [zapato.info.sntType] == 1;

[counts0, centers0] = hist(variable(ind0),50);
[counts1, centers1] = hist(variable(ind1),50);


% type 0 (proverb)
subplot(4,2,1:2:5)
    data0 = data(ind0,:);
    [~,I0] = sort(variable(ind0));
    data0 = data0(I0,:);
    hold on
        imagesc(t, 1:size(data0,1), data0, [-20 20])
        plot([0 0], [1 size(data0,1)], 'k--')
        xlim([min(t), max(t)])
        ylim([0 size(data0,1)])
    hold off
    
subplot(4,2,7)
    c0 = counts0/sum(counts0);
    hold on
        plot(centers0, c0,'b')
        legend('Proverb')
        plot([0 0], [0 max(c0)], 'k--')
        ylim([0 max(c0)])
        xlim([min(t), max(t)])
    hold off

% type 1 (Common)
subplot(4,2,2:2:6)
    data1 = data(ind1,:);
    [~,I1] = sort(variable(ind1));
    data1 = data1(I1,:);
    hold on
        imagesc(zapato.times, 1:size(data1,1), data1, [-20 20])
        plot([0 0], [1 size(data1,1)], 'k--')
        xlim([min(t), max(t)])
        ylim([0 size(data1,1)])
    hold off

subplot(4,2,8)
    c1 = counts1/sum(counts1);
    hold on
        plot(centers1, c1,'r')
        legend('Common')
        plot([0 0], [0 max(c1)], 'k--')
        ylim([0 max(c1)])
        xlim([min(t), max(t)])
    hold off
    
%%
close all
fix_p = [zapato.info.prevFixDur];
fix_n = [zapato.info.fixDur];

d = 100000;
i_cercanas = cercanas(fix_n, fix_p,d);

figure();
hold on
    scatter(fix_n(i_cercanas), fix_p(i_cercanas))
    xlabel('FixDur Pal N')
    ylabel('FixDur Pal N-1')
    xlim([0 600])
    ylim([0 600])

%%
s = zapato.info(i_cercanas);
s.words



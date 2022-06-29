classdef FP_EEGanalysis2020
    properties
    end
    methods(Static)

        % Pre Analysis
        function names = loadSubjects(folder, expe, pre, type)
            % Load filenames (*.type) from folder.
            % if expe: load suj.bdf/set
            % if pre: load suj_pre1/2.bdf/set
            names = dir(folder); 
            names = {names(:).name};
            
            names = names(1,3:end);
            names = names(cellfun(@(x) ~isempty(x),regexp(names, type)));
            
            if pre == 1 && expe == 1
                return
            elseif pre == 1 
                names = names(cellfun(@(x) ~isempty(x),regexp(names, '_pre')));
            elseif expe == 1 
                names = names(cellfun(@(x) isempty(x),regexp(names, '_pre')));
            end  
        end        
        function [BL_range, ERP, data_filtros, win, dist] = initializeVars()
            
            BL_range = [-75 25];
            ERP = []; 
            data_filtros = [];

            win = [];

            % Calculado a partir de los topoplots
            win.P200r.time  = [200 300];
            win.P200r.elect = [36:45];

            win.P200l.time  = [200 300];
            win.P200l.elect = [7:12 125:128];

            win.N400.time   = [300 450];
            win.N400.elect  = [1:6 17:21 30:35 50:53  65 66 87 97 98 110:114 124];
        %     win.N400.elect  = [1:6 18:20 31:35 50:53  66 87 97 98 110:114 124];


            win.P600.time   = [500 650];
        %     win.P600.time   = [600 700];
            win.P600.elect  = [3:9 16:22 29:32 34:38 44 45 47:51 112 113 121:126];
            win.P600.elect  = [3:9 16:22 29:32 35:38 44 45 125:126];


            dist.pred.preMJ_p = [];
            dist.pred.posMJ_p = [];
            dist.freq.preMJ_p = [];
            dist.freq.posMJ_p = [];
            dist.pos.preMJ_p  = [];
            dist.pos.posMJ_p  = [];
            dist.pred.preMJ_r = [];
            dist.pred.posMJ_r = [];
            dist.freq.preMJ_r = [];
            dist.freq.posMJ_r = [];
            dist.pos.preMJ_r  = [];
            dist.pos.posMJ_r  = [];
        end

        
        % Analysis
        function [BL_range, ERP, data_filtros, win, dist, zapato] = initVars()

            BL_range = [-200 0];
            ERP = []; 
            data_filtros = [];

            win = [];

            % Calculado a partir de los topoplots
            win.P200r.time  = [200 300];
            win.P200r.elect = [36:45];

            win.P200l.time  = [200 300];
            win.P200l.elect = [7:12 125:128];

            win.N400.time   = [300 450];
            win.N400.elect  = [1:6 17:21 30:35 50:53  65 66 87 97 98 110:114 124];
        %     win.N400.elect  = [1:6 18:20 31:35 50:53  66 87 97 98 110:114 124];


            win.P600.time   = [500 650];
        %     win.P600.time   = [600 700];
            win.P600.elect  = [3:9 16:22 29:32 34:38 44 45 47:51 112 113 121:126];
            win.P600.elect  = [3:9 16:22 29:32 35:38 44 45 125:126];


            dist.pred.preMJ_p = [];
            dist.pred.posMJ_p = [];
            dist.freq.preMJ_p = [];
            dist.freq.posMJ_p = [];
            dist.pos.preMJ_p  = [];
            dist.pos.posMJ_p  = [];
            dist.pred.preMJ_r = [];
            dist.pred.posMJ_r = [];
            dist.freq.preMJ_r = [];
            dist.freq.posMJ_r = [];
            dist.pos.preMJ_r  = [];
            dist.pos.posMJ_r  = [];
            
            zapato.data  = [];
            zapato.times = [];
            zapato.info  = [];
        end
        function [gramFilt, voltFilt, durFilt, filtData, indbadepoch, EEG] = filter(EEG, su)

            % Remove epochs by:
               % Grammatical: 
               %        - non-content words
               %        - < 2 chars
               %        - 1st word in sentence
               % Voltage:
               %        - >60mV in more than 5 electrodes
               % Fixation duration:
               %        - <50ms 

            filtData = [];   
           
            % Grammatical: 
            catsOpen = ['ANV' 'RS']; 
            m=min([EEG.epoch.pred]);

            gramFilt = (~ismember([EEG.epoch.catgram], catsOpen)) | ...
                                [EEG.epoch.pos] == 0 | ... % arranca en 0
                                [EEG.epoch.pos] == [EEG.epoch.SntcLength]-1 | ... % arranca en 0
                                [EEG.epoch.length] <=  2 ;%| ...
                                [EEG.epoch.pred] <= (m + abs(m*.05));

            % Voltage
            volt  = abs(EEG.data) > 80;   % Miro qué valores superan los 60uV
            sum_t = sum(volt, 2);         % Sumo en la dimensión tiempo
            sum_e = sum(sum_t > 0, 1);    % Sumo en la dimensión electrodos
            voltFilt = (squeeze(sum_e) > 10)'; % Trials con +5 electrodos

            % Fixation Duration
            for i = 1:length(EEG.epoch)
                ind = find([EEG.epoch(i).eventlatency{:}] == 0,1);
                EEG.epoch(i).dur = EEG.epoch(i).eventduration{ind};
            end
            %durFilt = [EEG.epoch.fixDur] < 200 | [EEG.epoch.fixDur] > 500;
            durFilt = [EEG.epoch.fixDur] < 50 | [EEG.epoch.fixDur] > 600 |...
                      [EEG.epoch.saccDur] < 5 | [EEG.epoch.saccDur] > 40;

            indbadepoch = voltFilt | gramFilt | durFilt;

            filtData.n_voltaje(su)    = sum(voltFilt & ~gramFilt);
            filtData.n_gramatical(su) = sum(gramFilt);
            filtData.n_total(su)      = sum(indbadepoch);

            % Éste índice lo uso para incluir, me quedo con el complemento
            indbadepoch = ~indbadepoch;
            
            for i = 1:length(EEG.epoch)
                EEG.epoch(i).filter = indbadepoch(i);
            end
        end
        function ERP = generateErp(EEG, DATA, ERP, su, indbadepoch, win, campo, voltFilt)
            
            %indbadepoch = logical(ones(1,length(EEG.epoch)));
            if strcmpi(campo, 'all')
                nombres_vars = {'all'};
                temp = [EEG.epoch.pos];
            elseif strcmpi(campo, 'pred')
                nombres_vars = {'q1','q2','q3'};
                temp = [EEG.epoch.pred];
                limites = quantile(temp(indbadepoch), [0 0.33 0.67 1]); 
                limites(1) = limites(1)-1;
            elseif strcmpi(campo, 'predPrev')
                nombres_vars = {'q1','q2','q3'};
                temp = [EEG.epoch.predPrev];
                limites = quantile(temp(indbadepoch), [0 0.33 0.67 1]);                
                limites(1) = limites(1)-1;
            elseif strcmpi(campo, 'predNext')
                nombres_vars = {'q1','q2','q3'};
                temp = [EEG.epoch.predNext];
                limites = quantile(temp(indbadepoch), [0 0.33 0.67 1]);                
                limites(1) = limites(1)-1;
            elseif strcmpi(campo, 'fixDur')
                nombres_vars = {'Short','Mid','Long'};
                temp = [EEG.epoch.fixDur];
                limites = [50 175 250 600]; 
            elseif strcmpi(campo, 'prevFixDur')
                nombres_vars = {'Short','Mid','Long'};
                temp = [EEG.epoch.prevFixDur];
                limites = [0 200 400 600];  
            elseif strcmpi(campo, 'saccDur')
                nombres_vars = {'Short','Mid','Long'};
                temp = [EEG.epoch.saccDur];
                limites = [0 16 32 60];               
            elseif strcmpi(campo, 'pred_type')
                nombres_vars = {'LoProv', 'HiProv', 'LoCommon','HiCommon'};
                tempPred = [EEG.epoch.pred];
                tempTipo = [EEG.epoch.sntType];
                limites = quantile(tempPred(indbadepoch), [0 0.48 0.52 1]);
                limites = quantile(tempPred(indbadepoch), [0 0.33 0.67 1]);                
                limites(1) = limites(1)-1;
                limites01 = [limites(1:2); zeros(1,2)];
                limites02 = [limites(3:4); zeros(1,2)];
                limites11 = [limites(1:2); ones(1,2) ];
                limites12 = [limites(3:4); ones(1,2) ];
                limites = cat(3, limites01, limites02, limites11, limites12);
            elseif strcmpi(campo, 'pred_type10')
                nombres_vars = {'Prov0','Prov1','Prov2','Prov3','Prov4',...
                                'Prov5', 'Prov6','Prov7','Prov8','Prov9',...
                                'Comm0','Comm1','Comm2','Comm3','Comm4',...
                                'Comm5', 'Comm6','Comm7','Comm8','Comm9'};
                nombres_vars = {'Prov0','Prov1','Prov2','Prov3','Prov4',...
                'Comm0','Comm1','Comm2','Comm3','Comm4'};

                tempPred = [EEG.epoch.pred];
                tempTipo = [EEG.epoch.sntType];
                limites = quantile(tempPred(indbadepoch), [0:.2:1]);                
                limites(1) = limites(1)-1;
                limites00 = [limites(1:2); zeros(1,2)];
                limites01 = [limites(2:3); zeros(1,2)];
                limites02 = [limites(3:4); zeros(1,2)];
                limites03 = [limites(4:5); zeros(1,2)];
                limites04 = [limites(5:6); zeros(1,2)];
%                 limites05 = [limites(6:7); zeros(1,2)];
%                 limites06 = [limites(7:8); zeros(1,2)];
%                 limites07 = [limites(8:9); zeros(1,2)];
%                 limites08 = [limites(9:10); zeros(1,2)];
%                 limites09 = [limites(10:11); zeros(1,2)];

                limites10 = [limites(1:2); ones(1,2)];
                limites11 = [limites(2:3); ones(1,2)];
                limites12 = [limites(3:4); ones(1,2)];
                limites13 = [limites(4:5); ones(1,2)];
                limites14 = [limites(5:6); ones(1,2)];
%                 limites15 = [limites(6:7); ones(1,2)];
%                 limites16 = [limites(7:8); ones(1,2)];
%                 limites17 = [limites(8:9); ones(1,2)];
%                 limites18 = [limites(9:10); ones(1,2)];
%                 limites19 = [limites(10:11); ones(1,2)];

                limites = cat(3, ...
                    limites00, limites01, limites02, limites03, limites04,...
                    limites10, limites11, limites12, limites13, limites14);
%                     limites05, limites06, limites07, limites08, limites09,
%                     limites15, limites16, limites17, limites18, limites19);

            elseif strcmpi(campo, 'sntc')
                nombres_vars = {'Memory','Common'};
                temp = [EEG.epoch.sntType] ;
                limites = [0 1]; 
            elseif strcmpi(campo, 'fixRank')    
%                 nombres_vars = {'fix1_m','fix2_m','fix3_m','fix4_m','fix5_m','fix6_m','fix7_m','fix8_m','fix9_m','fix10_m',...
%                                 'fix1_c','fix2_c','fix3_c','fix4_c','fix5_c','fix6_c','fix7_c','fix8_c','fix9_c','fix10_c'};
                nombres_vars = {'fix1_m','fix3_m','fix5_m','fix7_m','fix9_m',...
                                'fix1_c','fix3_c','fix5_c','fix7_c','fix9_c'};
                
                tempTipo = [EEG.epoch.sntType];
                tempFix = [EEG.epoch.pos];
                fix_prov = tempFix .* ~tempTipo;
                fix_comm = tempFix .*  tempTipo;
                tempPos = fix_prov + fix_comm;

                limites0 = [[1:2:10]' zeros(5,1)];
                limites1 = [[1:2:10]' ones(5,1)];
                limites = cat(1,limites0,limites1);
            elseif strcmpi(campo, 'PosAbs')
                nombres_vars = {'pos2','pos3','pos4','pos5','pos6',...
                                'pos7','pos8','pos9','pos10'};
                temp = [EEG.epoch.pos];
                limites = 1:10; 
            elseif strcmpi(campo, 'PosAbs_tipo')
                nombres_vars = {'pos2_m','pos3_m','pos4_m','pos5_m','pos6_m','pos7_m','pos8_m','pos9_m','pos10_m',...
                                'pos2_c','pos3_c','pos4_c','pos5_c','pos6_c','pos7_c','pos8_c','pos9_c','pos10_c'};
                
                tempTipo = [EEG.epoch.sntType];
                tempPos = [EEG.epoch.pos];
                
                nombres_vars = {'pos23_m','pos4_m','pos5_m','pos6_m','pos78_m','pos910_m',...
                                'pos23_c','pos4_c','pos5_c','pos6_c','pos78_c','pos910_c'};
                tempPos(tempPos==3)  = 2;
                tempPos(tempPos==8)  = 7;
                tempPos(tempPos==10) = 9;
                
                posAbs_prov = tempPos .* ~tempTipo;
                posAbs_comm = tempPos .*  tempTipo;
                tempPos = posAbs_prov + posAbs_comm;

                uqPos = unique(tempPos);
                uqPos = uqPos(uqPos>=2 & uqPos <=9);
                limites0 = [uqPos' zeros(length(uqPos),1)];
                limites1 = [uqPos' ones(length(uqPos),1)];
               
                limites = cat(1,limites0,limites1);
            elseif strcmpi(campo, 'PosRel')
                nombres_vars = {'pos_m2','pos_m1', ...
                                'pos_0', ...
                                'pos_1', 'pos_2'};
                temp = [EEG.epoch.posRelRP] + 100*[EEG.epoch.sntType];
                limites = -2:2; 
            elseif strcmpi(campo, 'PosRel_tipo')
                nombres_vars = {'pos_m4_m','pos_m3_m','pos_m2_m','pos_m1_m', 'pos_0_m', 'pos_1_m', 'pos_2_m', 'pos_3_m', 'pos_4_m','pos_5_m',...
                                'pos_m4_c','pos_m3_c','pos_m2_c','pos_m1_c', 'pos_0_c', 'pos_1_c', 'pos_2_c', 'pos_3_c', 'pos_4_c','pos_5_c'};
                
                tempTipo = [EEG.epoch.sntType];
                posRel_prov = [EEG.epoch.posRelRP] .* ~tempTipo;
                posRel_comm = [EEG.epoch.posRelRP] .*  tempTipo;
                tempPos = posRel_prov + posRel_comm ;
                                
                limites0 = [[-4:5]' zeros(10,1)];
                limites1 = [[1:10]' ones(10,1)];
                limites = cat(1,limites0,limites1);
            elseif strcmpi(campo, 'pre_post')
                nombres_vars = {'memory_pre','memory_post',...
                                'common_pre','common_post'};
                type     = [EEG.epoch.sntType];
                posRelRP = [EEG.epoch.posRelRP];
                ind_memory_pre  = type == 0 & posRelRP <= 0;
                ind_memory_post = type == 0 & posRelRP >  0;
                ind_common_pre  = type == 1 & posRelRP <= 4;
                ind_common_post = type == 1 & posRelRP >  4;
            else
                fprintf('Pifie %s \n', campo)
            end

            for c = 1:length(nombres_vars)
                variable = nombres_vars{c};
                if strcmpi(campo, 'all')
                    indtrial = find(temp>-1 & indbadepoch);
                elseif strcmpi(campo, 'pred') || ...
                   strcmpi(campo, 'predPrev') || ...
                   strcmpi(campo, 'predNext') 
                    if limites(c) == limites(c+1)
                        keyboard
                        indtrial = find(temp >= limites(c) & temp <= limites(c+1) & indbadepoch);
                    else
                        indtrial = find(temp >  limites(c) & temp <= limites(c+1) & indbadepoch);
                    end
                elseif  strcmpi(campo, 'PosRel_tipo') || ...
                        strcmpi(campo, 'PosAbs_tipo')
                    cond1 = tempPos == limites(c,1);
                    cond2 = tempTipo == limites(c,2);
                    indtrial = find(cond1 & cond2);
                elseif  strcmpi(campo, 'fixRank')                
                    cond1 = tempPos == limites(c,1) | tempPos == limites(c,1)+1;
                    cond2 = tempTipo == limites(c,2);
                    indtrial = find(cond1 & cond2);
                elseif strcmpi(campo, 'pred_type') || ...
                       strcmpi(campo, 'pred_type10')
                    cond1 = tempPred > limites(1,1,c) & tempPred <= limites(1,2,c);
                    cond2 = tempTipo == limites(2,1,c);
                    cond3 = indbadepoch;
                    indtrial = find(cond1 & cond2 & cond3);
                elseif strcmpi(campo, 'prevFixDur') || ...
                       strcmpi(campo, 'fixDur') || ...
                       strcmpi(campo, 'saccDur')
                    indtrial = find(temp >= limites(c) & temp <= limites(c+1) & indbadepoch);
                elseif strcmpi(campo, 'PosAbs') || ...
                       strcmpi(campo, 'PosRel') || ...
                       strcmpi(campo, 'sntc')   
                    indtrial = find(temp == limites(c) & indbadepoch);
                elseif strcmpi(campo, 'pre_post')
                    indtrial = eval(['ind_' variable]);
                end 

                ERP(su).(campo).(variable).m = nanmean(EEG.data(:,:,indtrial),3);
                ERP(su).(campo).(variable).n = size(EEG.data(:,:,indtrial),3); 
                ERP(su).(campo).(variable).trials = EEG.data(:,:,indtrial); 


                roi_names = fieldnames(win);
                for r = 1:length(roi_names)
                    roi = roi_names{r};
                    indelec = win.(roi).elect;
                    indtimes = (EEG.times >= win.(roi).time(1)) & (EEG.times <= win.(roi).time(2)) ;

                    ERP(su).([campo '_ROI' roi]).(variable).m = nanmean(nanmean(EEG.data(indelec,:,indtrial),1),3);
                    ERP(su).([campo '_ROI' roi]).(variable).n = size(EEG.data(indelec,:,indtrial),3);

                    ERP(su).([campo '_ROI_time' roi]).(variable).m = nanmean(nanmean(nanmean(EEG.data(indelec,indtimes,indtrial),1),2),3);
                    ERP(su).([campo '_ROI_time' roi]).(variable).n = size(EEG.data(indelec,indtimes,indtrial),3);
                end
            end
        end
        function zapato = generateZapato(EEG, elect, filtro, zapato)
            times = EEG.times>-100;
            data = squeeze(EEG.data(elect,times,filtro))'; 
            
            zapato.data  = [zapato.data; data];
            zapato.times = EEG.times(times);
            zapato.info  = [zapato.info EEG.epoch(filtro)];

        end
        function EEG = addIndProv(EEG)
            DATA = EEG.info.DATA;
            for i = 1:length(DATA)
                DATA(i).indProv = DATA(i).ind - (94 * DATA(i).tipo);
                
                indProv = find(DATA(i).indProv == [DATA.ind] & [DATA.tipo] == 0);
                if ~isempty(indProv);
                    DATA(i).posRelRP = DATA(indProv).posRelRP;
                else 
                    DATA(i).posRelRP = [];
                end
            end
            EEG.info.DATA = DATA;            
        end
    end%methods
end% classdef

function [line, lines] = findText(C, texts, first)
    % find one of the texts in file from line 1 to line maxlines
    % output: # of lines and text of line

    line  = [];
    lines = [];
    for i = 1:length(texts)% find each of str in texts
        if first == 1
            counter = 0;
%             keyboard
            while isempty(line) && counter < length(C)
                counter = counter+1;
                tline = C{counter};    
                if ~isempty(strfind(tline,texts{i})) 
                    line = [line counter];
                    lines = [lines tline];
                    return 
                end
                %    disp([num2str(contador) tline])
            end
        else
            a = cellfun(@(x) ~isempty(strfind(x,texts{i})), C);
            index = find(a);
            line  = [line index];
            lines = [lines C(line)];
        end                

    end
%     if isempty(line)
%         tline = [];
%     end
end
function all   = findBestEye(cal)

    fprintf('Finding best eye:\n')
    [~, lines] = findText(cal, {'!CAL VALIDATION'},0);
    
    splited = strsplit(lines{end});
    L = any(findText(splited, {'LEFT'},0));
    R = any(findText(splited, {'RIGHT'},0));
    
    if xor(L, R) 
        eyes = {'L' 'R'};
        all.eye     = eyes{[L R]};
        all.bestCal = all.eye;

        ind = find(cellfun(@(x) ~isempty(x), (strfind(splited,'ERROR'))));
        all.calErr  = str2num(splited{ind+1});
        
        fprintf('\t%s: %.2f\n',  all.eye,...
                                 all.calErr)
        fprintf('\tBest eye: %s\n', all.bestCal)

    elseif L && R
        all.eye = 'BOTH';
        
        indR = find(cellfun(@(x) ~isempty(x), (strfind(splited,'ERROR'))));
        abortedR = any(cellfun(@(x) ~isempty(x), (strfind(splited,'ABORTED'))));
        if abortedR; calR = inf; else; calR = str2num(splited{indR+1}); end

        indL = find(cellfun(@(x) ~isempty(x), (strfind(splited,'ERROR'))));
        abortedL = any(cellfun(@(x) ~isempty(x), (strfind(splited,'ABORTED'))));
        if abortedL; calL = inf; else; calL = str2num(splited{indL+1}); end

        if calR < calL
            all.bestCal = 'R';
            all.calErr  = calR;
        else
            all.bestCal = 'L';
            all.calErr  = calL;
        end


        fprintf('\tRIGHT: %.2f  LEFT: %.2f\n',  calR, calL)
        fprintf('\tBest eye: %s\n', all.bestCal)

    else 
        keyboard
        fprintf('Unknwon eye\n')
    end
end
function ppc   = calculatePPC(DATA)
    % length of sentence image / length (in chars) of sentence
    % mean of all the sentences
    tmp = [];
    for i = 1:length(DATA)
        tmp = [tmp size(DATA(1).imagen,2)/length(DATA(1).oracion)];
    end

    ppc = mean(tmp); % Pixel per character
end

%         function fixDur = calcFixDur(tEpoch)
%             % input: EEG.epoch
%             
%             % Find event of fix onset
%             iOnF = find([tEpoch.eventlatency{:}] == 0, 1, 'first');
%             eye = tEpoch.eventtype{iOnF}(1);
%             tOnF = 0;
%             tOfF = 1000; % Default: in case fixation longer than epoch, set on 500ms
% 
%             % Find fix offset (by finding first following saccade or blink)
%             iOfF = strcmp({tEpoch.eventtype{iOnF+1:end}}, 'R_saccade') + ...
%                    strcmp({tEpoch.eventtype{iOnF+1:end}}, 'L_saccade') + ...
%                    strcmp({tEpoch.eventtype{iOnF+1:end}}, 'R_blink') + ...
%                    strcmp({tEpoch.eventtype{iOnF+1:end}}, 'L_blink');
% 
%             if any(iOfF)
%                 iOfF = find(iOfF,1,'first') + iOnF;
%                 tOfF = tEpoch.eventlatency{iOfF}; 
%             end   
% 
%             % Find Saccade onset index
%             iOnS = strcmp({tEpoch.eventtype{1:iOnF}}, 'R_saccade') + ...
%                    strcmp({tEpoch.eventtype{1:iOnF}}, 'L_saccade');
%             iOnS = find(iOnS, 1, 'last');
% 
%             % Skip Fixations under 50ms
%             fixDur = (tOfF - tOnF );
%             
%         end
%         function cals  = findCalibratios(messages)
%             fprintf('Finding calibrations marks\n')
%             
%             C = messages;
%             [stLinesNum, ~]           = findText(C, {'CALIBRATION'}, 0);
%             [modeLinesNum, modeLines] = findText(C, {'!MODE'}, 0);
% 
%             headerLim  = [];
%             % For each block find the last calibration
%             for i = 1:length(modeLinesNum)
%                 mode     = strsplit(strtrim(modeLines{i}), ' '); 
%                 mode     = mode(length(mode));
%                 previous = find((stLinesNum  - modeLinesNum(i)) < 0);
% 
%                 if strcmpi(mode, 'L') || strcmpi(mode, 'R')
%                     % If MONO, I only care for the last calibration
%                     calStart = stLinesNum(length(previous));
%                 elseif strcmpi(mode, 'LR')
%                     % If BINO, I need the last 2 calibrations (L and R)
%                     calStart = stLinesNum(length(previous)-1);
%                 end
% 
%                 headerLim  = [headerLim  ; [calStart modeLinesNum(i)]];
%             end
% 
%             % if there was a fail start there is 2 calibrations for 1 block
%             diffStarts = diff(headerLim(:,1));
%             index      = find(diffStarts == 0);
%             headerLim(index,:)  = [];
%             
%             for i = 1:size(headerLim,1)            
%                 tCal = C(headerLim(i,1) : headerLim(i,2));
%                 cals(i) = findBestEye(tCal);
%             end
%         end
%         function DATA  = addFieldsEstim(DATA,ESTIM,campos)
%             try
%                 fprintf('\tAdd fields from ESTIM to DATA\n')
%                 for i=1:length(DATA)
%                     for indcampo=1:length(campos)
%                         DATA(i).(campos{indcampo})=ESTIM(DATA(i).ind+1).(campos{indcampo});
%                     end
%                 end
%             catch ME
%                 ME
%                 keyboard
%             end
%         end
%         function DATA  = addSpaces(DATA)
%             % Pixel for character
%             ppc = calculatePPC(DATA);
%             for trial=1:length(DATA)
%                 spaces = strfind(DATA(trial).oracion,' ');
%                 spaces = [0 spaces length(DATA(trial).oracion)+1];        
%                 if strcmp(DATA(trial).oracion(end),' ')
%                     spaces(end) = [];
%                 end
%                 str = DATA(trial).posoracion(1);
%                 DATA(trial).spaces = str +(spaces-.5)*ppc;
%             end
%             
% %             Max Vertical Pos for fixations
% %             mvp = init.height/2 + init.height*.2;
% %             all = ET.parsed(block);
% %             all = interpolateBlinksSaccades(all);            
% %             fprintf('\tBest eye: %s\n\n', all.bestCal)                      
% %             DATA = addFix(all,DATA,mvp, ppc);
% %             DATA = assignFixToWords(DATA, init);
% %             DATA = assignFixToChar(DATA, ppc);
% %             DATA = assignFixToFirstChar(all, DATA);
% %             DATA = addFixCross(all, DATA, init);
% %             DATA = addPupil(all, DATA, SUJ(iSuj));
% %             DATA = addMsgsToTrial(all, DATA);    
% %             DATA = addRT(all, DATA);
% %             DATA = addSpaces(DATA, ppc);
%                           
%         end
%         function DATA  = addLatency(DATA, EEG)
%             starts = arrayfun(@(x) ~isempty(x{1}),regexp({EEG.event.type}, '220')); 
%             starts = [EEG.event(starts).latency]/EEG.srate;
% 
%             ends = arrayfun(@(x) ~isempty(x{1}),regexp({EEG.event.type}, '221')); 
%             ends = [EEG.event(ends).latency]/EEG.srate;
%             
%             % If recording stated late and the first mark is an end mark
%             if length(starts)+1 == length(ends)
%                 ends = ends(2:end);
%             end
%             
%             inds  = [starts(4:end)'  ends(4:end)']; % Delete 3 sentences of practice 
%             for iData = 1:length(DATA)
%                 DATA(iData).latency = inds(iData,:);
%             end
%         end
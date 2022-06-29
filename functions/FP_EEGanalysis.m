classdef FP_EEGanalysis
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
       
        % Analysis
        function ERP = generateErp(EEG, ERP, su, campo)
            
            indbadepoch = logical(ones(1,length(EEG.epoch)));
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
            elseif strcmpi(campo, 'sntc')
                nombres_vars = {'Memory','Common'};
                temp = [EEG.epoch.sntType] ;
                limites = [0 1]; 
            elseif strcmpi(campo, 'fixRank')    
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
                elseif strcmpi(campo, 'pred_type') 
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

            end
        end
        function erpFig = dimigenFig(EEG, elects, erpFig)
            
            
            if ~isfield(erpFig,'im')
                erpFig.im= []; 
                erpFig.pfixDur=[]; 
                erpFig.fixDur=[]; 
                erpFig.nfixDur=[];
            end
            tmp = [EEG.epoch.prevFixDur];
            erpFig.pfixDur = [erpFig.pfixDur tmp];

            tmp = [EEG.epoch.fixDur];
            erpFig.fixDur = [erpFig.fixDur tmp];

            tmp = [EEG.epoch.nextFixDur];
            erpFig.nfixDur = [erpFig.nfixDur tmp];

            tmp = squeeze(nanmean(EEG.data(elects,:,:),1));
            erpFig.im = cat(2, erpFig.im, tmp);
        end
    end
end


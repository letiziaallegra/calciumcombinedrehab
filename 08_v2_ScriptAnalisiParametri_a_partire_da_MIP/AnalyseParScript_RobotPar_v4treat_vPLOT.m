%load dataMouseGCamp_RobotPar_Median ... -> MeanPar

if ~exist('MeanPar')
    error('LOAD the matrix "dataMouseGCamp_RobotPar_Median....mat" ')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%select what you want to plot
ParNameToPlot = {'Trial','Target Time','Force','Sub Mov','Attempts','Max Force'};
ROINameToPlot = {'Robot Par'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%select if plotting all of the animals or their mean (stroke vs control vs rehab)
% MP = 1; % ->stroke vs control vs rehab
MP = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%select if you want to pool the data of the one week
POOL = 1; % 
% POOL = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumAnimals = 1;
for iL=1:size(MeanPar,2)
    NumAnimals_cell = size(MeanPar{3,iL},1);
    if NumAnimals_cell >= NumAnimals
        
        NumAnimals_buf = sum(~cellfun('isempty',MeanPar{3,iL}(:,2)));
               
        if NumAnimals_buf > NumAnimals  %not empty
            indNumDay  = iL;
            AnimalNamesList = MeanPar{3,iL}( ~cellfun('isempty',MeanPar{3,iL}(:,2)),2);
            NumAnimals      = NumAnimals_buf;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Days to take
ListDaysToTake =  cell2mat(MeanPar(1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumPar  = size(MeanPar,1);
NumDay  = size(MeanPar,2);

MeanForceTot   = [];
StdForceTot    = [];
StdErrForceTot = [];

MeanForceWeekTot = [];
stdForceWeekTot  = [];
stdERRForceWeekTot =[];

%for the 4 parameters from Robot
Fig_all_days = figure('Name','Time course motor performance on the platfom');
Fig_weeks    = figure('Name','Time course weekly motor performance on the platfom');

for pi=1:length(ParNameToPlot)-1
     
    ParForce = pi+1; 
    ROI = 1;
    
    MeanForceTot = [];
    StdForceTot = [];
     
    for na=1:NumAnimals
        
        %animal name
        AnNam = AnimalNamesList{na};
        
        MeanForce = [];
        StdForce  = [];
        MeanForceTot_Wk   = [];
        StdErrForceTot_Wk = [];
        
        
        for nd=1:NumDay
            
            na_inDay = find(strcmp(MeanPar{ParForce,nd}(:,2),AnNam));
            
            
            if ~isempty(na_inDay)
                
                NameMouse  = AnNam;
                TreatMouse = MeanPar{ParForce,nd}{na_inDay,4};
                
                if sum(ismember(ListDaysToTake,nd))>0
                    
                    if ~isempty(MeanPar{ParForce,nd}{na_inDay,1})
                        MeanForce = [MeanForce; MeanPar{ParForce,nd}{na_inDay,1}(ROI,1)];
                        StdForce  = [StdForce;  MeanPar{ParForce,nd}{na_inDay,1}(ROI,2)];
                        

                    else
                        MeanForce = [MeanForce; NaN];
                        StdForce  = [StdForce;  NaN];
                        
                    end
                else
                    display('NO');
                end
            else
                %no animal in that day
                MeanForce = [MeanForce; NaN];
                StdForce  = [StdForce;  NaN];
 
            end
            
        end
        
        %tutti i valori medi e std del parametro nei vari giorni per quell'animale na
        %(diventano colonne delle matrici _Tot)
        MeanForceTot = [MeanForceTot MeanForce];
        
        NameMouseTot{na,1}   = NameMouse;
        TreatMouseTot(na,1)  = TreatMouse;
        
        
        
        if na==NumAnimals
            
            if MP
                
                MeanForce_Buf = MeanForceTot;
                
                clear MeanForceTot MeanFluoTot StdForceTot StdFluoTot
                
                for tr=0:3
                    
                    %nella matrice MeanForce_Buf prendo solo le colonne che si riferiscono a un certo treatment
                    StdForceTot(:,tr+1)     = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2);
                    StdErrForceTot(:,tr+1)  = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2)/sqrt(sum(TreatMouseTot==tr));
                    %                 MeanForceTot(:,tr+1) = nanmean(MeanForce_Buf(:,TreatMouseTot==tr),2);
                    MeanForceTot(:,tr+1) = nanmedian(MeanForce_Buf(:,TreatMouseTot==tr),2);
                                                             
                end
                
                Filename = {'control','stroke','rehab_Bont','rehab'};
                
            else
                MeanForce_Buf = MeanForceTot;
                Filename      = NameMouseTot;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Plot
            
            figure(Fig_all_days)
            
            if length(ListDaysToTake) ~= length(NumDay)
                t = ListDaysToTake';
            else
                t = [1:NumDay]';
            end
            
            subplot(2,3,pi)
            colorMP = [0 0 1;1 0 0;0 1 0; 1 0 1];
            if MP
                hold on                
                for mp=1:3
                    h_all = errorbar(t,MeanForceTot(:,mp),StdErrForceTot(:,mp),'Color',colorMP(mp,:),'LineWidth',2);
                    set(gca,'FontSize',14);
                end
            else
                %             plot(t,MeanForceTot,'-..','MarkerSize',15);
                for it=1:length(TreatMouseTot)
                    hold on
                    if TreatMouseTot(it)==0
                        h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(1,:));
                    elseif TreatMouseTot(it)==1
                        h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(2,:));
                    elseif TreatMouseTot(it)==2
                        h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(3,:));
                    elseif TreatMouseTot(it)==3
                        h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(4,:));
                    end
                    set(gca,'FontSize',14)
                end
            end
%             title([ParNameToPlot{ParForce},' ',ROINameToPlot{ROI}])
            ylabel(ParNameToPlot{ParForce})
            xlabel('Day')
            xlim([t(1)-1 t(end)+1]);
            set(gca,'xTick',t)
            set(gca,'xTickLabel',(num2str(t)))
            
            if   pi == 5
                legend(Filename,'Location','EastOutside')
            end
                        
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Plot settimanale  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Plot settimanale (tutte e 4 settimane) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(MeanForceTot,1)>5
                
                weekSpace = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20];
                weekNum   = [1:size(weekSpace,1)];
                
                %Plot settimanale (ultima settimana) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                weekSpace = [1 2 3 4 5];
                weekNum   = [1:size(weekSpace,1)];
                
            end
            
            figure(Fig_weeks)
                      
            
            hold on
            colorMP = [0 0 1;1 0 0;0 1 0; 1 0 1];
            
            %%%%initialization
            presTreat = unique(TreatMouseTot);
            lenTreat  = length(presTreat);
            medTot_ForceTot_Treat_Wk    = zeros(lenTreat,length(weekNum));
            stdTot_ForceTot_Treat_Wk    = zeros(lenTreat,length(weekNum));
            errstdTot_ForceTot_Treat_Wk = zeros(lenTreat,length(weekNum));
            %%%%%
            
            %treatment
            for tr=1:lenTreat
                
                %to save excel Stat
                ExcelSaveForce = [];
                
                %prendo solo le colonne che si riferiscono al treatment
                MeanForceTot_Treat = MeanForce_Buf(:,TreatMouseTot==presTreat(tr));
                
                for iw=1:length(weekNum)
                    
                    % % % % %
                    %media sui giorni (un valore per animale)
                    % % % % %
                    
                    %prendo solo le righe che si riferiscono alla week
                    MeanForceTot_Treat_Wk =  MeanForceTot_Treat(weekSpace(iw,:),:);
                    
                    %%%%median value for each animal in this week
                    if POOL == 1
                        med_ForceTot_Treat_Wk       =  reshape(MeanForceTot_Treat_Wk,size(MeanForceTot_Treat_Wk,1)*size(MeanForceTot_Treat_Wk,2),1);
                    else
                        %med_ForceTot_Treat_Wk       =  nanmedian(MeanForceTot_Treat_Wk,1)';
                        med_ForceTot_Treat_Wk       =  nanmean(MeanForceTot_Treat_Wk,1)';
                    end
                    
                    %total median value for the animals in this week
%                     medTot_ForceTot_Treat_Wk(tr,iw)       =  nanmedian(med_ForceTot_Treat_Wk);
                    medTot_ForceTot_Treat_Wk(tr,iw)       =  nanmean(med_ForceTot_Treat_Wk);
                    
                    stdTot_ForceTot_Treat_Wk(tr,iw)       =  nanstd(med_ForceTot_Treat_Wk,[]);
                    errstdTot_ForceTot_Treat_Wk(tr,iw)    =  stdTot_ForceTot_Treat_Wk(tr,iw) /sqrt(length(med_ForceTot_Treat_Wk));
                    
                    %to save excel Stat
                    ExcelSaveForce = [ExcelSaveForce, med_ForceTot_Treat_Wk];
                        
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%   Plot    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %force
                hold on
                subplot(2,3,pi)
                h_week = barwitherr(errstdTot_ForceTot_Treat_Wk , medTot_ForceTot_Treat_Wk);
                set(gca,'FontSize',14)
                title([ParNameToPlot{ParForce}])
                ylabel('Force Par')
                if weekNum == 1
                    set(gca,'xTick',[1 2 3])
                    set(gca,'xTickLabel',{'ctrl','strk','rehab'})
                    xlabel('Last Week')
                else
                    set(gca,'xTick',weekNum)
                    set(gca,'xTickLabel',(num2cell(weekNum )))
                    xlabel('Week')
                end
%                 legend(Filename)
                
                
                
                %to save excel Stat
                if size(weekSpace,1)>1
                    ExcelStr = 'Rehab_only';
                else
                    ExcelStr = 'Ctrl_Stk_Rehab';
                end
                
                %%%SAVE EXCEL
%                 xlswrite(['DataRobotStat_',ExcelStr,'_',ParNameToPlot{ParForce},'_'],ExcelSaveForce ,tr)
                %
                               
            end
            
            %%%SAVE FIGURES
%             saveas(h_all,'Fig_Robot_Parameter','fig');
%             saveas(h_all,'Fig_Robot_Parameter','jpg');
%             saveas(h_week,'Fig_Weekly_Robot_Parameter','fig');
%             saveas(h_week,'Fig_Weekly_Robot_Parameter','jpg');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end
        
    end
end






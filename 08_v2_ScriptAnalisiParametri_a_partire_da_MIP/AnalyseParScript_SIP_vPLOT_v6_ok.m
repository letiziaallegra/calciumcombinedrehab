%load dataMouseGCamp_Median -> MeanPar
clear all
close all
clc

PAR_INT = [5 6 11];
%PAR_INT = [3 4 7 8 9 10 11 12 13 14];

%salva excel con dati per la statistica
XLS_SAVE = 0;


for indexipar=1:length(PAR_INT)
    %     for ipar=13:13
    
    ipar = PAR_INT(indexipar)
    
    load('dataMouseGCamp_MEAN_SelROI_1_13-Mar-2018_10-44.mat')
    
%close all

if ~exist('MeanPar')
    error('LOAD the matrix "dataMouseGCamp_Median_Sel....mat" ')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SelROI=1;

if SelROI == 0
    %%%select what you want to plot -> MANUAL ROIs
    ParNameToPlot = {'Trial','Status','Diff Start Time','Duration','Max','Diff Max Time','FWHM','AUC','PtPAmp','nSubPeaks'};
    ROINameToPlot = {'Whole'};
    %% ROI of interest
    ROI = 1;
    %%%
elseif SelROI == 1
    %%%select what you want to plot -> MIP Based ROIs
    ParNameToPlot = {'Trial','Status','OnsetT','Duration','Max','MaxT','FWHM','AUC','PtPAmp','nRealPeaks','Slope','TimeToPeak','FirstMovT-OnsetT','FirstMovT-MaxT'};
    ROINameToPlot = {'MIP','SIP Rostral','SIP Caudal'};
    Folder_Where_Save = 'C:\Users\CNR-SSSUP\Google Drive_SL\Piattaforma Stefano\ELABORAZIONE DATA\08_v2_ScriptAnalisiParametri_a_partire_da_MIP\_Folder_Image_And_Stat_4wREHAB_STIM';
    %% ROI of interest
    ROI_List = [1 2 3];
    %%%
elseif SelROI == 2
    %%%select what you want to plot -> Centroid Based ROIs
    ParNameToPlot = {'Trial','Status','Diff Start Time','Duration','Max','Diff Max Time','FWHM','AUC','PtPAmp','nSubPeaks'};
    ROINameToPlot = {'Motor','Barrel','Occipital','Motor Small','Barrel Small','Occipital Small','Whole','Main ROI'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compare pair of different parameters
CompDiffPars = 0;
% CompDiffPars = 1;

if CompDiffPars == 0
    ParForce = ipar;
    ParFluo = ParForce;
else
    ParForce = 6;
    ParFluo  = 6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%select if plotting all of the animals or their mean (stroke vs control vs rehab)
MP = 1; % ->control vs stroke vs rehab 1w vs rehab 4w 
% MP = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%select if you want to pool the data of the one week
POOL = 1; % ->control vs stroke vs rehab 1w vs rehab 4w 
% POOL = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% number ID for each group
% control = 0; stroke = 1; robot+Bont_4w = 2; robot+Bont_1w = 3;
% robot_4w = 21; robot_1w = 31; robot+stim_4w=214; TOX =5
%treat_ID = [0 1 2 3 21 31 5];
treat_ID = [0 1 2 3 5];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NumAnimals = 1;
% for iL=1:size(MeanPar,2)
%     NumAnimals_cell = size(MeanPar{3,iL},1);
%     if NumAnimals_cell >= NumAnimals
%         
%         NumAnimals_buf = sum(~cellfun('isempty',MeanPar{3,iL}(:,2)));
%                
%         if NumAnimals_buf > NumAnimals  %not empty
%             indNumDay  = iL;
%             AnimalNamesList = MeanPar{3,iL}( ~cellfun('isempty',MeanPar{3,iL}(:,2)),2);
%             NumAnimals      = NumAnimals_buf;
%         end
%     end
% end


%cerca numero animali e nome animali (anche doppi)
for iL=1:size(MeanPar,2)
    
    if iL == 1
        AnimalNamesList = MeanPar{3,iL}(:,2);
        NumAnimals_cell = length(AnimalNamesList);
    else
        AnimalNamesList_Buf = MeanPar{3,iL}(:,2);
        NumAnimals_cell_Buf = length(AnimalNamesList_Buf);
        
        if NumAnimals_cell < NumAnimals_cell_Buf
            
        elseif NumAnimals_cell > NumAnimals_cell_Buf
            
        elseif NumAnimals_cell == NumAnimals_cell_Buf
            
            for inC=1:NumAnimals_cell    
                %sostituisci celle vuote
                if isempty(AnimalNamesList{inC})
                    AnimalNamesList{inC} = AnimalNamesList_Buf{inC};
                end

            end
            
            
        end
    end
    
    %ultimo day da check
    if iL == size(MeanPar,2)
%         AnimalNamesList
        NumAnimals = length(AnimalNamesList);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Days to take
ListDaysToTake =  cell2mat(MeanPar(1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumPar  = size(MeanPar,1);
NumDay  = size(MeanPar,2);

add=0;
%add=15;% se vuoi analizzare i giorni della quarta settimana

MeanForceTot   = [];
StdForceTot    = [];
StdErrForceTot = [];

MeanFluoTot    = [];
StdFluoTot     = [];
StdErrFluoTot  = [];

MeanForceWeekTot = [];
stdForceWeekTot  = [];
stdERRForceWeekTot =[];
MeanFluoWeekTot  = [];
stdFluoWeekTot   = [];
stdERRFluoWeekTot = [];



fig_Par      = figure('Name','Fig Parameters over the days of traing');
Fig_weeks    = figure('Name','Time course weekly Fluo and Force Parameters');



for nroi=1:length(ROI_List)
    
    ROI = ROI_List(nroi);
    
    %Load
    MeanForceTot = [];
    StdForceTot  = [];
    MeanFluoTot  = [];
    StdFluoTot   = [];
   
    %
    AnNam_Prev = [];
    
    %scorro animali
    for na=1:NumAnimals
        
        if na==NumAnimals-2
            ciao=1;
        end
        
        %animal name
        AnNam = AnimalNamesList{na};
        i_Check = 0;
        
        MeanForce = [];
        StdForce  = [];
        MeanFluo  = [];
        StdFluo   = [];
        

        %scorro i giorni per quell'animale
        for nd=1:NumDay
            
            na_inDay = find(strcmp(MeanPar{ParForce,nd}(:,2),AnNam));
            
            if ~isempty(na_inDay) && sum(ismember(na_inDay,na))
                
                NameMouse  = AnNam;

                %
                if length(na_inDay)>1 
                    if strcmp(AnNam,AnNam_Prev) & i_Check == 0
                        i_TR    = i_TR+1;
                        i_Check = 1;
                    elseif i_Check == 0
                        AnNam_Prev = AnNam;
                        i_TR    = 1;
                        i_Check = 1;
                    end
                else
                    i_TR = 1;
                end
                %
                
                TreatMouse = MeanPar{ParForce,nd}{na_inDay(i_TR),4};
                
                if sum(ismember(ListDaysToTake,nd+add))>0
                    
                    if ~isempty(MeanPar{ParForce,nd}{na_inDay(i_TR),1})
                        
                        MeanForce = [MeanForce; MeanPar{ParForce,nd}{na_inDay(i_TR),1}(ROI,1)];
                        StdForce  = [StdForce;  MeanPar{ParForce,nd}{na_inDay(i_TR),1}(ROI,2)];
                        MeanFluo  = [MeanFluo;  MeanPar{ParFluo,nd}{na_inDay(i_TR),1}(ROI,3)];
                        StdFluo   = [StdFluo;   MeanPar{ParFluo,nd}{na_inDay(i_TR),1}(ROI,4)];
                        
                    else
                        MeanForce = [MeanForce; NaN];
                        StdForce  = [StdForce;  NaN];
                        
                        MeanFluo  = [MeanFluo;  NaN];
                        StdFluo   = [StdFluo;   NaN];
                        
                    end
                else
                    display('NO');
                end
                
            else
                %no animal in that day
                MeanForce = [MeanForce; NaN];
                StdForce  = [StdForce;  NaN];
                
                MeanFluo  = [MeanFluo;  NaN];
                StdFluo   = [StdFluo;   NaN];
            end
            
        end
        
        
        %tutti i valori medi e std del parametro nei vari giorni per quell'animale na
        %(diventano colonne delle matrici _Tot)
        MeanForceTot = [MeanForceTot MeanForce];
        %     StdForceTot  = [StdForceTot  StdForce];
        
        MeanFluoTot = [MeanFluoTot MeanFluo];
        %     StdFluoTot  = [StdFluoTot  StdFluo];
        
        NameMouseTot{na,1}   = NameMouse;
        TreatMouseTot(na,1)  = TreatMouse;
        
        
        
        if na==NumAnimals
            
            if MP
                
                MeanForce_Buf = MeanForceTot;
                MeanFluo_Buf  = MeanFluoTot;
                
                clear MeanForceTot MeanFluoTot StdForceTot StdFluoTot
                
                for i_tr=1:length(treat_ID) 
                    
                    tr = treat_ID(i_tr);
                    
                    %nella matrice MeanForce_Buf prendo solo le colonne che si riferiscono a un certo treatment
                    StdForceTot(:,i_tr)     = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2);
                    StdErrForceTot(:,i_tr)  = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2)/sqrt(sum(TreatMouseTot==tr));
                    MeanForceTot(:,i_tr)    = nanmean(MeanForce_Buf(:,TreatMouseTot==tr),2);
%                     MeanForceTot(:,i_tr) = nanmedian(MeanForce_Buf(:,TreatMouseTot==tr),2);
                    
                    %nella matrice MeanFluo_Buf prendo solo le colonne che si riferiscono a un certo treatment
                    StdFluoTot(:,i_tr)     = nanstd(MeanFluo_Buf(:,TreatMouseTot==tr),[],2);
                    StdErrFluoTot(:,i_tr)  = nanstd(MeanFluo_Buf(:,TreatMouseTot==tr),[],2)/sqrt(sum(TreatMouseTot==tr));
                    MeanFluoTot(:,i_tr)    = nanmean(MeanFluo_Buf(:,TreatMouseTot==tr),2);
%                     MeanFluoTot(:,i_tr)  = nanmedian(MeanFluo_Buf(:,TreatMouseTot==tr),2);
                                                             
                end
                
                Filename = {'control','stroke','rehab 4w','rehab 1w','robot 4w','robot 1w','reSTI 4w'};
                Filename = {'control','stroke','rehab 4w','rehab 1w','only tox'};
                
            else
                Filename = NameMouseTot;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Plot
            
            
            %%%%%%%%
            if 1
                figure(fig_Par)
                
                if length(ListDaysToTake) ~= length(NumDay)
                    t = ListDaysToTake';
                else
                    t = [1:NumDay]';
                end
                %             colorMP = [0 0 1;1 0 0;0 1 0]*ROI/length(ROI_List);
                %             colorMP = [0 0 1;1 0 0;0 1 0]/length(ROI_List);
                colorMP = [0 0 1;1 0 0;0 1 0; 0 1 1; 1 0 1; 0.5 0.5 1; 1 0.5 0.5];
                
                subplot(length(ROI_List),3,(nroi-1)*length(ROI_List)+1)
                if MP
                    hold on
                    for mp=1:length(treat_ID)
                        h_all = errorbar(t,MeanForceTot(:,mp),StdErrForceTot(:,mp),'Color',colorMP(mp,:),'LineWidth',2);
                        set(gca,'FontSize',14)
                    end
                else
                    %             plot(t,MeanForceTot,'-..','MarkerSize',15);
                    for it=1:length(TreatMouseTot)
                        hold on
                        if TreatMouseTot(it)   == treat_ID(1)
                            h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(1,:));
                        elseif TreatMouseTot(it)== treat_ID(2)
                            h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(2,:));
                        elseif TreatMouseTot(it)== treat_ID(3)
                            h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(3,:));
                        elseif TreatMouseTot(it)== treat_ID(4)
                            h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(4,:));
                        elseif TreatMouseTot(it)== treat_ID(5)
                            h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(5,:));
                        elseif TreatMouseTot(it)== treat_ID(6)
                            h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(6,:));
                        elseif TreatMouseTot(it)== treat_ID(7)
                            h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(7,:));
                        end
                        set(gca,'FontSize',14)
                    end
                end
                title([ParNameToPlot{ParForce},' ',ROINameToPlot{ROI}])
                ylabel('Force Par')
                xlabel('Day')
                xlim([t(1)-1 t(end)+1]);
                set(gca,'xTick',t)
                set(gca,'xTickLabel',(num2str(t)))
                
                
                subplot(length(ROI_List),3,(nroi-1)*length(ROI_List)+2)
                if MP
                    hold on
                    for mp=1:length(treat_ID)
                        h_all = errorbar(t,MeanFluoTot(:,mp),StdErrFluoTot(:,mp),'Color',colorMP(mp,:),'LineWidth',2)
                        set(gca,'FontSize',14)
                    end
                else
                    %             plot(t,MeanFluoTot,'-..','MarkerSize',15);
                    for it=1:length(TreatMouseTot)
                        hold on
                        if TreatMouseTot(it)==treat_ID(1)
                            h_all = plot(t,MeanFluoTot(:,it),'-..','MarkerSize',15,'Color',colorMP(1,:));
                        elseif TreatMouseTot(it)==treat_ID(2)
                            h_all = plot(t,MeanFluoTot(:,it),'-..','MarkerSize',15,'Color',colorMP(2,:));
                        elseif TreatMouseTot(it)==treat_ID(3)
                            h_all = plot(t,MeanFluoTot(:,it),'-..','MarkerSize',15,'Color',colorMP(3,:));
                        elseif TreatMouseTot(it)==treat_ID(4)
                            h_all = plot(t,MeanFluoTot(:,it),'-..','MarkerSize',15,'Color',colorMP(4,:));
                        elseif TreatMouseTot(it)== treat_ID(5)
                            h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(5,:));
                        elseif TreatMouseTot(it)== treat_ID(6)
                            h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(6,:));
                        elseif TreatMouseTot(it)== treat_ID(7)
                            h_all = plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(7,:));
                        end
                        set(gca,'FontSize',14)
                    end
                end
                
                title([ParNameToPlot{ParFluo},' ',ROINameToPlot{ROI}])
                ylabel('Fluo Par')
                xlabel('Day')
                xlim([t(1)-1 t(end)+1]);
                set(gca,'xTick',t)
                set(gca,'xTickLabel',(num2str(t)))
                
                
                subplot(length(ROI_List),3,(nroi-1)*length(ROI_List)+3)
                if MP
                    hold on
                    for mp=1:length(treat_ID)
                        erProp = ((StdErrForceTot(:,mp)./MeanForceTot(:,mp)) + ((StdErrFluoTot(:,mp)./MeanFluoTot(:,mp)))) .* (MeanForceTot(:,mp)./MeanFluoTot(:,mp));
                        h_all = errorbar(t,MeanForceTot(:,mp)./MeanFluoTot(:,mp),erProp ,'Color',colorMP(mp,:),'LineWidth',2);
                        set(gca,'FontSize',14)
                    end
                else
                    h_all = plot(t,MeanForceTot./MeanFluoTot,'-..','MarkerSize',15);
                end
                title([ParNameToPlot{ParForce},'/',ParNameToPlot{ParFluo},'',ROINameToPlot{ROI}])
                ylabel('Force/Fluo')
                xlabel('Day')
                xlim([t(1)-1 t(end)+1]);
                set(gca,'xTick',t)
                set(gca,'xTickLabel',(num2str(t)))
            end
%             legend(Filename)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Plot settimanale  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if MP            
            %Plot settimanale (tutte e 4 settimane) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(MeanForce_Buf,1)>5
                
                weekSpace = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20];
                weekNum   = [1:size(weekSpace,1)];
                
                %Plot settimanale (ultima settimana) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                weekSpace = [1 2 3 4 5];
                weekNum   = [1:size(weekSpace,1)];
                
            end
            
            figure(Fig_weeks);
            
            hold on
            colorMP = [0 0 1;1 0 0;0 1 0; 0 1 1; 1 0 1; 0.5 0.5 1; 1 0.5 0.5];
                        
            %%%%initialization
            %presTreat = unique(TreatMouseTot);
            presTreat = treat_ID;
            lenTreat  = length(presTreat);
            
            medTot_ForceTot_Treat_Wk    = zeros(lenTreat,length(weekNum));
            stdTot_ForceTot_Treat_Wk    = zeros(lenTreat,length(weekNum));
            errstdTot_ForceTot_Treat_Wk = zeros(lenTreat,length(weekNum));
            
            medTot_FluoTot_Treat_Wk     = zeros(lenTreat,length(weekNum));
            stdTot_FluoTot_Treat_Wk     = zeros(lenTreat,length(weekNum));
            errstdTot_FluoTot_Treat_Wk  = zeros(lenTreat,length(weekNum));
            %%%%%
            
            
            %treatment
            for tr=1:lenTreat
                
                %to save excel Stat
                ExcelSaveForce = [];
                ExcelSaveFluo  = [];
                
                %prendo solo le colonne che si riferiscono al treatment
                MeanForceTot_Treat = MeanForce_Buf(:,TreatMouseTot==presTreat(tr));
                %prendo solo le colonne che si riferiscono al treatment
                MeanFluoTot_Treat  = MeanFluo_Buf(:,TreatMouseTot==presTreat(tr));
                
                for iw=1:length(weekNum)
                    
                    % % % % %
                    %media sui giorni (un valore per animale)
                    % % % % %
                    
                    %%%prendo solo le righe che si riferiscono alla week
                    MeanForceTot_Treat_Wk =  MeanForceTot_Treat(weekSpace(iw,:),:);
                    %%%prendo solo le righe che si riferiscono alla week
                    MeanFluoTot_Treat_Wk  =  MeanFluoTot_Treat(weekSpace(iw,:),:);
                    
                    %%%%median value for each animal in this week
                    if POOL == 1
                        med_ForceTot_Treat_Wk       =  reshape(MeanForceTot_Treat_Wk,size(MeanForceTot_Treat_Wk,1)*size(MeanForceTot_Treat_Wk,2),1);
                    else
%                         med_ForceTot_Treat_Wk       =  nanmedian(MeanForceTot_Treat_Wk,1)';
                        med_ForceTot_Treat_Wk       =  nanmean(MeanForceTot_Treat_Wk,1)';
                    end
                    
                    %%%total median value for the animals in this week
%                     medTot_ForceTot_Treat_Wk(tr,iw)       =  nanmedian(med_ForceTot_Treat_Wk);
                    medTot_ForceTot_Treat_Wk(tr,iw)       =  nanmean(med_ForceTot_Treat_Wk);
                    stdTot_ForceTot_Treat_Wk(tr,iw)       =  nanstd(med_ForceTot_Treat_Wk,[]);
                    errstdTot_ForceTot_Treat_Wk(tr,iw)    =  stdTot_ForceTot_Treat_Wk(tr,iw) /sqrt(length(med_ForceTot_Treat_Wk));
                    
                    
                    %%%median value for each animal in this week
                    if POOL == 1
                        med_FluoTot_Treat_Wk    =  reshape( MeanFluoTot_Treat_Wk,size( MeanFluoTot_Treat_Wk,1)*size(MeanFluoTot_Treat_Wk,2),1);
                    else
%                         med_FluoTot_Treat_Wk    =  nanmedian(MeanFluoTot_Treat_Wk,1)';
                        med_FluoTot_Treat_Wk    =  nanmean(MeanFluoTot_Treat_Wk,1)';
                    end
                    
                    %%%total median value for the animals in this week
%                     medTot_FluoTot_Treat_Wk(tr,iw)        =  nanmedian(med_FluoTot_Treat_Wk);
                    medTot_FluoTot_Treat_Wk(tr,iw)        =  nanmean(med_FluoTot_Treat_Wk);
                    stdTot_FluoTot_Treat_Wk(tr,iw)        =  nanstd(med_FluoTot_Treat_Wk,[]);
                    errstdTot_FluoTot_Treat_Wk(tr,iw)     =  stdTot_FluoTot_Treat_Wk(tr,iw) /sqrt(length(med_ForceTot_Treat_Wk));
                    
                    
                    %to save excel Stat
                    ExcelSaveForce = [ExcelSaveForce, med_ForceTot_Treat_Wk];
                    ExcelSaveFluo  = [ExcelSaveFluo,  med_FluoTot_Treat_Wk];
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%   Plot    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %force
                hold on
                subplot(length(ROI_List),2,(nroi-1)*(length(ROI_List)-1)+1)
                h = barwitherr(errstdTot_ForceTot_Treat_Wk , medTot_ForceTot_Treat_Wk);
                set(gca,'FontSize',14)
                title([ParNameToPlot{ParForce},' ',ROINameToPlot{ROI}])
                ylabel('Force Par')
                if weekNum == 1
                    set(gca,'xTick',[1 2 3 4 5])
                    %set(gca,'xTickLabel',{'ctrl','strk','reh4w','reh1w','rbt4w','rbt1w','Sti4w'})
                    set(gca,'xTickLabel',{'ctrl','strk','reh4w','reh1w','tox'})
                    xlabel('Week')                    
                else
                    set(gca,'xTick',weekNum)
                    set(gca,'xTickLabel',(num2cell(weekNum )))
                    xlabel('Week')
                end
                    
                %fluo
                hold on
                subplot(length(ROI_List),2,(nroi-1)*(length(ROI_List)-1)+2)
                h_week = barwitherr(errstdTot_FluoTot_Treat_Wk , medTot_FluoTot_Treat_Wk);
                set(gca,'FontSize',14)
                title([ParNameToPlot{ParFluo},' ',ROINameToPlot{ROI}])
                ylabel('Fluo Par')
                if weekNum == 1
                    set(gca,'xTick',[1 2 3 4 5])
                    %set(gca,'xTickLabel',{'ctrl','strk','reh4w','reh1w','rbt4w','rbt1w','Sti4w'})
                    set(gca,'xTickLabel',{'ctrl','strk','reh4w','reh1w','tox'})
                    xlabel('Week')                    
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
                
                if XLS_SAVE
                xlswrite([Folder_Where_Save,'\','DataForceStat_',ExcelStr,'_',ParNameToPlot{ParForce},'_',ROINameToPlot{nroi}],ExcelSaveForce ,tr)
                xlswrite([Folder_Where_Save,'\','DataFluoStat_',ExcelStr,'_',ParNameToPlot{ParForce},'_',ROINameToPlot{nroi}],ExcelSaveFluo ,tr)
%                 xlswrite(['C:\Users\Stefano\Desktop\Folder_Image_And_Stat_OGGI_v2\','DataForceStat_',ExcelStr,'_',ParNameToPlot{ParForce},'_',ROINameToPlot{nroi}],ExcelSaveForce ,tr)
%                 xlswrite(['C:\Users\Stefano\Desktop\Folder_Image_And_Stat_OGGI_v2\','DataFluoStat_',ExcelStr,'_',ParNameToPlot{ParForce},'_',ROINameToPlot{nroi}],ExcelSaveFluo ,tr)
                end
                %
                
            end
             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
    end %end NumAnimals
    
end %end ROI_List

% saveas(h_all, [Folder_Where_Save,'\','Fig_',ParNameToPlot{ParForce}],'fig');
% saveas(h_all, [Folder_Where_Save,'\','Fig_',ParNameToPlot{ParForce}],'jpg');
% saveas(h_week,[Folder_Where_Save,'\','Fig_Weekly_',ParNameToPlot{ParForce}],'fig');
% saveas(h_week,[Folder_Where_Save,'\','Fig_Weekly_',ParNameToPlot{ParForce}],'jpg');

clearvars -except indexipar PAR_INT XLS_SAVE
clc

end


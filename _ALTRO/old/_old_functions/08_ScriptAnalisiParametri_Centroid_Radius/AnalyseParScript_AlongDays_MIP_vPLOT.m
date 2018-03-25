%load dataMouseGCamp_Median -> MeanPar

if ~exist('MeanPar')
    error('LOAD the matrix "dataMouseGCamp_Median_Sel....mat" ')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SelROI=1;

if SelROI == 0
    %%%select what you want to plot -> MANUAL ROIs
    ParNameToPlot = {'Trial','Status','Diff Start Time','Duration','Max','Diff Max Time','FWHM','AUC','PtPAmp','nSubPeaks'};
    ROINameToPlot = {'Motor','Barrel','Occipital','Motor Small','Barrel Small','Occipital Small','Whole'};
    %% ROI of interest
    ROI = 6;
    %%%
elseif SelROI == 1
    %%%select what you want to plot -> MIP Based ROIs
    ParNameToPlot = {'Trial','Status','Diff Start Time','Duration','Max','Diff Max Time','FWHM','AUC','PtPAmp','nSubPeaks'};
    ROINameToPlot = {'Max MIP','Sum MIP'};
    %% ROI of interest
    ROI_List = [1 2];
    %%%
elseif SelROI == 2
    %%%select what you want to plot -> Centroid Based ROIs
    ParNameToPlot = {'Trial','Status','Diff Start Time','Duration','Max','Diff Max Time','FWHM','AUC','PtPAmp','nSubPeaks'};
    ROINameToPlot = {'Motor','Barrel','Occipital','Motor Small','Barrel Small','Occipital Small','Whole','Main ROI'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compare pair of different parameters
CompDiffPars = 1;
% CompDiffPars = 0;
if CompDiffPars == 0
    ParForce = 3;
    ParFluo = ParForce;
else
    ParForce = 3;
    ParFluo  = 8;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%select if plotting all of the animals or their mean (stroke vs control vs rehab)
% MP = 1; % ->stroke vs control vs rehab
MP = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumAnimals = size(MeanPar{3,1},1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Days to take
ListDaysToTake = [1 2 3 4 5];
% ListDaysToTake = [2 5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumPar  = size(MeanPar,1);
NumDay  = size(MeanPar,2);




fig_Par = figure('Name','Fig Parameters over the days ot traing');

for nroi=1:length(ROI_List)
    
    ROI = ROI_List(nroi);
    
    %Load
    MeanForceTot = [];
    StdForceTot  = [];
    MeanFluoTot  = [];
    StdFluoTot   = [];
    
    for na=1:NumAnimals
        
        MeanForce = [];
        StdForce  = [];
        MeanFluo  = [];
        StdFluo   = [];
        
        for nd=1:NumDay
            
            if sum(ismember(ListDaysToTake,nd))>0
                
                if ~isempty(MeanPar{ParForce,nd}{na,1})
                    MeanForce = [MeanForce; MeanPar{ParForce,nd}{na,1}(ROI,1)];
                    StdForce  = [StdForce;  MeanPar{ParForce,nd}{na,1}(ROI,2)];
                    
                    MeanFluo  = [MeanFluo;  MeanPar{ParFluo,nd}{na,1}(ROI,3)];
                    StdFluo   = [StdFluo;   MeanPar{ParFluo,nd}{na,1}(ROI,4)];
                    
                    NameMouse  = MeanPar{ParForce,nd}{na,2};
                    TreatMouse = MeanPar{ParForce,nd}{na,4};
                else
                    MeanForce = [MeanForce; NaN];
                    StdForce  = [StdForce;  NaN];
                    
                    MeanFluo  = [MeanFluo;  NaN];
                    StdFluo   = [StdFluo;   NaN];
                    
                end
            end
            
        end
        
        %tutti i valori medi e std del parametro nei vari giorni per quell'animale na
        %(diventano colonne delle matrici _Tot)
        MeanForceTot = [MeanForceTot MeanForce];
        %     StdForceTot  = [StdForceTot  StdForce];
        
        MeanFluoTot = [MeanFluoTot MeanFluo];
        %     StdFluoTot  = [StdFluoTot  StdFluo];
        
        NameMouseTot{na,1} = NameMouse;
        TreatMouseTot(na,1)  = TreatMouse;
        
        
        
        if na==NumAnimals
            
            if MP
                
                MeanForce_Buf = MeanForceTot;
                MeanFluo_Buf  = MeanFluoTot;
                
                clear MeanForceTot MeanFluoTot StdForceTot StdFluoTot
                
                for tr=0:2
                    
                    StdForceTot(:,tr+1)     = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2);
                    StdErrForceTot(:,tr+1)  = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2)/sqrt(sum(TreatMouseTot==tr));
                    %                 MeanForceTot(:,tr+1) = nanmean(MeanForce_Buf(:,TreatMouseTot==tr),2);
                    MeanForceTot(:,tr+1) = nanmedian(MeanForce_Buf(:,TreatMouseTot==tr),2);
                    
                    StdFluoTot(:,tr+1)      = nanstd(MeanFluo_Buf(:,TreatMouseTot==tr),[],2);
                    StdErrFluoTot(:,tr+1)  = nanstd(MeanFluo_Buf(:,TreatMouseTot==tr),[],2)/sqrt(sum(TreatMouseTot==tr));
                    %                 MeanFluoTot(:,tr+1)  = nanmean(MeanFluo_Buf(:,TreatMouseTot==tr),2);
                    MeanFluoTot(:,tr+1)  = nanmedian(MeanFluo_Buf(:,TreatMouseTot==tr),2);
                    
                    
                end
                
                Filename = {'control','stroke','rehab'};
                
            else
                Filename = NameMouseTot;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Plot
            
            figure(fig_Par)
            
            if length(ListDaysToTake) ~= length(NumDay)
                t = ListDaysToTake';
            else
                t = [1:NumDay]';
            end
            colorMP = [0 0 1;1 0 0;0 1 0]*ROI/2;
            
            subplot(1,3,1)
            if MP
                hold on
                for mp=1:3
                    errorbar(t,MeanForceTot(:,mp),StdErrForceTot(:,mp),'Color',colorMP(mp,:),'LineWidth',2)
                    set(gca,'FontSize',14)
                end
            else
                %             plot(t,MeanForceTot,'-..','MarkerSize',15);
                for it=1:length(TreatMouseTot)
                    hold on
                    if TreatMouseTot(it)==0
                        plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(1,:));
                    elseif TreatMouseTot(it)==1
                        plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(2,:));
                    elseif TreatMouseTot(it)==2
                        plot(t,MeanForceTot(:,it),'-..','MarkerSize',15,'Color',colorMP(3,:));
                    end
                    set(gca,'FontSize',14)
                end
            end
            title([ParNameToPlot{ParForce},' ',ROINameToPlot{ROI}])
            ylabel('Force')
            xlabel('Day')
            xlim([t(1)-1 t(end)+1]);
            set(gca,'xTick',t)
            set(gca,'xTickLabel',(num2str(t)))
            
            
            subplot(1,3,2)
            if MP
                hold on
                for mp=1:3
                    errorbar(t,MeanFluoTot(:,mp),StdErrFluoTot(:,mp),'Color',colorMP(mp,:),'LineWidth',2)
                    set(gca,'FontSize',14)
                end
            else
                %             plot(t,MeanFluoTot,'-..','MarkerSize',15);
                for it=1:length(TreatMouseTot)
                    hold on
                    if TreatMouseTot(it)==0
                        plot(t,MeanFluoTot(:,it),'-..','MarkerSize',15,'Color',colorMP(1,:));
                    elseif TreatMouseTot(it)==1
                        plot(t,MeanFluoTot(:,it),'-..','MarkerSize',15,'Color',colorMP(2,:));
                    elseif TreatMouseTot(it)==2
                        plot(t,MeanFluoTot(:,it),'-..','MarkerSize',15,'Color',colorMP(3,:));
                    end
                    set(gca,'FontSize',14)
                end
            end
            
            title([ParNameToPlot{ParFluo},' ',ROINameToPlot{ROI}])
            ylabel('Fluo')
            xlabel('Day')
            xlim([t(1)-1 t(end)+1]);
            set(gca,'xTick',t)
            set(gca,'xTickLabel',(num2str(t)))
            
            
            subplot(1,3,3)
            if MP
                hold on
                for mp=1:3
                    erProp = ((StdErrForceTot(:,mp)./MeanForceTot(:,mp)) + ((StdErrFluoTot(:,mp)./MeanFluoTot(:,mp)))) .* (MeanForceTot(:,mp)./MeanFluoTot(:,mp));
                    errorbar(t,MeanForceTot(:,mp)./MeanFluoTot(:,mp),erProp ,'Color',colorMP(mp,:),'LineWidth',2)
                    set(gca,'FontSize',14)
                end
            else
                plot(t,MeanForceTot./MeanFluoTot,'-..','MarkerSize',15);
            end
            title([ParNameToPlot{ParForce},'/',ParNameToPlot{ParFluo},'',ROINameToPlot{ROI}])
            ylabel('Force/Fluo')
            xlabel('Day')
            xlim([t(1)-1 t(end)+1]);
            set(gca,'xTick',t)
            set(gca,'xTickLabel',(num2str(t)))
            
            legend(Filename)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
    end %end NumAnimals
end %end ROI_List






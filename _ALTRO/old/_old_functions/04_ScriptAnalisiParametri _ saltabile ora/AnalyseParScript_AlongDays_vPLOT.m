%load dataMouseGCamp_Mean -> MeanPar

%%%select what you want to plot
ParNameToPlot = {'Trial','Status','Diff Start Time','Duration','Max','Diff Max Time','FWHM','AUC','PtPAmp','nSubPeaks'};
ROINameToPlot = {'Motor','Barrel','Occipital','Motor Small','Barrel Small','Occipital Small','Whole'};
%%% compare pair od differene parameters
CompDiffPars = 1;
% CompDiffPars = 0;
%%% compare pair of different parameters
if CompDiffPars == 0
    ParForce = 3;
    ParFluo = ParForce;
else
    ParForce = 8;
    ParFluo  = 8;    
end
%% ROI of interest
ROI = 4;
%%%

%%%select if plotting all of the animals or their mean (stroke vs control)
MP = 1;
% MP = 0;
%%%
%%%
NumAnimals = 6;
%%%

NumPar  = size(MeanPar,1);
NumDay  = size(MeanPar,2);


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
        
        if ~isempty(MeanPar{ParForce,nd}{na,1})
            MeanForce = [MeanForce; MeanPar{ParForce,nd}{na,1}(ROI,1)];
            StdForce  = [StdForce;  MeanPar{ParForce,nd}{na,1}(ROI,2)];
            
            MeanFluo  = [MeanFluo;  MeanPar{ParFluo,nd}{na,1}(ROI,3)];
            StdFluo   = [StdFluo;   MeanPar{ParFluo,nd}{na,1}(ROI,4)];
            
            %         MeanFluo  = [MeanFluo;  MeanPar{Par,nd}{na,1}(ROI,3)./MeanPar{Par,nd}{na,1}(ROI,1)];
            %         StdFluo   = [StdFluo;   MeanPar{Par,nd}{na,1}(ROI,4)];
            
            
            NameMouse  = MeanPar{ParForce,nd}{na,2};
            TreatMouse = MeanPar{ParForce,nd}{na,4};
        else
            MeanForce = [MeanForce; NaN];
            StdForce  = [StdForce;  NaN];
            
            MeanFluo  = [MeanFluo;  NaN];
            StdFluo   = [StdFluo;   NaN];
            
        end
        
    end
    
    MeanForceTot = [MeanForceTot MeanForce];
    StdForceTot  = [StdForceTot  StdForce];
    
    MeanFluoTot = [MeanFluoTot MeanFluo];
    StdFluoTot  = [StdFluoTot  StdFluo];
    
    NameMouseTot{na,1} = NameMouse;
    TreatMouseTot(na,1)  = TreatMouse;
    
    

    if na==NumAnimals
        
        if MP
            
            MeanForce_Buf = MeanForceTot;
            MeanFluo_Buf  = MeanFluoTot;
            
            clear MeanForceTot MeanFluoTot StdForceTot StdFluoTot
            
            for tr=0:1
                
                StdForceTot(:,tr+1)     = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2);
                StdErrForceTot(:,tr+1)  = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2)/sqrt(sum(TreatMouseTot==tr));
%                 MeanForceTot(:,tr+1) = nanmean(MeanForce_Buf(:,TreatMouseTot==tr),2);
                MeanForceTot(:,tr+1) = nanmedian(MeanForce_Buf(:,TreatMouseTot==tr),2);
            
                StdFluoTot(:,tr+1)      = nanstd(MeanFluo_Buf(:,TreatMouseTot==tr),[],2);
                StdErrFluoTot(:,tr+1)  = nanstd(MeanFluo_Buf(:,TreatMouseTot==tr),[],2)/sqrt(sum(TreatMouseTot==tr));
%                 MeanFluoTot(:,tr+1)  = nanmean(MeanFluo_Buf(:,TreatMouseTot==tr),2);
                MeanFluoTot(:,tr+1)  = nanmedian(MeanFluo_Buf(:,TreatMouseTot==tr),2);
                
                
            end
            
            Filename = {'control','stroke'};
            
        else
            Filename = NameMouseTot;
        end
        
        %Plot
        figure
        
        t = [1:NumDay]';
        subplot(1,3,1)
        if MP
            hold on
            colorMP = [0 0 1;1 0 0];
            for mp=1:2
                errorbar(t,MeanForceTot(:,mp),StdErrForceTot(:,mp),'Color',colorMP(mp,:))
            end
        else
            plot(t,MeanForceTot,'-..','MarkerSize',15);
        end
        title([ParNameToPlot{ParForce},' ',ROINameToPlot{ROI}])
        ylabel('Force')
        xlabel('Day')
        xlim([t(1)-1 t(end)+1]);
        
        subplot(1,3,2)
        if MP
            hold on
            colorMP = [0 0 1;1 0 0];
            for mp=1:2
                errorbar(t,MeanFluoTot(:,mp),StdErrFluoTot(:,mp),'Color',colorMP(mp,:))
            end
        else
            plot(t,MeanFluoTot,'-..','MarkerSize',15);
        end
        
        title([ParNameToPlot{ParFluo},' ',ROINameToPlot{ROI}])
        ylabel('Fluo')
        xlabel('Day')
        xlim([t(1)-1 t(end)+1]);      
        
        
        subplot(1,3,3)
        if MP
            hold on
            colorMP = [0 0 1;1 0 0];
            for mp=1:2
                erProp = ((StdErrForceTot(:,mp)./MeanForceTot(:,mp)) + ((StdErrFluoTot(:,mp)./MeanFluoTot(:,mp)))) .* (MeanForceTot(:,mp)./MeanFluoTot(:,mp));
                errorbar(t,MeanForceTot(:,mp)./MeanFluoTot(:,mp),erProp ,'Color',colorMP(mp,:))
            end
        else
            plot(t,MeanForceTot./MeanFluoTot,'-..','MarkerSize',15);
        end
        title([ParNameToPlot{ParForce},'/',ParNameToPlot{ParFluo},'',ROINameToPlot{ROI}])
        ylabel('Force/Fluo')
        xlabel('Day')
        xlim([t(1)-1 t(end)+1]);     
       
        legend(Filename)
        
    end

end






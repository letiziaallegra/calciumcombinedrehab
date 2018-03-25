%load dataMouseGCamp_Mean -> MeanPar

%%%select what you want to plot
ParNameToPlot = {'Trial','Target Time','Force','Sub Mov','Attempts'};
ROINameToPlot = {'Robot Par'};

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

%for paramenters
figure
for pi=1:4
     
    ParForce = pi+1; 
    ROI = 1;
    
    MeanForceTot = [];
    StdForceTot = [];
     
    for na=1:NumAnimals
        
        MeanForce = [];
        StdForce  = [];        
        
        for nd=1:NumDay
            
            if ~isempty(MeanPar{ParForce,nd}{na,1})
                
                MeanForce = [MeanForce; MeanPar{ParForce,nd}{na,1}(ROI,1)];
                StdForce  = [StdForce;  MeanPar{ParForce,nd}{na,1}(ROI,2)];
                
                
                %         MeanFluo  = [MeanFluo;  MeanPar{Par,nd}{na,1}(ROI,3)./MeanPar{Par,nd}{na,1}(ROI,1)];
                %         StdFluo   = [StdFluo;   MeanPar{Par,nd}{na,1}(ROI,4)];
                
                
                NameMouse  = MeanPar{ParForce,nd}{na,2};
                TreatMouse = MeanPar{ParForce,nd}{na,4};
            else
                MeanForce = [MeanForce; NaN];
                StdForce  = [StdForce;  NaN];
                
            end
            
        end
        
        MeanForceTot = [MeanForceTot MeanForce];
        StdForceTot  = [StdForceTot  StdForce];
        
        NameMouseTot{na,1} = NameMouse;
        TreatMouseTot(na,1)  = TreatMouse;
        
        
        
        if na==NumAnimals
            
            if MP
                
                MeanForce_Buf = MeanForceTot;
                
                clear MeanForceTot MeanFluoTot StdForceTot StdFluoTot
                
                for tr=0:1
                    
                    StdForceTot(:,tr+1)     = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2);
                    StdErrForceTot(:,tr+1)  = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2)/sqrt(sum(TreatMouseTot==tr));
                    %                 MeanForceTot(:,tr+1) = nanmean(MeanForce_Buf(:,TreatMouseTot==tr),2);
                    MeanForceTot(:,tr+1) = nanmedian(MeanForce_Buf(:,TreatMouseTot==tr),2);
                    
                    
                end
                
                Filename = {'control','stroke'};
                
            else
                Filename = NameMouseTot;
            end
            
            
            %Plot
            subplot(2,2,pi)
            t = [1:NumDay]';
            if MP
                hold on
                colorMP = [0 0 1;1 0 0];
                for mp=1:2
                    errorbar(t,MeanForceTot(:,mp),StdErrForceTot(:,mp),'Color',colorMP(mp,:))
                end
            else
                plot(t,MeanForceTot,'-..','MarkerSize',15);
            end
%             title([ParNameToPlot{ParForce},' ',ROINameToPlot{ROI}])
            ylabel(ParNameToPlot{ParForce})
            xlabel('Day')
            xlim([t(1)-1 t(end)+1]);
            
            if   pi == 4
                legend(Filename)
            end
            
        end
        
    end
end






%load dataMouseGCamp_RobotPar_Median ... -> MeanPar

if ~exist('MeanPar')
    error('LOAD the matrix "dataMouseGCamp_RobotPar_Median....mat" ')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%select what you want to plot
ParNameToPlot = {'Trial','Target Time','Force','Sub Mov','Attempts'};
ROINameToPlot = {'Robot Par'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%select if plotting all of the animals or their mean (stroke vs control vs rehab)
MP = 1; % ->stroke vs control vs rehab
% MP = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumAnimals = size(MeanPar{2,1},1);
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumPar  = size(MeanPar,1);
NumDay  = size(MeanPar,2);
MeanForceTot = [];
StdForceTot  = [];

%for the 4 parameters from Robot
Fig_all_days = figure('Name','Time course motor performance on the platfom');
Fig_weeks    = figure('Name','Time course weekly motor performance on the platfom');
for pi=1:4
     
    ParForce = pi+1; 
    ROI = 1;
    
    MeanForceTot = [];
    StdForceTot = [];
     
    for na=1:NumAnimals
        
        MeanForce = [];
        StdForce  = [];
        MeanForceTot_Wk   = [];
        StdErrForceTot_Wk = [];
        
        for nd=1:NumDay
            
            if ~isempty(MeanPar{ParForce,nd})
                
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
                
                for tr=0:2
                    
                    %nella matrice MeanForce_Buf prendo solo le colonne che si riferiscono a un certo treatment
                    StdForceTot(:,tr+1)     = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2);
                    StdErrForceTot(:,tr+1)  = nanstd(MeanForce_Buf(:,TreatMouseTot==tr),[],2)/sqrt(sum(TreatMouseTot==tr));
                    %                 MeanForceTot(:,tr+1) = nanmean(MeanForce_Buf(:,TreatMouseTot==tr),2);
                    MeanForceTot(:,tr+1) = nanmedian(MeanForce_Buf(:,TreatMouseTot==tr),2);
                    
                    
                end
                
                Filename = {'control','stroke','rehab'};
                
            else
                Filename = NameMouseTot;
            end
            
            
            %Plot
            figure(Fig_all_days)
            subplot(2,2,pi)
            t = [cell2mat(MeanPar(1,1)):cell2mat(MeanPar(1,end))]';
            if MP
                hold on
                colorMP = [0 0 1;1 0 0;0 1 0];
                for mp=1:3
                    errorbar(t,MeanForceTot(:,mp),StdErrForceTot(:,mp),'Color',colorMP(mp,:),'LineWidth',2)
                    set(gca,'FontSize',14)
                end
            else
                plot(t,MeanForceTot,'-..','MarkerSize',15);
            end
%             title([ParNameToPlot{ParForce},' ',ROINameToPlot{ROI}])
            ylabel(ParNameToPlot{ParForce})
            xlabel('Day')
            xlim([t(1)-1 t(end)+1]);
            set(gca,'xTick',t)
            set(gca,'xTickLabel',(num2str(t)))
            
            if   pi == 4
                legend(Filename)
            end
            
            
            
            %Plot settimanale (più di una settimana) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(MeanForceTot,1)>5
                
                weekSpace = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20];
                weekNum   = [1:size(weekSpace,1)];
                
                figure(Fig_weeks)
                subplot(2,2,pi)
                              
                
                if MP
                    hold on
                    colorMP = [0 0 1;1 0 0;0 1 0];
                    for mp=1:3
                        
                        for iw = 1:weekNum(end)
                            
                            MeanForceTot_Wk(iw,mp)   = [nanmean(MeanForceTot(weekSpace(iw,:),mp))];
%                             MeanForceTot_Wk(iw,mp)   = [nanmedian(MeanForceTot(weekSpace(iw,:),mp))];
                            
                            StdErrForceTot_Wk(iw,mp) = [nanstd((MeanForceTot(weekSpace(iw,:),mp)),[])/sqrt(length(MeanForceTot(weekSpace(iw,:),mp)))];
                                                        
                        end
                        
                        errorbar(weekNum, MeanForceTot_Wk(:,mp),StdErrForceTot_Wk(:,mp),'Color',colorMP(mp,:),'LineWidth',2)
                        set(gca,'FontSize',14)
                    end
                else
%                     plot(t,MeanForceTot,'-..','MarkerSize',15);
                end
                %             title([ParNameToPlot{ParForce},' ',ROINameToPlot{ROI}])
                ylabel(ParNameToPlot{ParForce})
                xlabel('Week')
                xlim([weekNum(1)-1 weekNum(end)+1]);
                set(gca,'xTick',weekNum)
                set(gca,'xTickLabel',(num2cell(weekNum )))
                                
                if   pi == 4
                    legend(Filename)
                end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
            end
        end
        
    end
end






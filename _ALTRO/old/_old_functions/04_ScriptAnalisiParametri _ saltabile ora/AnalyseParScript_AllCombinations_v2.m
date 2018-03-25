%Script To analyse Parameters from Peaks
%AnalyseParScript
%Analyis performed on the parameters of dataGCamp.PeaksPar_Fx_Fluo
%        each column of dataGCamp.PeaksPar{:,1} refers to a peak:
%                   column 1 -> trial
%                   column 2 -> status
%                   column 3 -> peak start point 
%                   column 4 -> peak duration (fixed window)
%                   column 5 -> max value peak
%                   column 6 -> max value peak point (IstAmpMaxWSig)
%                   column 7 -> full width at half maximum (FWHM)
%                   column 8 -> area under the peak curve (AUC)
%                   column 9 -> peak-to-peak amplitude (PtPAmp)
%                   column 10 -> number of sub-peaks (nSubPeaks)

% close all
clc


%consider the previous list
ParNameToTake = {'Trial','Status','Start Time','Duration','Max','Max Time','FWHM','AUC','PtPAmp','nSubPeaks'};
% ParToTake     = [3 4 5 6 7 8 9 10];
ParToTake     = [3 5 8 10];
ColorsStart   = [0 0 1; 0 1 0; 1 0 0; 1 0 1; 0 1 1; 1 1 0; 0 0 0];
NumPar        = length(ParToTake);
% %ROI
% ROINameToTake = {'Motor','Barrel','Occipital','Motor Small','Barrel Small','Occipital Small','Whole'};
% ParToTake     = [1 2 3 4 5 6 7];


if isfield(dataGCamp,'PeaksPar_Fx_Fluo')
    
    %num signals to analyse: force + fluo ROIs
    NumSig        = length(dataGCamp.PeaksPar_Fx_Fluo);
    
    %for each of the parameters (one figure for each parameter) -> Par of the force
    for ip=1:NumSig
        
        %info Par
        Curr_Par      = ParToTake(ip);
        Curr_Par_Name = ParNameToTake{Curr_Par};
        
        figure('Name','Parameters from Peaks')
        
        %for each Par-couple of the FORCE and i-th FLUO ROI
        for id=2:NumSig
            
            if ~isempty(dataGCamp.PeaksPar_Fx_Fluo{1,id-1})
                
                %force par (force is fixed) -> the peaks selected can change
                ForcePar = dataGCamp.PeaksPar_Fx_Fluo{1,id-1}(:,Curr_Par);
                
                for ipROI=1:NumPar
                    
                    %fluo par
                    CurrPar_ROI      = ParToTake(ipROI);
                    CurrPar_ROI_Name = ParNameToTake{CurrPar_ROI};
                    FluoPar          = dataGCamp.PeaksPar_Fx_Fluo{2,id-1}(:,CurrPar_ROI);
                    
                    hold on
                    subplot(3,3,ipROI)
                    ParToPlot = sortrows( [ForcePar, FluoPar],1);
                    plot(ParToPlot(:,1), ParToPlot(:,2), 'Color',ColorsStart(id-1,:),'Marker','o','MarkerSize',6,'MarkerFaceColor',ColorsStart(id-1,:))
                    %                 scatter(ForcePar, FluoPar,10, ColorsStart(id-1,:),'filled')
                    ylabel([CurrPar_ROI_Name,' (Fluo)'])
                    xlabel([Curr_Par_Name,' (Force)'])
                end
                legend(dataGCamp.Info.ROI)
            end
            
        end
        
    end
    
else
    display('Computation of parameters has to be done before the analysis')
end
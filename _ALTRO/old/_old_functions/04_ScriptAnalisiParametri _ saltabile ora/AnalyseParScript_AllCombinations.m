%Script To analyse Parameters from Peaks
%AnalyseParScript
%Analyis performed on the parameters of dataGCamp.PeaksPar{:,1}
%        each column of dataGCamp.PeaksPar{:,1} refers to a peak:
%                   column 1 -> trial
%                   column 2 -> status
%                   column 3 -> peak start point
%                   column 4 -> peak duration
%                   column 5 -> max value peak
%                   column 6 -> max value peak point (IstAmpMaxWSig)
%                   column 7 -> full width at half maximum (FWHM)
%                   column 8 -> area under the peak curve (AUC)
%                   column 9 -> peak-to-peak amplitude (PtPAmp)

% close all
clc


%consider the previous list
ParNameToTake = {'Trial','Status','Start Time','Duration','Max','Max Time','FWHM','AUC','PtPAmp'};
ParToTake     = [3 4 5 6 7 8 9];
ColorsStart   = [0 0 1; 0 1 0; 1 0 0; 1 0 1; 0 1 1; 1 1 0; 0 0 0];
NumPar        = length(ParToTake);


if isfield(dataGCamp,'PeaksPar') && isfield(dataGCamp,'CorrespondanceForceFluoPeaks')
    
    %num signals to analyse: force + fluo ROIs
    NumSig        = length(dataGCamp.PeaksPar);
    
        
    %for each of the parameters (one figure for each parameter) -> Par of the force
    for ip=1:NumPar
        
        %info Par
        Curr_Par      = ParToTake(ip);
        Curr_Par_Name = ParNameToTake{Curr_Par};
        
        figure('Name','Parameters from Peaks')
        
        %for each couple Par of the force - all fluoROI Pars
        for id=2:NumSig
            
            %force peaks to take account of
            ind_ForcePeaksToTake = dataGCamp.CorrespondanceForceFluoPeaks{1,id-1}(:,1);
            %fluo peaks to take account of
            ind_FLuoPeaksToTake  = dataGCamp.CorrespondanceForceFluoPeaks{1,id-1}(:,2);
            
            %force par (force is fixed) -> the peaks selected can change
            ForcePar = dataGCamp.PeaksPar{1,1}(ind_ForcePeaksToTake,Curr_Par);
            
            for ipROI=1:NumPar
                %fluo par
                CurrPar_ROI      = ParToTake(ipROI);
                CurrPar_ROI_Name = ParNameToTake{CurrPar_ROI};
                FluoPar  = dataGCamp.PeaksPar{1,id}(ind_FLuoPeaksToTake,CurrPar_ROI);
                
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
    
else
    display('Computation of parameters has to be done before the analysis')
end
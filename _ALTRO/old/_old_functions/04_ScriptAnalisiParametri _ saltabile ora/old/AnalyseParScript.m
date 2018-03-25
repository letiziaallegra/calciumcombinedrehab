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

close all
clc


%consider the previous list
ParNameToTake = {'Trial','Status','Start Time','Duration','Max','Max Time','FWHM','AUC','PtPAmp'};
ParToTake     = [3 4 5 6 7 8 9];
ColorsStart   = [0 0 1; 0 1 0; 1 0 0; 1 0 1; 0 1 1; 1 1 0; 0 0 0];
NumPar        = length(ParToTake);

%load
% load('dataMouseGCamp_GCaMP9_stroke_03')

if isfield(dataGCamp,'PeaksPar') && isfield(dataGCamp,'CorrespondanceForceFluoPeaks')
   
    %num signals to analyse: force + fluo ROIs
    NumSig        = length(dataGCamp.PeaksPar);
    
    figure('Name','Parameters from Peaks')
    
    %for each of the parameters (one figure for each parameter)
    for ip=1:NumPar 
        
        %info Par
        Curr_Par_Name = ParNameToTake{ParToTake(ip)};
        Curr_Par      = ParToTake(ip);        
        
        subplot(3,3,ip)
%         figure('Name',Curr_Par_Name)
        
        %for each couple force-fluoROI (the loop is for the ROI (force is fixed)
        for id=2:NumSig
            
            %force peaks to take account of
            ind_ForcePeaksToTake = dataGCamp.CorrespondanceForceFluoPeaks{1,id-1}(:,1);
            %fluo peaks to take account of
            ind_FLuoPeaksToTake  = dataGCamp.CorrespondanceForceFluoPeaks{1,id-1}(:,2);
            
            %force par (force is fixed) -> the peaks selected can change
            ForcePar = dataGCamp.PeaksPar{1,1}(ind_ForcePeaksToTake,Curr_Par);
            %force par
            FluoPar  = dataGCamp.PeaksPar{1,id}(ind_FLuoPeaksToTake,Curr_Par);
            
            hold on
            
%             plot(ForcePar, FluoPar, 'Color',ColorsStart(id-1,:))
            scatter(ForcePar, FluoPar,10, ColorsStart(id-1,:),'filled')

        end
        
        title(Curr_Par_Name)
        ylabel('Par from FLUO Peaks')
        xlabel('Par from FORCE Peaks')
    
        
    end
    legend(dataGCamp.Info.ROI)
    
else
    display('Computation of parameters has to be done before the analysis')
end
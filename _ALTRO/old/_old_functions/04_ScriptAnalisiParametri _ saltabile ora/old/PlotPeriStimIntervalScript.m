%Script To plot peristimulus interval (the 0 is the star of the force peak)
%AnalyseParScript
%Analyis performed on the parameters of dataGCamp.PeaksPar{:,1}
%        each column of dataGCamp.PeaksPar{:,1} refers to a peak:

%                   column 3 -> peak start point
%                   column 4 -> peak duration


close all
clc

%load
% load('dataMouseGCamp_GCaMP9_stroke_03')

if isfield(dataGCamp,'PeaksPar') && isfield(dataGCamp,'CorrespondanceForceFluoPeaks')
    
    %num signals to analyse: force + fluo ROIs
    NumSig        = length(dataGCamp.PeaksPar);
    Fs            = dataGCamp.Info.Fs;
    
    %
    StartTime_Curr_Par = 3;
    Duration_Curr_Par  = 4;
    
    f = figure('Name','Parameters from Peaks')
    
    %for each force - fluoROI
    for id=2:NumSig
        
        %force peaks to take account of
        ind_ForcePeaksToTake = dataGCamp.CorrespondanceForceFluoPeaks{1,id-1}(:,1);
        %fluo peaks to take account of
        ind_FLuoPeaksToTake  = dataGCamp.CorrespondanceForceFluoPeaks{1,id-1}(:,2);
        
        %force par (force is fixed) -> the peaks selected can change
        StartTime_ForcePar = dataGCamp.PeaksPar{1,1}(ind_ForcePeaksToTake,StartTime_Curr_Par);
        
        %fluo start
        StartTime_FluoPar = dataGCamp.PeaksPar{1,id}(ind_FLuoPeaksToTake,StartTime_Curr_Par);
        %fluo duration
        Duration_FluoPar  = dataGCamp.PeaksPar{1,id}(ind_FLuoPeaksToTake,Duration_Curr_Par);
        
        
        hold on
        subplot(3,3,id-1)
            %%Fig whole Sig
            ss  = figure;
        
            %for each peaks
        for iNPeaks=1:length(StartTime_ForcePar)
            
            %start F
            stF = (StartTime_ForcePar(iNPeaks))*Fs;
            %start Fluo
            stFluo = (StartTime_FluoPar(iNPeaks))*Fs;
            %end Fluo
            enFluo = (StartTime_FluoPar(iNPeaks)+Duration_FluoPar(iNPeaks)-1)*Fs;
            
            Diff_start = stF-stFluo;
            if Diff_start<0
                stFluo = stFluo-abs(Diff_start)+1;
            elseif Diff_start>0
                stFluo = stFluo+abs(Diff_start)-1;
            end
            
            %fluo sig
            FluoROI_Sig = dataGCamp.fluoROI(:,id-1);
                %%Fig whole Sig 
                figure(ss)
                if iNPeaks==1
                    plot(FluoROI_Sig)
                end
            %Fluo interval to plot
            IntervalFluoROI_Sig = FluoROI_Sig(stFluo:enFluo);
                hold on
                plot(stFluo:enFluo,FluoROI_Sig(stFluo:enFluo),'r')
            %time
            t = [0:length(IntervalFluoROI_Sig)-1]'/Fs;
            
            figure(f)
            hold on
            plot(t,IntervalFluoROI_Sig)
            
        end
        ylabel('Fluo')
        xlabel('Time')
        ylim([0 20])
        xlim([-0.5 3])
        title(dataGCamp.Info.ROI{id-1})
        
    end
    
    
else
    display('Computation of parameters has to be done before the analysis')
end
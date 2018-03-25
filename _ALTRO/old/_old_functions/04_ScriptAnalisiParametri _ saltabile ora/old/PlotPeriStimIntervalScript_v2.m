%Script To plot peristimulus interval (the 0 is the star of the force peak)
%AnalyseParScript
%Analyis performed on the parameters of dataGCamp.PeaksPar{:,1}
%        each column of dataGCamp.PeaksPar{:,1} refers to a peak:

%                   column 3 -> peak start point
%                   column 4 -> peak duration


% close all
clc

%%%
Alignement = 'startForce';
% Alignement = 'maxForce';
%%%

%%%
% PlotMed = 1;
PlotMed = 0;
Med_store = [];
%%%


%load
% load('dataMouseGCamp_GCaMP9_stroke_03')

if isfield(dataGCamp,'PeaksPar') && isfield(dataGCamp,'CorrespondanceForceFluoPeaks')
    
    %num signals to analyse: force + fluo ROIs
    NumSig        = length(dataGCamp.PeaksPar);
    Fs            = dataGCamp.Info.Fs;
    %
    ForceSignal   = -dataGCamp.fx;
    
    %
    StartTime_Curr_Par = 3;
    Duration_Curr_Par  = 4;
    MaxTime_Curr_Par  = 6;
    
    f = figure('Name','Parameters from Peaks');
    
    %for each force - fluoROI
    for id=2:NumSig
        
        %force peaks to take account of
        ind_ForcePeaksToTake = dataGCamp.CorrespondanceForceFluoPeaks{1,id-1}(:,1);
        %fluo peaks to take account of
        ind_FLuoPeaksToTake  = dataGCamp.CorrespondanceForceFluoPeaks{1,id-1}(:,2);
        
        %force par (force is fixed) -> the peaks selected can change
        StartTime_ForcePar = dataGCamp.PeaksPar{1,1}(ind_ForcePeaksToTake,StartTime_Curr_Par);
        %force duration
        Duration_ForcePar  = dataGCamp.PeaksPar{1,1}(ind_ForcePeaksToTake,Duration_Curr_Par);
        %time of the max
        MaxTime_ForcePar  = dataGCamp.PeaksPar{1,1}(ind_ForcePeaksToTake,MaxTime_Curr_Par);
        
        %fluo start
        StartTime_FluoPar = dataGCamp.PeaksPar{1,id}(ind_FLuoPeaksToTake,StartTime_Curr_Par);
        %fluo duration
        Duration_FluoPar  = dataGCamp.PeaksPar{1,id}(ind_FLuoPeaksToTake,Duration_Curr_Par);
        %fluo sig
        FluoROI_Sig = dataGCamp.fluoROI(:,id-1);       
        %fluo sig
        FluoROI_Sig = cheb2LPfilt(FluoROI_Sig,6,2,Fs);

        
        hold on
        subplot(3,3,id-1)
            %%Fig whole Sig
            ss  = figure;
            subplot(211)
            plot(ForceSignal)
            subplot(212)
            plot(FluoROI_Sig)
        
            %for each peaks
        for iNPeaks=1:length(StartTime_ForcePar)
            
            %start Force
            stF = (StartTime_ForcePar(iNPeaks))*Fs;
            %end Force
            enF = (StartTime_ForcePar(iNPeaks)+Duration_ForcePar(iNPeaks))*Fs-1;
            %max Force
            maxF = (MaxTime_ForcePar (iNPeaks))*Fs;
                        
            %start Fluo
            stFluo = (StartTime_FluoPar(iNPeaks))*Fs;
            %end Fluo
            enFluo = (StartTime_FluoPar(iNPeaks)+Duration_FluoPar(iNPeaks))*Fs-1;
            
            %to align in agreement with 0 (t start of ForcePeak)
            if strcmp(Alignement,'startForce')
                Diff_start = stF-stFluo;
            elseif strcmp(Alignement,'maxForce')
                stF = maxF;
                Diff_start = stF-stFluo;                
            end
            %Force to plot
            IntervalForce_Sig = ForceSignal(stF:enF);
            %Force time
            tF = [0:length(IntervalForce_Sig)-1]'/Fs;

            
            %Fluo interval to plot
%             IntervalFluoROI_Sig = FluoROI_Sig(stFluo:enFluo);
                IntervalFluoROI_Sig = FluoROI_Sig(stF-(1*100):stF+(3*100));
                figure(ss)
                subplot(211)
                hold on
                plot(stF:enF,ForceSignal(stF:enF),'r')
                subplot(212)
                hold on
                plot(stFluo:enFluo,FluoROI_Sig(stFluo:enFluo),'r')
                
            %time Fluo
%             tF = [0:length(IntervalFluoROI_Sig)-1]'/Fs;
            tFluo = [0-Diff_start:0-Diff_start+length(IntervalFluoROI_Sig)-1]'/Fs;
            
            
         
            
            if PlotMed == 1
                %store dei vettori
                Med_store = [Med_store IntervalFluoROI_Sig];
                               
                %plot
                if iNPeaks==length(StartTime_ForcePar)
                    med   = median(Med_store,2);
                    devst = std(Med_store,[],2);
                    med_der = derivative(med,0.04);
                    med_der = sgolayfilt(med_der,3,11);
                    
                    figure(f)
                    plot(med)
                    hold on
                    plot(med+devst);
                    plot(med-devst);
                    plot(med_der,'r');
                    
                    Med_store = [];
                
                end
            else
                figure(f)
%                 hold on
%                 plot(tF,IntervalForce_Sig*10,'r')
                hold on
                plot(tFluo,IntervalFluoROI_Sig)
                plot(IntervalFluoROI_Sig)
            end
                        
        end
        
        
        ylabel('Fluo')
        xlabel('Time')
%         ylim([0 20])
%         xlim([-0.5 3])
        title(dataGCamp.Info.ROI{id-1})
        
    end
    
    
else
    display('Computation of parameters has to be done before the analysis')
end
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

%%% -> plt mediana
PlotMed = 1;
% PlotMed = 0;
Med_store = [];
%%%

if isfield(dataGCamp,'PeaksPar')
    
    %num signals to analyse: force + fluo ROIs
    NumSig        = length(dataGCamp.fluoROI);
    Fs            = dataGCamp.Info.Fs;
    %
    ForceSignal   = -dataGCamp.fx;
%     ForceSignal   = sgolayfilt(ForceSignal,3,21);
    
    %
    StatusOfInterest   = 3;
    %columns
    Status_Curr_Par    = 2;
    StartTime_Curr_Par = 3;
    Duration_Curr_Par  = 4;
    MaxTime_Curr_Par   = 6;
    
    f = figure('Name','Parameters from Peaks');
    
    %for each force - fluoROI
    for id=2:NumSig
%     for id=5:5
        
        
        %%% force  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %force peaks to take account of
        ind_ForcePeaksToTake = find(dataGCamp.PeaksPar{1,1}(:,Status_Curr_Par) == StatusOfInterest);
        
        %force par (force is fixed) -> the peaks selected can change
        StartTime_ForcePar = dataGCamp.PeaksPar{1,1}(ind_ForcePeaksToTake,StartTime_Curr_Par);
        %force duration
        Duration_ForcePar  = dataGCamp.PeaksPar{1,1}(ind_ForcePeaksToTake,Duration_Curr_Par);
        %time of the max
        MaxTime_ForcePar  = dataGCamp.PeaksPar{1,1}(ind_ForcePeaksToTake,MaxTime_Curr_Par);        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%% fluo   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %fluo sig
        FluoROI_Sig = dataGCamp.fluoROI(:,id-1);       
        %fluo sig -> Low Pass 6 Hz
        FluoROI_Sig = cheb2LPfilt(FluoROI_Sig,6,2,Fs);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        hold on
        subplot(3,3,id-1)
            %%%% Fig whole Sig
            ss  = figure;
            subplot(211)
            plot(ForceSignal)
            subplot(212)
            plot(FluoROI_Sig)
            %%%%%%%%%%%%%%%
        
        %for each peaks
        for iNPeaks=1:length(StartTime_ForcePar)
            
            %start Force
            stF = (StartTime_ForcePar(iNPeaks))*Fs;
            %end Force
            enF = (StartTime_ForcePar(iNPeaks)+Duration_ForcePar(iNPeaks))*Fs-1;
            %max Force
            maxF = (MaxTime_ForcePar (iNPeaks))*Fs;

            
            %to align in agreement with 0 (t start of ForcePeak)
            if strcmp(Alignement,'startForce')
                stF_real = stF;
            elseif strcmp(Alignement,'maxForce')
                stF_real = maxF;        
            end
            
            
            %Force to plot
            IntervalForce_Sig = ForceSignal(stF_real:enF);
            %Force time
            tF     = [0:length(IntervalForce_Sig)-1]; % [points]
            tF_sec = tF/Fs; %[sec]

            
            %Fluo interval to plot
            IntFluoToPlot_sec   = [1 3]'; %[sec]
            IntFluoToPlot       = [IntFluoToPlot_sec(1) IntFluoToPlot_sec(2)]*Fs; % [points]
            IntervalFluoROI_Sig = FluoROI_Sig(     stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))     );
            
            IntervalFluoROI_Sig_Pre_Abs = (FluoROI_Sig(     stF_real-(IntFluoToPlot(1))    :   stF_real)   );
            
            t_FluoROI_Sig_Pre_Abs = [tF_sec(1) - IntFluoToPlot_sec(1)   : 1/Fs: tF_sec(1)] ;
            
                %%%% Fig whole Sig
%                 figure(ss)
%                 subplot(211)
%                 hold on
%                 plot(stF_real:enF,IntervalForce_Sig,'r')
%                 subplot(212)
%                 hold on
%                 IntervalFluoROI_Sig_plot = FluoROI_Sig( stF_real:enF );
%                 plot(stF_real:enF,IntervalFluoROI_Sig_plot,'r')
                %%%%%%%%%%%%%%%
                   
            
            %Fluo time
            tFluo_sec = [tF_sec(1) - IntFluoToPlot_sec(1)   : 1/Fs:   tF_sec(1)+ IntFluoToPlot_sec(2)]'; %[sec]
                                
         
            %Force to plot long
            IntervalForce_Sig_Long = ForceSignal(  stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))   );
            Force_Pre = ForceSignal(  stF_real-(IntFluoToPlot(1)): stF_real);
            Fluo_Pre  = FluoROI_Sig(  stF_real-(IntFluoToPlot(1)): stF_real);
            
            derForcePre = median(derivative(Force_Pre ,0.04))+ 1*std(derivative(Force_Pre ,0.04));
            derFluoPre  = median(derivative(Fluo_Pre, 0.04))+ 1*std(derivative(Force_Pre ,0.04));

            
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
                    plot(tFluo_sec,med)
                    hold on
                    plot(tFluo_sec,med+devst);
                    plot(tFluo_sec,med-devst);
                    plot(tFluo_sec,med_der,'r');
                    
                    Med_store = [];
                
                end
            else
                figure(f)
                hold on
%                 plot(tFluo_sec,(IntervalForce_Sig_Long- Force_offset)*10,'g')
                plot(tFluo_sec,IntervalForce_Sig_Long*10,'g')
                
                derF = derivative(IntervalForce_Sig_Long,0.04);
                plot(tFluo_sec,derF,'m')
                derFF = derivative(sgolayfilt(derF,3,21),0.04);
                plot(tFluo_sec,derFF,'c')
                                
                
                plot(tF_sec,IntervalForce_Sig*10,'b')
%                 derF = derivative(IntervalForce_Sig,0.04);
%                 derF = sgolayfilt(derF,3,11);
%                 plot(tF_sec,derF*10,'k')
%                 der2F = derivative(derF ,0.04);
%                 der2F = sgolayfilt(der2F,3,11);
%                 plot(tF_sec,der2F,'k')  

                    plot(tFluo_sec,IntervalFluoROI_Sig,'r')
                    derFluo = derivative(IntervalFluoROI_Sig,0.04);
                    plot(tFluo_sec,derFluo,'k')
                    derFluoF = derivative(sgolayfilt(derFluo,3,21),0.04);
                    plot(tFluo_sec,derFluoF,'b')
                
%                 [minimum imin] = min(der(1:101))
%                 der = sgolayfilt(der,3,11);

%                 plot( t_FluoROI_Sig_Pre_Abs, IntervalFluoROI_Sig_Pre_Abs,'k')

                               
                ciao=1
                
            end
                        
        end
        
        
        ylabel('Fluo')
        xlabel('Time')
        title(dataGCamp.Info.ROI{id-1})
        
        
%         %%%% find info Peak Fluo
%         Threshold  = std(Sig(RestInterval(1):RestInterval(2)))*3;
        
        
        
        
        
    end
    
    
else
    display('Computation of force parameters has to be done before the analysis')
end
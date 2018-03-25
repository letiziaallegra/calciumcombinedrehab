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
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%Folder Data
User = getenv('username');
UsbPortHD  = 'I';
dirAllDataGCamp = [UsbPortHD,':\LENS\_data_MAT_GCamp\'];
%list animals
list_dirDataGCamp = dir(dirAllDataGCamp);

%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Select what you want to plot (along the days of trials%%%
%%% Parameters
%%% 
ParNameToTake = {'Trial','Status','Start Time','Duration','Max','Max Time','FWHM','AUC','PtPAmp','nSubPeaks'};
ParToTake     = [3 4 5 6 7 8 9 10];
% ParToTake     = [8 9 10 3];
NumPar        = length(ParToTake);
%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% ROIs
%%% ROIs manual
ROINameToTakeManual = {'Motor','Barrel','Occipital','Motor Small','Barrel Small','Occipital Small','Whole'};
ROIToTakeManual     = [1 2 3 4 5 6 7];
% ROIToTakeManual     = [4];
NumROIManual       = length(ROIToTakeManual );
ColorsStartManual   = [0 0 1; 0 1 0; 1 0 0; 1 0 1; 0 1 1; 1 1 0; 0 0 0];
%%%
%%% ROIs MIP based
ROINameToTakeMIP = {'Max MIP', 'Sum MIP'};
ROIToTakeMIP     = [1 2];
NumROIMIP        = length(ROIToTakeMIP);
ColorsStartMIP   = [0 0 0.5; 0 0.5 0; 0.5 0 0];
%%%
%%% SELECT WHICH ROI YOU WANT TO PLOT
%  ROI_Manual                -> SelROI = 0
%               ROI_Centroid -> SelROI = 1
%  ROI_Manual + ROI_Centroid -> SelROI = 2
SelROI = 1;
switch SelROI
    case 0
        %ROI to take
        ROINameToTake = ROINameToTakeManual;
        NumROI        = NumROIManual;
        ROIToTake     = ROIToTakeManual;
        ColorsStart   = ColorsStartManual;
        ListROIs      = zeros(NumROIManual,1);
    case 1
        %ROI to take
        ROINameToTake = ROINameToTakeMIP;
        NumROI        = NumROIMIP;
        ROIToTake     = ROIToTakeMIP;
        ColorsStart   = ColorsStartMIP;
        ListROIs      = ones(NumROIMIP,1);
    case 2
        %ROI to take
        ROINameToTake = [ROINameToTakeManual, ROINameToTakeMIP];
        NumROI        = NumROIManual+NumROIMIP;
        ROIToTake     = [ROIToTakeManual ROIToTakeMIP];
        ListROIs      = [zeros(NumROIManual,1); ones(NumROIMIP,1)];
end

%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% PLOT

%%% PLOT figure
% Plot_Par_Ok = 1;
Plot_Par_Ok = 0;

%%% SubPlot
sPr=2;sPc=5;
% sPr=2;sPc=2;

%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% MEAN

%%% MeanPar
num_tot_day=5;
%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Rehab -> giorni scelti
%%%%
% num_Rehab_day=[]; %-> all days
num_Rehab_day=5;
%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%%
MeanPar = cell(length(ParNameToTake),num_tot_day);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dLGC=3:length(list_dirDataGCamp) %for animals
    
    %num animal
    n_animal = dLGC-2;
    
    NameAnimal =  list_dirDataGCamp(dLGC,1).name;
    
    if strfind(NameAnimal,'control')
        Treat = 0;
    elseif strfind(NameAnimal,'stroke')
        Treat = 1;
        if strfind(NameAnimal,'stroke_BoNT')
            Treat = 2;
        end
    end
    
    %list days of one animal
    dirNameGCamp = dir([dirAllDataGCamp,'\',NameAnimal]);
    
    %Rehab -> scelgo gli ultimi num_Rehab_day giorni
    if Treat == 2 && ~isempty(num_Rehab_day)
            dirNameGCamp  = dirNameGCamp([1:2,length(dirNameGCamp)-num_Rehab_day+1:length(dirNameGCamp)]);
    end
        
    
    for dNGC=3:length(dirNameGCamp) % for days
        
        %num day
        n_day = dNGC-2;
        
        NameDay =  dirNameGCamp(dNGC,1).name;        
        
        Path_GCamp_Day      = [dirAllDataGCamp,'\',NameAnimal,'\',NameDay];
%         NameFile_GCamp_Day  = ['dataMouseGCamp_',NameAnimal,'_',NameDay,'_Par'];  
        NameFile_GCamp_Day  = ['dataMouseGCamp_',NameAnimal,'_',NameDay,'_Par_MIP_Par']; 
        
        if exist([Path_GCamp_Day,'\',NameFile_GCamp_Day,'.mat']) %if exist
            
            %load GCamp
            load([Path_GCamp_Day,'\',NameFile_GCamp_Day,'.mat'])
            
            if isfield(dataGCamp,'ROI_MIP_PeaksPar_Fx_Fluo') %if PeaksPar_Fx_Fluo
                
                %num signals to analyse: force + fluo ROIs                
                switch SelROI
                    case 0
                        %ROI to take
                        NumSig        = size(dataGCamp.PeaksPar_Fx_Fluo,2);
                    case 1
                        %ROI Centroid based
                        NumSig        = size(dataGCamp.ROI_MIP_PeaksPar_Fx_Fluo,2);
                    case 2
                        %ROI + ROI Centroid Based
                        NumSig        = size(dataGCamp.PeaksPar_Fx_Fluo,2) + size(dataGCamp.ROI_MIP_PeaksPar_Fx_Fluo,2);
                end
                
                if NumSig ~= length(ROINameToTake)
                    error('check the ROIs')
                end
               
                if Plot_Par_Ok
                    F_Par = figure('Name',['Parameters from Peaks in animal_',NameAnimal,'_day',NameDay]);
                end
                
                %for each of the Par (one figure for each parameter) -> Par of the force
                for ip=1:NumPar  %if NumSig

                    %info Par
                    Curr_Par      = ParToTake(ip);
                    Curr_Par_Name = ParNameToTake{Curr_Par};
                    
                    
                    MeanForcePar = [];
                    StdForcePar  = [];
                    MeanFluoPar  = [];
                    StdFluoPar   = [];
                                       
                    %for each Par-couple of the FORCE and i-th FLUO ROI
                    for iR=1:NumROI
                        
                        %info ROI
                        Curr_ROI       = ROIToTake(iR);
                        Curr_ROI_Name  = ROINameToTake(ROIToTake(iR));
                        
                        
                        
                        if ListROIs(iR) == 0
                            %Force Par
                            ForcePar = dataGCamp.PeaksPar_Fx_Fluo{1,Curr_ROI}(:,Curr_Par);
                            %Fluo Par
                            FluoPar  = dataGCamp.PeaksPar_Fx_Fluo{2,Curr_ROI}(:,Curr_Par);
                        elseif ListROIs(iR) == 1
                            %Force Par
                            ForcePar = dataGCamp.ROI_MIP_PeaksPar_Fx_Fluo{1,Curr_ROI}(:,Curr_Par);
                            %Fluo Par
                            FluoPar  = dataGCamp.ROI_MIP_PeaksPar_Fx_Fluo{2,Curr_ROI}(:,Curr_Par);
                        end
                        
                        
                        %time start
                        if Curr_Par == 3 || Curr_Par == 6
                            
                            ForcePar_Buf = ForcePar-FluoPar;
                            FluoPar_Buf  = FluoPar-ForcePar;
                            
                            ForcePar = ForcePar_Buf;
                            FluoPar  = FluoPar_Buf;
                            
                        end
                        
                        
                        if Plot_Par_Ok
                            figure(F_Par)
                            hold on
                            subplot(sPr,sPc,ip)
                            ParToPlot = sortrows( [ForcePar, FluoPar],1);
                            plot(ParToPlot(:,1), ParToPlot(:,2), 'Color',ColorsStart(iR,:),'Marker','o','MarkerSize',6,'MarkerFaceColor',ColorsStart(iR,:))
                            %                 scatter(ForcePar, FluoPar,10, ColorsStart(id-1,:),'filled')
                        end
                        
                        %Mean Force Par
                        %                             MeanForcePar = mean(ForcePar);
                        MeanForcePar = median(ForcePar);
                        StdForcePar  = std(ForcePar,[]);
                        
                        %Mean Fluo Par
                        %                             MeanFluoPar  = mean(FluoPar);
                        MeanFluoPar  = median(FluoPar);
                        StdFluoPar   = std(FluoPar,[]);
                        
                        %store mean data
                        MeanPar{Curr_Par,n_day}{n_animal,1}(Curr_ROI,:) = [MeanForcePar StdForcePar MeanFluoPar StdFluoPar];
                        MeanPar{Curr_Par,n_day}{n_animal,2} = [NameAnimal];
                        MeanPar{Curr_Par,n_day}{n_animal,3} = [NameDay];
                        MeanPar{Curr_Par,n_day}{n_animal,4} = [Treat];
                        
                        
                        if Plot_Par_Ok
                            figure(F_Par)
                            subplot(sPr,sPc,ip)
                            title(Curr_Par_Name)
                            ylabel('Fluo')
                            xlabel('Force')
                            if ip==NumPar %legend at the end
                                legend(ROINameToTake{ROIToTake},'Location','EastOutside')
                            end
                        end
                        
                        
                    end
                    
                end %end if NumSig
                
            else
                display('Computation of parameters has to be done before the analysis')
            end %end if PeaksPar_Fx_Fluo
            
        else
            display(['missing DataGCamp_Par in,' Path_GCamp_Day]);
        end      %end if exist
        
    end %end for days
    
end %end for animals



if ~isempty(MeanPar)  
    
    %%%%%%
    Curnt = clock;
    H = Curnt(4);
    M = Curnt(5);
    %es. 30-Jan-2014
    DateString = date;
    
    Filename = ['dataMouseGCamp_Median','_SelROI_',num2str(SelROI),'_',DateString,'_',num2str(H),'-',num2str(M)];
    save(Filename,'MeanPar')
    %%%%%%
end

























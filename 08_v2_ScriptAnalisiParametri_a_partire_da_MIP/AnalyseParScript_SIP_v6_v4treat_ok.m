%Script To analyse Parameters from Peaks
%AnalyseParScript
%Analyis performed on the parameters of dataGCamp.PeaksPar_Fx_Fluo
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
%                   column 10 -> NumberOfRealPeak inside the peak     %non più il num of subpeaks inside the main peak (numPks)
%                   column 11 -> SlopeInitial (slope between onset _ max value)
%                   column 12 -> Time to Peak ( [max value peak point - peak start point] )
%                   column 13 -> Time Distance bw first movement - onset  
%                   column 14 -> Time Distance bw first movement - time of max 
% close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%%
%Folder Data
User = getenv('username');
UsbPortHD  = 'M';
% dirAllDataGCamp = [UsbPortHD,':\LENS\_data_MAT_GCamp\'];
dirAllDataGCamp = ['C:\Users\asus\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\_data_MAT_GCamp_Store_Analysis_Par'];
%list animals
list_dirDataGCamp = dir(dirAllDataGCamp);

%%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Select what you want to plot (along the days of trials%%%
%%% Parameters
%%%
ParNameToTake = {'Trial','Status','OnsetT','Duration','Max','MaxT','FWHM','AUC','PtPAmp','nRealPeaks','Slope','TimeToPeak','FirstMovT-OnsetT','FirstMovT-MaxT'};
ParToTake     = [3 4 5 6 7 8 9 10 11 12 13 14];
%ParToTake     = [5];
NumPar        = length(ParToTake);
%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% ROIs
%%% ROIs manual
% ROINameToTakeManual = {'Motor','Barrel','Occipital','Motor Small','Barrel Small','Occipital Small','Whole'};
% ROIToTakeManual     = [1 2 3 4 5 6 7];
ROINameToTakeManual = {'Whole'};
ROIToTakeManual     = [1];
% ROIToTakeManual     = [4];
NumROIManual       = length(ROIToTakeManual );
ColorsStartManual   = [0 0 1; 0 1 0; 1 0 0; 1 0 1; 0 1 1; 1 1 0; 0 0 0];
%%%
%%% ROIs MIP based
ROINameToTakeMIP = {'MIP','SIP_Rostral','SIP_Caudal'};
ROIToTakeMIP     = [1 2 3];
NumROIMIP        = length(ROIToTakeMIP);
ColorsStartMIP   = [0 0 0.5; 0 0.5 0; 0.5 0 0];
%%%
%%% SELECT WHICH ROI YOU WANT TO PLOT
%  ROI_Manual                -> SelROI = 0
%               ROI_MIP      -> SelROI = 1
%  ROI_Manual + ROI_MIP      -> SelROI = 2
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
%%% MEAN or MEDIAN
chStat_MEAN_MEDIAN = 'MEAN';
% chStat_MEAN_MEDIAN = 'MEDIAN';

%%% MeanPar
% num_Rehab_day=[]; %-> all day
num_tot_day=5; %-> in this wat I will take the first and the last of of rehab
%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Rehab -> giorni scelti
% num_Rehab_day=[]; %-> all days of Rehab (only Rehab are selected)
num_Rehab_day=5; %-> in this way I will take the first and the last of of rehab
%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%%
MeanPar = cell(length(ParNameToTake),num_tot_day);

%num animal
n_animal_count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dLGC=3:length(list_dirDataGCamp) %for animals
    
    n_animal_count = n_animal_count+1;
    Treat_old = [];
       
    NameAnimal =  list_dirDataGCamp(dLGC,1).name;
    
    
    %%%% select animal
    
%     if strfind(NameAnimal,'control')
%         Treat = 0;
%         
%     elseif strfind(NameAnimal,'stroke')
%         Treat = 1;
%         
%         if strfind(NameAnimal,'stroke_BoNT')
%             if isempty(num_Rehab_day)
%                 Treat   = [2];
%                 
%             else
%                 Treat   = [2 3]; %-> LAST WEEK and FIRST WEEK
%                 
%             end
%             
%         elseif strfind(NameAnimal,'stroke_Rehab')
%             if isempty(num_Rehab_day)
%                 Treat   = [21];
%                 
%             else
%                 Treat   = [21 31]; %-> LAST WEEK and FIRST WEEK
%                 
%             end
%             
%         elseif strfind(NameAnimal,'stroke_STIM')
%             
%                 Treat   = [214]; %stroke reabilitato 4 settimane STIM
%                 %Treat   = [21 31]; %-> LAST WEEK and FIRST WEEK
%                 
%                       
%         end
%     end
    if contains(NameAnimal,'control')
        Treat = 0;
        
    elseif contains(NameAnimal,'stroke')
        Treat = 1;
        
        if contains(NameAnimal,'stroke_BoNT')
            if isempty(num_Rehab_day)
                Treat   = [2];
                
            else
                Treat   = [2 3]; %-> LAST WEEK and FIRST WEEK
                
            end
            
        elseif contains(NameAnimal,'stroke_Rehab')
            if isempty(num_Rehab_day)
                Treat   = [21];
                
            else
                Treat   = [21 31]; %-> LAST WEEK and FIRST WEEK
                
            end
            
        elseif contains(NameAnimal,'stroke_STIM')
            
                Treat   = [214]; %stroke reabilitato 4 settimane STIM
                %Treat   = [21 31]; %-> LAST WEEK and FIRST WEEK
                
                      
        end
    elseif contains(NameAnimal,'TOX')
        
        Treat   = [5];
    end
    %%%%
    
    
    %list days of one animal
    PathFile_GCamp = [dirAllDataGCamp,'\',NameAnimal];
    dirNameGCamp   = dir(PathFile_GCamp);
    
    %     %Rehab -> scelgo gli ultimi num_Rehab_day giorni
    %     if Treat == 2 && ~isempty(num_Rehab_day)
    %             last_five_days = [16, 17, 18, 19, 20];
    %             dirNameGCamp  = dirNameGCamp([1:2,length(dirNameGCamp)-num_Rehab_day+1:length(dirNameGCamp)]);
    %     end
  
    
    
    for dNGC=3:length(dirNameGCamp) % for days
        
        if (isempty(num_Rehab_day) && (  (Treat == 0) || (Treat == 1)) )  
            break %no rehab are excluded
        end
        
        %num day
        n_day = dNGC-2;
        
        NameFile_GCamp_Day =  dirNameGCamp(dNGC,1).name;
                
        if contains(NameFile_GCamp_Day,'_Par_MIPSIP_Par') %if exist
            
            %load GCamp
            load([PathFile_GCamp,'\',NameFile_GCamp_Day])
            DayFile = dataGCamp.Info.Date;
            
            for itr=1:length(Treat) %for length Treat
                
                EndDay = 20;
                if length(Treat)>1 %FIRST and LAST REHAB WEEK
                    
                    
                    %%% Robot + BoNT
                    if Treat(itr) == 2 %LAST REHAB WEEK
                        
                        ListDayToUse = [EndDay-num_Rehab_day+1:EndDay];
                        
                        if ismember(str2num(DayFile),ListDayToUse)
                            Tk = 1;
                        else
                            Tk = 0;
                        end
                        RehabNameYN  = '';
                        
                    elseif Treat(itr) == 3 %First REHAB WEEK
                        
                        ListDayToUse = 1:5;
                        
                        if ismember(str2num(DayFile),ListDayToUse)
                            Tk = 1;
                        else
                            Tk = 0;
                        end
                        RehabNameYN  = '';
                    %%%
                    %%% Robot 
                    elseif Treat(itr) == 21 || Treat(itr) == 214 %LAST REHAB WEEK
                       
                        ListDayToUse = [EndDay-num_Rehab_day+1:EndDay];
                        
                        if ismember(str2num(DayFile),ListDayToUse)
                            Tk = 1;
                        else
                            Tk = 0;
                        end
                        RehabNameYN  = '';
                        
                    elseif Treat(itr) == 31 %First REHAB WEEK
                        
                        ListDayToUse = 1:5;
                        
                        if ismember(str2num(DayFile),ListDayToUse)
                            Tk = 1;
                        else
                            Tk = 0;
                        end
                        RehabNameYN  = '';
                        
                    end
                    
                else %other
                    
                    if Treat(itr) == 2 && isempty(num_Rehab_day)
                        Tk = 1;
                        ListDayToUse = 1:EndDay;
                        RehabNameYN  = 'Rehab_';
                        
                           
                    else
                        
                        Tk = 1;
                        ListDayToUse = 1:5;
                        RehabNameYN  = '';
                    end
                end
                
                
                if Tk %if check
                    
                    if isfield(dataGCamp,'ROI_MIP_SIP_PeaksPar_Fx_Fluo') %if PeaksPar_Fx_Fluo
                        
                        %num signals to analyse: force + fluo ROIs
                        switch SelROI
                            case 0
                                %ROI to take
                                NumSig        = size(dataGCamp.PeaksPar_Fx_Fluo,2);
                            case 1
                                %ROI Centroid based
                                NumSig        = size(dataGCamp.ROI_MIP_SIP_PeaksPar_Fx_Fluo,2);
                            case 2
                                %ROI + ROI Centroid Based
                                NumSig        = size(dataGCamp.PeaksPar_Fx_Fluo,2) + size(dataGCamp.ROI_MIP_SIP_PeaksPar_Fx_Fluo,2);
                        end
                        
                        if NumSig ~= length(ROINameToTake)
                            error('check the ROIs')
                        end
                        
                        if Plot_Par_Ok
                            F_Par = figure('Name',['Parameters from Peaks in animal_',NameAnimal,'_day',NameFile_GCamp_Day]);
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
                                    
                                    if ~isempty(dataGCamp.ROI_MIP_SIP_PeaksPar_Fx_Fluo{1,Curr_ROI})
                                        
                                        %Force Par
                                        ForcePar = dataGCamp.ROI_MIP_SIP_PeaksPar_Fx_Fluo{1,Curr_ROI}(:,Curr_Par);
                                        %Fluo Par
                                        FluoPar  = dataGCamp.ROI_MIP_SIP_PeaksPar_Fx_Fluo{2,Curr_ROI}(:,Curr_Par);
                                        
                                    else
                                        ForcePar = NaN;
                                        FluoPar  = NaN;
                                        
                                    end
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
                                
                                %Mean Par
                                if strcmp(chStat_MEAN_MEDIAN,'MEAN')
                                    %FORCE
                                    MeanForcePar = mean(ForcePar);
                                    %FLUO
                                    MeanFluoPar  = mean(FluoPar);
                                    %Median Par
                                elseif strcmp(chStat_MEAN_MEDIAN,'MEDIAN')
                                    %FORCE
                                    MeanForcePar = median(ForcePar);
                                    %FLUO
                                    MeanFluoPar  = median(FluoPar);
                                end
                                StdForcePar  = std(ForcePar,[]);
                                StdFluoPar   = std(FluoPar,[]);
                                
                                %update num animal
                                if Treat== 2 | Treat == 3 | Treat == 21 | Treat == 31 | Treat == 214
                                    
                                    if isempty(Treat_old)
                                        Treat_old = Treat(itr);
                                        
                                    elseif Treat_old ~= Treat(itr)
                                        n_animal_count = n_animal_count +1;
                                        Treat_old = Treat(itr);
                                        
                                    end
                                end
                                
                                %store mean data
                                Day_int = str2num(DayFile);
                                i_mP = find(ListDayToUse==Day_int);
                                
                                %store mean data
                                MeanPar{1,i_mP}                                       = Day_int;
                                MeanPar{Curr_Par,i_mP}{n_animal_count,1}(Curr_ROI,:)  = [MeanForcePar StdForcePar MeanFluoPar StdFluoPar];
                                MeanPar{Curr_Par,i_mP}{n_animal_count,2}              = [NameAnimal];
                                MeanPar{Curr_Par,i_mP}{n_animal_count,3}              = [NameFile_GCamp_Day];
                                MeanPar{Curr_Par,i_mP}{n_animal_count,4}              = [Treat(itr)];
                                
                                
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
                    display('Day from Rehab not to take (check Tk==0)')
                end %end check
                
            end %end for length Treat

            
        elseif strfind(NameFile_GCamp_Day,'.ini') %if exist
            display([NameFile_GCamp_Day, ' is not to be taken']);
            
        else
            display(['missing DataGCamp_Par in,' Path_GCamp_Day]);
        end      %end if exist
        
    end %end for days
    
end %end for animals



% %check of all the animals are present very day (even if in that day the
% %animal did not performed the trial session (I will put NaN)
% for chDay=1:size(MeanPar,2)
%
%     ListAnimalMeanParDay = MeanPar
%
%     for lg=3:length(list_dirDataGCamp)
%         %num animal
%         n_animal = dLGC-2;
%         NameAnimal =  list_dirDataGCamp(lg,1).name;
%
%         if ~sum(strcmp(MeanPar{3,chDay}(:,2),NameAnimal))
%            ciao = 1
%         end
%
%     end
%
% end




if ~isempty(MeanPar)
    
    %%%%%%
    Curnt = clock;
    H = Curnt(4);
    M = Curnt(5);
    %es. 30-Jan-2014
    DateString = date;
    
    Filename = ['dataMouseGCamp_',RehabNameYN,chStat_MEAN_MEDIAN,'_SelROI_',num2str(SelROI),'_',DateString,'_',num2str(H),'-',num2str(M)];
    save(Filename,'MeanPar')
    %%%%%%
end

























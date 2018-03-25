%Script To analyse Parameters from ROBOT -> 4 groups
%AnalyseParScript
%Analyis performed on the parameters of dataGCamp.PeaksPar_Fx_Fluo
%        each column of dataGCamp.PeaksPar_Fx_Robot{:,1} refers to a peak:

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%%
%Folder Data
User = getenv('username');
UsbPortHD = 'F';

Username = getenv('username');
% dirAllDataGCamp = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_data_MAT_GCamp_Store_Rehab\'];
dirAllDataGCamp = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\_data_MAT_GCamp_Store_Analysis_Par'];
%list animals
list_dirDataGCamp = dir(dirAllDataGCamp);

%%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Select what you want to plot (along the days of trials%%%
%%% Parameters
%%%
ParNameToTake = {'Trial','Target Time','Force','Sub Mov','Attempts','Max Force'};
ParToTake     = [2 3 4 5 6];
NumPar        = length(ParToTake);

%%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Stored in
ROINameToTake = {'Robot Par'};
ROIToTake     = [1];
NumROI        = length(ROIToTake);
ColorsStart   = [0 0 1; 0 1 0; 1 0 0; 1 0 1; 0 1 1; 1 1 0; 0 0 0];


%%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% PLOT5
%%% PLOT figure
% Plot_Par_Ok = 1;
Plot_Par_Ok = 0;

%%% SubPlot
sPr=2;sPc=2;
%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% MEAN or MEDIAN
chStat_MEAN_MEDIAN = 'MEAN';
% chStat_MEAN_MEDIAN = 'MEDIAN';

%%% MeanPar
% num_Rehab_day=[]; %-> all day
num_tot_day=5;
%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Rehab -> giorni scelti
num_Rehab_day=[]; %-> all days of Rehab (only stroke_BoNT or stroke_Rehab are selected)
% num_Rehab_day=5;
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
        elseif strfind(NameAnimal,'stroke_Rehab')
            Treat = 3;
        end
        
    end
    
    %list days of one animal
    PathFile_GCamp = [dirAllDataGCamp,'\',NameAnimal];
    dirNameGCamp = dir(PathFile_GCamp);
    
    %     %Rehab -> scelgo gli ultimi num_Rehab_day giorni
    %     if Treat == 2 && ~isempty(num_Rehab_day)
    %         dirNameGCamp  = dirNameGCamp([1:2,length(dirNameGCamp)-num_Rehab_day+1:length(dirNameGCamp)]);
    %     end
    
    for dNGC=3:length(dirNameGCamp) % for days
        
        %caso di tutti  gionri (solo stroke_BoNT o stroke_Rehab)
        if isempty(num_Rehab_day) && ((Treat == 0) || (Treat == 1))
            break
        end

        %num day
        n_day = dNGC-2;
        
        NameFile_GCamp_Day =  dirNameGCamp(dNGC,1).name;
        
        
        if strfind(NameFile_GCamp_Day,'_Par_MIPSIP_Par') %if exist
            
            %load GCamp
            load([PathFile_GCamp,'\',NameFile_GCamp_Day])
            DayFile = dataGCamp.Info.Date;
            
            EndDay = 20;
            if (Treat == 2 && ~isempty(num_Rehab_day)) ||...
               (Treat == 3 && ~isempty(num_Rehab_day))
                                
                ListDayToUse = [EndDay-num_Rehab_day+1:EndDay];
                
                if ismember(str2num(DayFile),ListDayToUse)
                    Tk = 1;
                else
                    Tk = 0;
                end
                RehabNameYN  = '';
                
            elseif (Treat == 2 && isempty(num_Rehab_day)) ||...
                   (Treat == 3 && isempty(num_Rehab_day))
                Tk = 1;
                ListDayToUse = 1:EndDay;
                RehabNameYN  = 'Rehab_';
                
            else
                Tk = 1;
                ListDayToUse = 1:5;
                RehabNameYN  = '';
                
            end
            
            
            
            if Tk %if check
                
                                
                if isfield(dataGCamp,'PeaksPar_Fx_Robot') %if PeaksPar_Fx
                    
                    %num signals to analyse: force + fluo ROIs
                    NumSig        = 1;
                    
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
                        
                        %in this case there is just one NumROI 
                        for iR=1:NumROI
                            
                            %info ROI
                            Curr_ROI       = ROIToTake(iR);
                            Curr_ROI_Name  = ROINameToTake(ROIToTake(iR));
                            
                            if ~isempty(dataGCamp.PeaksPar_Fx_Robot)
                                
                                %Par From ForcePeak
                                tr                = dataGCamp.PeaksPar_Fx_Robot(:,1);
                                Par_FromForcePeak = dataGCamp.PeaksPar_Fx_Robot(:,Curr_Par);
                                
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                if Plot_Par_Ok
                                    figure(F_Par)
                                    hold on
                                    subplot(sPr,sPc,ip)
                                    plot(tr, Par_FromForcePeak, 'Color',ColorsStart(iR,:),'Marker','o','MarkerSize',6,'MarkerFaceColor',ColorsStart(iR,:))
                                    %                 scatter(ForcePar, FluoPar,10, ColorsStart(id-1,:),'filled')
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                %Mean Force Par
                                if strcmp(chStat_MEAN_MEDIAN,'MEAN')
                                    MeanForcePar = nanmean(Par_FromForcePeak);
                                elseif strcmp(chStat_MEAN_MEDIAN,'MEDIAN')
                                    MeanForcePar = nanmedian(Par_FromForcePeak);
                                end
                                StdForcePar  = nanstd(Par_FromForcePeak,[]);
                                
                                
                                %store mean data                                
                                Day_int = str2num(DayFile);
                                i_mP = find(ListDayToUse==Day_int);
                                
                                MeanPar{1,i_mP}                                 = Day_int;
                                MeanPar{Curr_Par,i_mP }{n_animal,1}(Curr_ROI,:) = [MeanForcePar StdForcePar];
                                MeanPar{Curr_Par,i_mP }{n_animal,2}             = [NameAnimal];
                                MeanPar{Curr_Par,i_mP }{n_animal,3}             = [NameFile_GCamp_Day];
                                MeanPar{Curr_Par,i_mP }{n_animal,4}             = [Treat];
                                
                            end
                            
                            if Plot_Par_Ok
                                figure(F_Par)
                                subplot(sPr,sPc,ip)
                                ylabel(Curr_Par_Name)
                                xlabel('trials')
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
            
        elseif strfind(NameFile_GCamp_Day,'.ini') %if exist
            display([NameFile_GCamp_Day, ' is not to be taken']);
        else
            display(['missing DataGCamp_Par in,' PathFile_GCamp]);
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
    
    Filename = ['dataMouseGCamp_RobotPar_',RehabNameYN,chStat_MEAN_MEDIAN,'_',DateString,'_',num2str(H),'-',num2str(M)];
    save(Filename,'MeanPar')
    %%%%%%
end

























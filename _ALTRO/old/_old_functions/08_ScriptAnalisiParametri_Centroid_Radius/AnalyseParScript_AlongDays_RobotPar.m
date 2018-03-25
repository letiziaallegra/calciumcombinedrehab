%Script To analyse Parameters from ROBOT
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
%Choice of the animal and trial day
UsbPortHD = 'I';

Username = getenv('username');
dirAllDataGCamp = [UsbPortHD,':\LENS\_data_MAT_GCamp\'];
%list animals
list_dirDataGCamp = dir(dirAllDataGCamp);

%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Select what you want to plot (along the days of trials%%%
%%% Parameters
%%% 
ParNameToTake = {'Trial','Target Time','Force','Sub Mov','Attempts'};
ParToTake     = [2 3 4 5];
NumPar        = length(ParToTake);

%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Stored in
ROINameToTake = {'Robot Par'};
ROIToTake     = [1];
NumROI        = length(ROIToTake);
ColorsStart   = [0 0 1; 0 1 0; 1 0 0; 1 0 1; 0 1 1; 1 1 0; 0 0 0];


%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% PLOT
%%% PLOT figure
% Plot_Par_Ok = 1;
Plot_Par_Ok = 0;

%%% SubPlot
sPr=2;sPc=2;
%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% MEAN

%%% MeanPar
num_tot_day=5;
%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Rehab -> giorni scelti
%%%%
num_Rehab_day=[]; %-> all days
% num_Rehab_day=5;
%%%

% %%% %%% %%% %%% %%% %%% %%% %%% %%%
% MeanPar = cell(length(ParNameToTake),num_tot_day);


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
        NameFile_GCamp_Day  = ['dataMouseGCamp_',NameAnimal,'_',NameDay,'_Par']; 
        
        if exist([Path_GCamp_Day,'\',NameFile_GCamp_Day,'.mat']) %if exist
            
            %load GCamp
            load([Path_GCamp_Day,'\',NameFile_GCamp_Day,'.mat'])
            
            %Day
            Day_str = dataGCamp.Info.Date;
            Day_int = str2num(dataGCamp.Info.Date);
            
            if isfield(dataGCamp,'PeaksPar_Fx_Robot') %if PeaksPar_Fx_Fluo
                
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
                    MeanFluoPar  = [];
                    StdFluoPar   = [];
                                       
                    %for each Par-couple of the FORCE and i-th FLUO ROI
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
                            MeanForcePar = nanmean(Par_FromForcePeak);
%                             MeanForcePar = nanmedian(Par_FromForcePeak);
                            StdForcePar  = nanstd(Par_FromForcePeak,[]);
                            
                            
                            %store mean data
                            MeanPar{1,Day_int}                                = Day_int;
                            MeanPar{Curr_Par,Day_int}{n_animal,1}(Curr_ROI,:) = [MeanForcePar StdForcePar];
                            MeanPar{Curr_Par,Day_int}{n_animal,2} = [NameAnimal];
                            MeanPar{Curr_Par,Day_int}{n_animal,3} = [NameDay];
                            MeanPar{Curr_Par,Day_int}{n_animal,4} = [Treat];
                                                       
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
    
    Filename = ['dataMouseGCamp_RobotPar_Median','_',DateString,'_',num2str(H),'-',num2str(M)];
    save(Filename,'MeanPar')
    %%%%%%
end

























Username = getenv('username');
dirAllDataGCamp = ['C:\Users\',Username,'\Dropbox\LENS\ELABORAZIONE DATA\_data_MAT_GCamp'];
%list animals
list_dirDataGCamp = dir(dirAllDataGCamp);


for dLGC=3:length(list_dirDataGCamp) %for animals
    
    %num animal
    n_animal = dLGC-2;
    
    NameAnimal =  list_dirDataGCamp(dLGC,1).name;
    
    
    %list days of one animal
    dirNameGCamp = dir([dirAllDataGCamp,'\',NameAnimal]);
    
    for dNGC=3:length(dirNameGCamp) % for days
        
       %num day
        n_day = dNGC-2;
        
        NameDay =  dirNameGCamp(dNGC,1).name;        
        
        Path_GCamp_Day      = [dirAllDataGCamp,'\',NameAnimal,'\',NameDay];
        NameFile_GCamp_Day  = ['dataMouseGCamp_',NameAnimal,'_',NameDay]; 
        NameFile_GCamp_Day  = ['dataMouseGCamp_',NameAnimal,'_',NameDay,'_Par'];  
        
        %load GCamp
        load([Path_GCamp_Day,'\',NameFile_GCamp_Day,'.mat'])
        
        
        %%%%% AGGIUNGI QUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t_target  = [];
        ForceMean = [];
        SubMov    = [];
        Attempts  = [];
        
        %%%
        t_target =  [dataGCamp.InfoTrial.Trials_Start_End_St3(:,2) - dataGCamp.InfoTrial.Trials_Start_End_St3(:,1)] /dataGCamp.Info.Fs;
        
        NumTrials = dataGCamp.InfoTrial.NumTrials;
        for ntf=1:NumTrials
            %find total num force peaks over threshold in that trial
            i_F_OT_Tr = find(dataGCamp.PeaksForceNoCleanPeakPar_OverNoiseThreshold(:,Trial_Curr_Par) == ntf);
            %find total num force peaks over threshold in that trial in status 3
            i_F_OT_Tr_St3  = find(dataGCamp.PeaksForceNoCleanPeakPar_OverNoiseThreshold(i_F_OT_Tr,Status_Curr_Par) == StatusOfInterest);
            %total num force peaks
            NumTotForcePeaksOT = length(i_F_OT_Tr_St3);
            
            
            %find trial for peaks ok
            i_F_Tr = find(dataGCamp.PeaksForceNoCleanPeakPar(:,Trial_Curr_Par) == ntf);
            if ~isempty(i_F_Tr)
                for nFtr=1:length(i_F_Tr)
                    %find status in trial for peaks ok
                    i_F_Tr_St3 = find(dataGCamp.PeaksForceNoCleanPeakPar(i_F_Tr,Status_Curr_Par) == StatusOfInterest);
                    
                    %%%
                    ForceMean(ntf,1) =  mean(dataGCamp.PeaksForceNoCleanPeakPar(i_F_Tr(i_F_Tr_St3),5));
                    
                    %%%
                    SubMov(ntf,1)    =  length(dataGCamp.PeaksForceNoCleanPeakPar(i_F_Tr(i_F_Tr_St3),5));
                    
                    %%%
                    Attempts(ntf,1)  =  NumTotForcePeaksOT-SubMov(ntf,1);
                    
                end
            else
                %%%
                ForceMean(ntf,1) =  NaN;
                
                %%%
                SubMov(ntf,1)    =  NaN;
                
                %%%
                Attempts(ntf,1)  =  NaN;
            end
        end
        rmfield(dataGCamp,'PeaksPar_Fx_Robot');
        dataGCamp.PeaksPar_Fx_Robot = [[1:NumTrials]', t_target, ForceMean, SubMov, Attempts];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         IndexStatusMENO1 = find(dataGCamp.status==-1);
%         %punto di start di ogni prova
%         WholeTrials_Index = [IndexStatusMENO1(find(diff(IndexStatusMENO1)>1)); IndexStatusMENO1(end)];
%         % number of completed trials
%         NumTrials = length(WholeTrials_Index);
%         dataGCamp.InfoTrial.NumTrials = NumTrials;
%         
%         TrialsVector = zeros(length(dataGCamp.status),1);
%         WholeTrials_Index = [1; WholeTrials_Index];
%         for nt=1:length(WholeTrials_Index)
%             if nt<length(WholeTrials_Index)
%                 TrialsVector(WholeTrials_Index(nt): WholeTrials_Index(nt+1)-1) =  nt;
%             else
%                 TrialsVector(WholeTrials_Index(nt): end) =  nt;
%             end
%         end
%         % division of the array based on the trials (1-> trial 1; 2-> trial 2 etc.)
%         dataGCamp.InfoTrial.TrialsVector  = TrialsVector;
%         
%         IndexStatus3 = find(dataGCamp.status == 3);
%         %inizio dei trials (status3) nel vettore IndexTrials3
%         Trials_In = [1; find(diff(IndexStatus3)>1)+1];
%         %fine dei trials (status3) nel vettore IndexTrials3
%         Trials_End = [Trials_In(2:end)-1; length(IndexStatus3)];
%         
%         % start and end points of the status 3 of the trials
%         dataGCamp.InfoTrial.Trials_Start_End_St3  = [IndexStatus3(Trials_In(1:NumTrials)) IndexStatus3(Trials_End(1:NumTrials))];
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         dataGCamp.Info.Fs_FluoImages = 25;
%         if strcmp(NameFile_GCamp_Day,'dataMouseGCamp_GCaMP14_stroke_02')
%             dataGCamp.Info.Fs_FluoImages = 100;
%         end
%        
%         if isfield(dataGCamp.Info,'IntervalToFluoForceMean')
%             Int = dataGCamp.Info.IntervalToFluoForceMean;
%         else
%             Int = dataGCamp.Info.IntervalToFluoForceMean_PointFluoFreq;
%         end
%         
%         dataGCamp.Info.IntervalToFluoForceMean_PointFluoFreq  = Int;
%         dataGCamp.Info.IntervalToFluoForceMean_sec            = Int/dataGCamp.Info.Fs_FluoImages;
%         dataGCamp.Info.IntervalToFluoForceMean_PoinForceFreq  = (Int/dataGCamp.Info.Fs_FluoImages)*dataGCamp.Info.Fs;
%         
%         if isfield(dataGCamp.Info,'IntervalToFluoForceMean')
%             dataGCamp.Info = rmfield(dataGCamp.Info,'IntervalToFluoForceMean');
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %save GCamp
        save([Path_GCamp_Day,'\',NameFile_GCamp_Day,'.mat'],'dataGCamp')
        
        
    end
end
        
        
        
        
   
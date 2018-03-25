%data animal load
function [folderTASK_FLUO folderTASK_ForceFileName] = fun_loadGCaMPanimal(UsbPort,Animal_Name,TrialDay)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   control    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% GCaMP4 %%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Animal_Name,'GCaMP4_control')
    
     switch TrialDay
         case '01'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150525\'];
             folderTASK_FLUO = [folderTASK,'150525_GCaMP4_trialx15_8bits - ok dal trialx15_0274\ok'];
             folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
         case '02'
         case '03'
         case '04'
         case '05'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150529\'];
             folderTASK_FLUO = [folderTASK,'150529_GCaMP4_trialx15_8bits'];
             folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
     end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   stroke    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% GCaMP9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP9_stroke')
    
    switch TrialDay
        case '01'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150727\'];
            folderTASK_FLUO = [folderTASK,'150727_GCaMP9_trailx15_Surperficial_8bits\ok'];
            folderTASK_ForceFileName = [folderTASK,'trialx15.txt'];
        case '02'
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150729\'];
            folderTASK_FLUO = [folderTASK,'150729_GCaMP9_trialx15_8bits'];
            folderTASK_ForceFileName = [folderTASK,'trialx15.txt'];
        case '04'
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150731\'];
            folderTASK_FLUO = [folderTASK,'150731_GCaMP9_trailx15_8bits'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
    end
end
%%%% 

           
end
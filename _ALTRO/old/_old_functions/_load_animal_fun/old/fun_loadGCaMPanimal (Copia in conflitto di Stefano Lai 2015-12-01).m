%data animal load
function [folderTASK_FLUO folderTASK_ForceFileName] = fun_loadGCaMPanimal(UsbPort,Animal_Name,TrialDay)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   control    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% GCaMP3 %%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Animal_Name,'GCaMP3_control')  
    
    switch TrialDay
        case '01'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150511\'];
            folderTASK_FLUO = [folderTASK,'sequenza'];
            folderTASK_ForceFileName = [folderTASK,'trial4_sync.txt'];
        case '02'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150512\'];
            folderTASK_FLUO = [folderTASK,'sequenza'];
            folderTASK_ForceFileName = [folderTASK,'trial complete x15_TMP_sync.txt'];
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150513\'];
            folderTASK_FLUO = [folderTASK,'sequenza'];
            folderTASK_ForceFileName = [folderTASK,'TC15_F_sync.txt'];
        case '04'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150514\'];
            folderTASK_FLUO = [folderTASK,'sequenza'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150515\'];
            folderTASK_FLUO = [folderTASK,'sequenza'];
            folderTASK_ForceFileName = [folderTASK,'trial15_ok_sync.txt'];
    end   
    
    
%%%%% GCaMP4 %%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP4_control')
    
     switch TrialDay
         case '01'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150525\'];
             folderTASK_FLUO = [folderTASK,'150525_GCaMP4_trialx15_8bits - ok dal trialx15_0274\ok'];
             folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
         case '02'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150526\'];
             folderTASK_FLUO = [folderTASK,'trialx15_8bits'];
             folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
         case '03'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150527\'];
             folderTASK_FLUO = [folderTASK,'trialx15con resistenza'];
             folderTASK_ForceFileName = [folderTASK,'trialx15completoOK_sync.txt'];
         case '04'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150528\'];
             folderTASK_FLUO = [folderTASK,'trialx15'];
             folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
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
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '02'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150728\'];
            folderTASK_FLUO = [folderTASK,'sequenza'];
            folderTASK_ForceFileName = [folderTASK,'trialx15vero_sync.txt'];
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150729\'];
            folderTASK_FLUO = [folderTASK,'150729_GCaMP9_trialx15_8bits'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '04'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150730\'];
            folderTASK_FLUO = [folderTASK,'sequenza'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150731\'];
            folderTASK_FLUO = [folderTASK,'150731_GCaMP9_trailx15_8bits'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
    end
%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   stroke    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% GCaMP10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP10_stroke')
    
     switch TrialDay
         case '01'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150803\'];
             folderTASK_FLUO = [folderTASK,'sequenza'];
             folderTASK_ForceFileName = [folderTASK,'2_primi4-5trials_poi si e sfilata la zampa_sync.txt'];
         case '02'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150804\'];
%              folderTASK_FLUO = [folderTASK,'sequenza2'];
%              folderTASK_ForceFileName = [folderTASK,'2_trialx6_deep_sync.txt'];
             folderTASK_FLUO = [folderTASK,'sequenza1'];
             folderTASK_ForceFileName = [folderTASK,'1_trialx8+1_sync.txt'];
         case '03'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150805\'];
             folderTASK_FLUO = [folderTASK,'sequenza'];
             folderTASK_ForceFileName = [folderTASK,'1_trialx12_sync.txt'];
         case '04'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150806\'];
             folderTASK_FLUO = [folderTASK,'sequenza'];
             folderTASK_ForceFileName = [folderTASK,'2_trialx9_deep_sync.txt'];
         case '05'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150807\'];
             folderTASK_FLUO = [folderTASK,'sequenza'];
             folderTASK_ForceFileName = [folderTASK,'1_trialx8_sync.txt'];
     end
     
     
     
%%%% GCaMP11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP11_stroke')
         
         switch TrialDay
             case '01'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150808\'];
                 folderTASK_FLUO = [folderTASK,'sequenza'];
                 folderTASK_ForceFileName = [folderTASK,'1_trialx13_deep_sync.txt'];
             case '02'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150809\'];
                 folderTASK_FLUO = [folderTASK,'sequenza'];
                 folderTASK_ForceFileName = [folderTASK,'1_trialx15_1mmdeep_sync.txt'];
             case '03'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150810\'];
                 folderTASK_FLUO = [folderTASK,'sequenza'];
                 folderTASK_ForceFileName = [folderTASK,'2_trialx10_1mmdeep_sync.txt'];
             case '04'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150811\'];
                 folderTASK_FLUO = [folderTASK,'sequenza'];
                 folderTASK_ForceFileName = [folderTASK,'1_trialx10_1mmdeep_sync.txt'];
             case '05'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150812\'];
                 folderTASK_FLUO = [folderTASK,'sequenza'];
                 folderTASK_ForceFileName = [folderTASK,'2_trialx7_1mmdeep_sync.txt'];
         end
   
         
%%%% GCaMP14 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP14_stroke')
         
         switch TrialDay
             case '01'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151102\'];
                 folderTASK_FLUO = [folderTASK,'sequence_trialx15_deep'];
                 folderTASK_ForceFileName = [folderTASK,'trialx15deep_sync.txt'];
             case '02'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151103\'];
                 folderTASK_FLUO = [folderTASK,'sequenza_trialx15_deep'];
                 folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
             case '03'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151104\'];
                 folderTASK_FLUO = [folderTASK,'sequenza_trialx15_deep'];
                 folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
             case '04'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151105\'];
                 folderTASK_FLUO = [folderTASK,'sequenza_trialx15_1mmdeep'];
                 folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
             case '05'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151106\'];
                 folderTASK_FLUO = [folderTASK,'sequenza_trialx15_1mmdeep'];
                 folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
         end

           
end
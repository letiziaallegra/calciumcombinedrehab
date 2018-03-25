%data animal load
function [folderTASK_FLUO folderTASK_ForceFileName] = fun_loadGCaMPanimal(UsbPort,Animal_Name,TrialDay)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   control   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% GCaMP3 %%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Animal_Name,'GCaMP3_control')  
    
    switch TrialDay
%         case '01'
%             folderTASK      = ['C:\Users\CNR-SSSUP\Desktop\LENS\150511_GCaMP3\'];
%             folderTASK_FLUO = [folderTASK,'Trial_3_gcamp3'];
% %             folderTASK_ForceFileName = [folderTASK,'1_sync.txt'];
%             folderTASK_ForceFileName = [folderTASK,'trial3_sync.txt'];
        case '02'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150512\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trial complete x15_TMP_sync.txt'];
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150513\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'TC15_F_sync.txt'];
        case '04'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150514\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150515\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trial15_ok_sync.txt'];
    end   
    
    
%%%%% GCaMP4 %%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP4_control')
    
     switch TrialDay
%          case '01'
%              folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150525\'];
%              folderTASK_FLUO = [folderTASK,'150525_GCaMP4_trialx15_8bits - ok dal trialx15_0274\tutte'];
%              folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
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
     
     
%%%%% GCaMPChR2_1_control %%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMPChR2_1_control')
    
    switch TrialDay
        case '01'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151123\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_1mmdeep_sync.txt'];
        case '02'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151124\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151125\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '04'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151126\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151127\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
    end
    
%%%%% GCaMPChR2_7_control %%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMPChR2_7_control')
    
    switch TrialDay
        case '01'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151207\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '02'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151208\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151209\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '04'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151210\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151211\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx14_deep_sync.txt'];
    end

    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   stroke    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% GCaMP9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP9_stroke')
    
    switch TrialDay
%         case '01'
%             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150727\'];
%             folderTASK_FLUO = [folderTASK,'150727_GCaMP9_trailx15_Surperficial_8bits\ok'];
%             folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '02'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150728\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15vero_sync.txt'];
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150729\'];
            folderTASK_FLUO = [folderTASK,'150729_GCaMP9_trialx15_8bits'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '04'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150730\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150731\'];
            folderTASK_FLUO = [folderTASK,'150731_GCaMP9_trailx15_8bits'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
    end
%%%% 


%%%% GCaMP10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP10_stroke')
    
     switch TrialDay
         case '01'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150803\'];
             folderTASK_FLUO = [folderTASK,'MAT_trial'];
             folderTASK_ForceFileName = [folderTASK,'2_primi4-5trials_poi si e sfilata la zampa_sync.txt'];
         case '02'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150804\'];
%              folderTASK_FLUO = [folderTASK,'MAT_trial2'];
%              folderTASK_ForceFileName = [folderTASK,'2_trialx6_deep_sync.txt'];
             folderTASK_FLUO = [folderTASK,'MAT_trial'];
             folderTASK_ForceFileName = [folderTASK,'1_trialx8+1_sync.txt'];
         case '03'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150805\'];
             folderTASK_FLUO = [folderTASK,'MAT_trial'];
             folderTASK_ForceFileName = [folderTASK,'1_trialx12_sync.txt'];
         case '04'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150806\'];
             folderTASK_FLUO = [folderTASK,'MAT_trial'];
             folderTASK_ForceFileName = [folderTASK,'2_trialx9_deep_sync.txt'];
         case '05'
             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150807\'];
             folderTASK_FLUO = [folderTASK,'MAT_trial'];
             folderTASK_ForceFileName = [folderTASK,'1_trialx8_sync.txt'];
     end
     
     
     
%%%% GCaMP11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP11_stroke')
         
         switch TrialDay
             case '01'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150808\'];
                 folderTASK_FLUO = [folderTASK,'MAT_trial'];
                 folderTASK_ForceFileName = [folderTASK,'1_trialx13_deep_sync.txt'];
             case '02'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150809\'];
                 folderTASK_FLUO = [folderTASK,'MAT_trial'];
                 folderTASK_ForceFileName = [folderTASK,'1_trialx15_1mmdeep_sync.txt'];
             case '03'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150810\'];
                 folderTASK_FLUO = [folderTASK,'MAT_trial'];
                 folderTASK_ForceFileName = [folderTASK,'2_trialx10_1mmdeep_sync.txt'];
             case '04'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150811\'];
                 folderTASK_FLUO = [folderTASK,'MAT_trial'];
                 folderTASK_ForceFileName = [folderTASK,'1_trialx10_1mmdeep_sync.txt'];
             case '05'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150812\'];
                 folderTASK_FLUO = [folderTASK,'MAT_trial'];
                 folderTASK_ForceFileName = [folderTASK,'2_trialx7_1mmdeep_sync.txt'];
         end
   
         
%%%% GCaMP14 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP14_stroke')
         
         switch TrialDay
             case '01'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151102\'];
                 folderTASK_FLUO = [folderTASK,'MAT_trial'];
                 folderTASK_ForceFileName = [folderTASK,'trialx15deep_sync.txt'];
             case '02'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151103\'];
                 folderTASK_FLUO = [folderTASK,'MAT_trial'];
                 folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
             case '03'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151104\'];
                 folderTASK_FLUO = [folderTASK,'MAT_trial'];
                 folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
             case '04'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151105\'];
                 folderTASK_FLUO = [folderTASK,'MAT_trial'];
                 folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
             case '05'
                 folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151106\'];
                 folderTASK_FLUO = [folderTASK,'MAT_trial'];
                 folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
         end
         
%%%% GCaMP15 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
elseif strcmp(Animal_Name,'GCaMP15_stroke')
    
    switch TrialDay
        case '01'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151102\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '02'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151103\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151104\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '04'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151105\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151106\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
    end
    
    
    
%%%% GCaMPChR2_8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMPChR2_8_stroke')
    
    switch TrialDay
        case '01'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_160111\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '02'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_160112\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx11_3_sync.txt'];
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_160113\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '04'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_160114\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx14_sync.txt'];
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_160115\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
    end
    
%%%% GCaMPChR2_9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMPChR2_9_stroke')
    
    switch TrialDay
        case '01'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_160111\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '02'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_160112\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_160113\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '04'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_160114\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_160115\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx12_sync.txt'];
    end
    
         
         
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   stroke + BonT/E  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% GCaMP12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP12_stroke_BoNT')
    
    switch TrialDay
        case '01'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150928\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '02'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150929\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx10_1mmdeep_sync.txt'];
        case '03'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_150930\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx10_1mmdeep_sync.txt'];
        case '04'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151001\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx10_1mmdeep_sync.txt'];
        case '05'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151002\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'tialx15_sync.txt'];
        case '06'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151005\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '07'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151006\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15deep_bis_sync.txt'];            
        case '08'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151007\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15deep_sync.txt'];
        case '09'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151008\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15deep_sync.txt'];
        case '10'
            %post staccato
            
        case '11'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151012\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15deep_sync.txt'];
        case '12'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151013\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_1mmdeep_sync.txt'];
        case '13'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151014\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '14'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151015\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_1mmdeep_sync.txt'];
        case '15'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151016\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '16'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151019\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_1mmdeep_sync.txt'];            
        case '17'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151020\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '18'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151021\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_1mmdeep_sync.txt'];
        case '19'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151022\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '20'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151023\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_1mmdeep_sync.txt'];           
    end
    
    
%%%% GCaMP16 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP16_stroke_BoNT')
    
    switch TrialDay
        case '15'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151204\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '16'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151207\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '17'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151208\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx13_sync.txt'];
        case '18'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151209\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15bis_sync.txt'];
        case '19'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151210\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx13_sync.txt'];
        case '20'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151211\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
    end
    
%%%% GCaMP17 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP17_stroke_BoNT')
    
    switch TrialDay
%         case '16'
%             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151207\'];
%             folderTASK_FLUO = [folderTASK,'MAT_trial'];
%             folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
%         case '17'
%             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151208\'];
%             folderTASK_FLUO = [folderTASK,'MAT_trial'];
%             folderTASK_ForceFileName = [folderTASK,'trialx10_sync.txt'];
%         case '18'
%             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151209\'];
%             folderTASK_FLUO = [folderTASK,'MAT_trial'];
%             folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
%         case '19'
%             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151210\'];
%             folderTASK_FLUO = [folderTASK,'MAT_trial'];
%             folderTASK_ForceFileName = [folderTASK,'trialx8_sync.txt'];
%         case '20'
%             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151211\'];
%             folderTASK_FLUO = [folderTASK,'MAT_trial'];
%             folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
    end
    
%%%% GCaMP18 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMP18_stroke_BoNT')
    
    switch TrialDay
%         case '16'
%             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151207\'];
%             folderTASK_FLUO = [folderTASK,'MAT_trial'];
%             folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
%         case '17'
%             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151208\'];
%             folderTASK_FLUO = [folderTASK,'MAT_trial'];
%             folderTASK_ForceFileName = [folderTASK,'trialx8_sync.txt'];
        case '18'
            folderTASK      = [UsbPort,':\LENS\Animals Data _ da sistemare\',Animal_Name,'\',TrialDay,'_151209\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
%         case '19'
%             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151210\'];
%             folderTASK_FLUO = [folderTASK,'MAT_trial'];
%             folderTASK_ForceFileName = [folderTASK,'trialx14_sync.txt'];
%         case '20'
%             folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151211\'];
%             folderTASK_FLUO = [folderTASK,'MAT_trial'];
%             folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
    end
    
%%%% GCaMPChR2_3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Animal_Name,'GCaMPChR2_3_stroke_BoNT')
    
    switch TrialDay
        case '15'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151204\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx12_deep_sync.txt'];
        case '16'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151207\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '17'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151208\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_deep_sync.txt'];
        case '18'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151209\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '19'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151210\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx15_sync.txt'];
        case '20'
            folderTASK      = [UsbPort,':\LENS\Animals Data\',Animal_Name,'\',TrialDay,'_151211\'];
            folderTASK_FLUO = [folderTASK,'MAT_trial'];
            folderTASK_ForceFileName = [folderTASK,'trialx12_deep1_sync.txt'];
    end
    
     
    
    
end %end function


























%Script To compute Parameters from Force Peaks
%ComputeForceParScript

% clear
% close all
clc

%%%%  Folder Info
UsbPort = 'I';
User = getenv('username');
%%%

%%%%  Choice of the animal and trial day
%%%%%%%%%%%
% ListAnimalTogether = {'GCaMP3_control', 'GCaMP4_control', ...
%                       'GCaMP9_stroke','GCaMP10_stroke', 'GCaMP11_stroke', 'GCaMP14_stroke', 'GCaMP15_stroke',...
%                       'GCaMPChR2_1_control'...
%                       'GCaMPChR2_7_control', 'GCaMPChR2_17_control', 'GCaMPChR2_18_control',...
%                         'GCaMPChR2_20_control', 'GCaMPChR2_21_control', 'GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                         'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke',...
%                         'GCaMPChR2_3_stroke_BoNT',...
%                         'GCaMP16_stroke_BoNT','GCaMP17_stroke_BoNT','GCaMP18_stroke_BoNT','GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                         'GCaMPChR2_13_stroke_BoNT',...
%                         'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};
                    

%%%%
AnimalMainDir        = [UsbPort,':\LENS\Animals Data'];
ListPos = [];
CurrDir = pwd;
for LATo = 1:length(ListAnimalTogether) %list animals
    
    Animal_Name                    = [ListAnimalTogether{LATo}];
    
    AnimalCurrDir                  = [AnimalMainDir,'\',Animal_Name];
    
    %%%%%%%%
    
    
    %%%%%%%%
    ListFolderAnimalDir = dir(AnimalCurrDir);
    
    
    %for days
    for lfcd=3:length(ListFolderAnimalDir)
        
        DayCurrDir           = ListFolderAnimalDir(lfcd,1).name;
        display([Animal_Name,'_',DayCurrDir])
        ListFolderDayCurrDir = dir([AnimalCurrDir,'\',DayCurrDir]);
        
        
        %load(Pos)
        fileTxT_OK = 0;
        if exist([AnimalCurrDir,'\',DayCurrDir,'\','Pos.mat'])
            
            %load POS
            load([AnimalCurrDir,'\',DayCurrDir,'\','Pos.mat'])
            ListPos{LATo,lfcd} = [Animal_Name,'\',DayCurrDir,'_ok'];
            
            

            %files inside days
            for fl=3:length(ListFolderDayCurrDir)
                
                NameFile = ListFolderDayCurrDir(fl,1).name;
                 
                %dataGCamp
                if strfind(NameFile,'dataMouseGCamp')
                    
                    if strcmp(NameFile,['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),'_Par_MIP.mat'])
                        
                        load([AnimalCurrDir,'\',DayCurrDir,'\',NameFile])
                        M_Par_MIP = dataGCamp;
                        clear dataGCamp
                        
                    elseif strcmp(NameFile,['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),'_Par.mat'])
                        
                        load([AnimalCurrDir,'\',DayCurrDir,'\',NameFile])
                        M_Par = dataGCamp;
                        clear dataGCamp
                        
                        
                    elseif strcmp(NameFile, ['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),'.mat'])
                        
                        load([AnimalCurrDir,'\',DayCurrDir,'\',NameFile])
                        M = dataGCamp;
                        clear dataGCamp
                    end
                    
                end
                
                %file TXT _sync
                if strfind(NameFile,'_sync')
                    
                    dataTXT = importdata([AnimalCurrDir,'\',DayCurrDir,'\',NameFile]);
                    dataTXT_name = NameFile;
                    dataTXT = dataTXT(:,1:8);
                    fileTxT_OK = 1;
                
                end
                
            end
                        
            if fileTxT_OK
                
                cd(CurrDir)
                dataTXT_updated = synchronize_vWebCam_LENS_v3_for_Update_Position_Sync(res, dataTXT,8);
                clear res dataTXT
                
                cd([AnimalCurrDir,'\',DayCurrDir])
                
                save(dataTXT_name, 'dataTXT_updated', '-ascii');
                
                %allineamento con data in dataGCamp
                zeroStart       = find(dataTXT_updated(:,7)==0);
                dataTXT_updated(1:zeroStart(1)-1,:) = [];
                difflen         = size(dataTXT_updated,1) - size(M_Par_MIP.fx,1);
                dataTXT_updated = dataTXT_updated(1:end-difflen,:);
                
                %_Par_MIP
                M_Par_MIP.pos   = dataTXT_updated(:,9);
                M_Par_MIP.speed = dataTXT_updated(:,10);
                M_Par_MIP.acc   = dataTXT_updated(:,11);                
                dataGCamp = M_Par_MIP;
                save(['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),'_Par_MIP.mat'],'dataGCamp')
                if 1
                    figure('Name','Fig_Sig_MIP_Area')
                    hold on
                    plot(dataGCamp.t,dataGCamp.status/10,'r');
                    plot(dataGCamp.t,dataGCamp.fx,'b');
                    plot(dataGCamp.t,dataGCamp.pos/10,'g');
                    plot(dataGCamp.t,dataGCamp.speed/100,'m');
                    plot(dataGCamp.t,dataGCamp.ROI_MIP.ROI_Signal(:,1)/10,'k');
                    plot(dataGCamp.t,dataGCamp.ROI_MIP.ROI_Signal(:,2)/10,'g');
                    plot(dataGCamp.t,dataGCamp.ROI_MIP.ROI_Signal(:,3)/10,'c');
                    legend({'Status','Force','Pos/10','Speed/100','max MIP/10','sum A MIP/10','sum P MIP/10'})
                    saveas(gca,['Fig_Sig_MIP_',Animal_Name,'_',DayCurrDir(1:2)]);
                    close
                    
                    figure('Name','Fig_Sig_Whole_Area')
                    hold on
                    plot(dataGCamp.t,dataGCamp.status/10,'r');
                    plot(dataGCamp.t,dataGCamp.fx,'b');
                    plot(dataGCamp.t,dataGCamp.pos/10,'g');
                    plot(dataGCamp.t,dataGCamp.speed/100,'m');
                    plot(dataGCamp.t,dataGCamp.fluoROI/10,'g');
                    legend({'Status','Force','Pos/10','Speed/100','Whole Area/10'})
                    saveas(gca,['Fig_whole_sig_',Animal_Name,'_',DayCurrDir(1:2)]);
                    close
                end
                clear M_Par_MIP
                
                %_Par
                M_Par.pos   = dataTXT_updated(:,9);
                M_Par.speed = dataTXT_updated(:,10);
                M_Par.acc   = dataTXT_updated(:,11);                
                dataGCamp = M_Par;
                save(['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),'_Par.mat'],'dataGCamp')
                clear M_Par
                
                %_Par
                M.pos   = dataTXT_updated(:,9);
                M.speed = dataTXT_updated(:,10);
                M.acc   = dataTXT_updated(:,11);                
                dataGCamp = M;
                save(['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),'.mat'],'dataGCamp')
                clear M               
                
                clear dataTXT_updated
                
            else                
                ListPos{LATo,lfcd} = [Animal_Name,'\',DayCurrDir,'_okPos_NO_SYNC'];                
            end
                           
            
        else
            display([AnimalCurrDir,'\',DayCurrDir,'\','Pos.mat'])
            ListPos{LATo,lfcd} = [Animal_Name,'\',DayCurrDir,'_NO_POS'];
        end
                
    end
  
   
end %end LATo

cd(CurrDir)
save('ListPos','ListPos')












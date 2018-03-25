%carica il ImDelay_Day_M_All_NOMEANIMALE_control_MEDIA_TOT
%carica il DataIm_Mean_Delay_All_Days_NOMEANIMALE
%e poi lancia questo script

% %
UserName = 'CNR-SSSUP';
UsbPort = 'I';
%
% UserName = 'Stefano';
% UsbPort = 'F';

%MEAN or MEDIAN
% MeaMed = 1; %MEAN
MeaMed = 2; %MEDIAN
 

%%%%%%%% load data File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if MeaMed == 1
    %%%%%%% loading data folder
    DirData      = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\_Analisi_Corr_Iso_MIP\Isochrones\Data_Isochrones_SelSeq_MAX\Data_Delay_Mean of Trials'];       
    %%%%%%% save data folder
    SaveDir      = [UsbPort,':\LENS\Isocrone_SelSeq_MAX\Figures_Weekly_Mean of Trials'];
    
elseif MeaMed == 2
    %%%%%%% loading data folder
    DirData      = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\_Analisi_Corr_Iso_MIP\Isochrones\Data_Isochrones_SelSeq_MAX\Data_Delay_Median of Trials'];       
    %%%%%%% save data folder
    SaveDir      = [UsbPort,':\LENS\Isocrone_SelSeq_MAX\Figures_Weekly_Median of Trials'];
    
end
%%%%%%% directory to find right list of consecutive days (ci possono essere dei days mancanti in alcuni animali)
DirCons      = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\Script_Flip_Find_References\MAT_Rot_Trans'];


for n_i = 1:3
    
    switch n_i
        
        case 1
            %control
            TreatDir = '\Control';
            DirDataTreat = [DirData,TreatDir];
            num_week         = 1;
            num_day_per_week = 5;
            list_days = [1:5];
        case 2
            %stroke
            TreatDir = '\Stroke';
            DirDataTreat = [DirData,TreatDir];
            num_week         = 1;
            num_day_per_week = 5;
            list_days = [1:5];
        case 3
            %rehab
            TreatDir = '\Rehab';
            DirDataTreat = [DirData,TreatDir];
            num_week         = 4;
            num_day_per_week = 5;
            %list real day
            list_days = [1:5;6:10;11:15;16:20];
    end
    
    NameFileDir = dir(DirDataTreat);
    
    for nfd_i = 3:length(NameFileDir)
        
        FileData = NameFileDir(nfd_i,1).name;
        
        if isempty(strfind(FileData,'.ini')) % problem google drive
            
            load([DirDataTreat,'\',FileData])
            
            %%%%check name
            if exist('ImDelay_Day_M_All','var') %if exist
                
                i_name = strfind(FileData,'GCaMP');
                AnimalName = FileData(i_name:end-4);
                
                for i_w=1:num_week % if week
                    
                    weekly_index = [];
                    
                    load([DirCons,'\',AnimalName,'_Rot_Trans_Par'])
                    list_real_days = rot_transl(:,1);
                    
                    for i_rdy = 1:num_day_per_week
                        weekly_index = [weekly_index; find(list_real_days == list_days(i_w,i_rdy))]; 
                    end

                    %exception
                    if strcmp(AnimalName,'GCaMPChR2_25_stroke')
                        weekly_index = 1:4;
                    end
                    
                    if MeaMed == 1
                        IsoMean     = nanmean(ImDelay_Day_M_All(:,:,weekly_index),3);
                    elseif MeaMed == 2
                        IsoMean     = nanmedian(ImDelay_Day_M_All(:,:,weekly_index),3);
                    end                 
                    
                    IsoMean_Big = imresize(IsoMean, 40);
                    
                    %%name
                    fig_Name = [AnimalName,'_Week_',num2str(i_w),'_5_days_average'];
                    
                    %%plot
                    fig_Mean = figure('Name',fig_Name);
                    imagesc(IsoMean_Big)
                    changeLabel_MIP(size(IsoMean_Big,1),size(IsoMean_Big,2));
%                     caxis([-0.1 0.3])
                    colorbar
                    
                    %save Images Isochrones
                    saveas(gca,[SaveDir,TreatDir,'\',fig_Name],'fig')
                    saveas(gca,[SaveDir,TreatDir,'\',fig_Name],'jpeg')
                    close
                    
                    if n_i==3
                        
                        if i_w == 1
                            fig4 = figure('Name',[AnimalName,'_Weeks_5_days_average']);
                        end
                        
                        figure(fig4)
                        subplot(2,2,i_w)
                        imagesc(IsoMean_Big)
                        changeLabel_MIP(size(IsoMean_Big,1),size(IsoMean_Big,2));
                        caxis([-0.1 0.6])
                        colorbar
                        
                        if i_w == 4
                            saveas(gca,[SaveDir,'\',AnimalName,'_Weeks_5_days_average'],'fig')
                            saveas(gca,[SaveDir,'\',AnimalName,'_Weeks_5_days_average'],'jpeg')
                            close
                        end
                        
                    end
                    
                    
                end % end if week
                
            end %end if exist
            %%%%
            
        end %end problem google drive
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


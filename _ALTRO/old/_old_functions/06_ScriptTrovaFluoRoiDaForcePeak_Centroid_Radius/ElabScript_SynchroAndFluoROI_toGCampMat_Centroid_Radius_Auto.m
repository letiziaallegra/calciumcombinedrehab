%
% script to syncrhonize fluo and force signals and
% define the ROIs signals based on the Centroid and the Radius assessed
%

close all
clear all
clc
CurrDir = cd;


User = 'CNR-SSSUP';
UsbPortHD  = 'I';
% Animal = [];   %all the animals
Animal = 'GCaMP16_stroke_BoNT';   

%%%%%%%%%%% dataGCamp (ROI arbitrarie) -> folder where saving data
MainDir = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_data_MAT_GCamp\'];

%%%%%%%%%%% centroid
MainDirCentroid = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_data_MAT_GCamp_Centroid\'];
NumAnimalCentroidFolder = dir(MainDirCentroid);

if ~isempty(Animal)     
    NumAnimalCentroidFolder = [1 2 3];
end


%%%%%% START %%%%%%%%%%%

for nAn_i=3:length(NumAnimalCentroidFolder)
    
    %%%%%%%%%%%%
    %name animal
    if isempty(Animal)  
        Animal_name = NumAnimalCentroidFolder(nAn_i,1).name;
    else
        Animal_name = Animal;
    end
    
    %%%%%%%%%%%%
    %main centroid matrix
    DirCentr = [MainDirCentroid,Animal_name,'\Centroid_Matrix'];
    %load Centroid Matrix (centroid for each day) -> X_Y_Weighted_Centroid_DAY
    load([DirCentr,'\','X_Y_Weighted_Centroid_DAY'])
    %%%%%%%%%%%%%
    
    %list days to take
    ListDays = X_Y_Weighted_Centroid_DAY(:,1);
    
    for nDay_i=1:length(ListDays)
        
        currDay = ListDays(nDay_i);
        if currDay<=9
            currDay_str = ['0',num2str(currDay)];
        else
            currDay_str = [num2str(currDay)];
        end
        
        %load single dataGCamp matrix
        DirGCampMatrix      = [MainDir,'\',Animal_name,'\',currDay_str ];
        FileNameGCampMatrix = ['dataMouseGCamp_',Animal_name,'_',currDay_str,'_Par'];
        
        %load matrix GCamp
        load([DirGCampMatrix,'\',FileNameGCampMatrix ]);
        
        
        %%% Initial Check %%%
        if exist('dataGCamp')
            if isfield(dataGCamp,'PeaksPar_Fx_Fluo')
                display('Check OK')
            else
                error('PeaksPar_Fx_Fluo missing: perform the parameters extraction')
            end
        else
            error('load dataGCamp_Par')
        end
        %%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Radius assessed
        Radius        = round(167.51/2); %[pixels] <- FISSO
        %Centroids
        Centroid_List = round(X_Y_Weighted_Centroid_DAY(nDay_i,2:3)); %[pixels] 
        ROI_num_Centroid = size(Centroid_List,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%% Info images store %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        UsbPort = UsbPortHD;
        %
        Animal_Name = dataGCamp.Info.Name;
        TrialDay    = dataGCamp.Info.Date;
        %interval of frames to compute f0
        IntervalToFluoMean = dataGCamp.Info.IntervalToFluoForceMean_PointFluoFreq;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%% data folder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cd ..
        cd('_load_animal_fun')
        [folderTASK_FLUO ~] = fun_loadGCaMPanimal(UsbPort,Animal_Name,TrialDay);
        cd(CurrDir)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%% Time and Frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t      = dataGCamp.t;
        Fs     = dataGCamp.Info.Fs;
        FsFluo = dataGCamp.Info.Fs_FluoImages; %[Hz]
        ResampleParRobotFluo = Fs/FsFluo;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        %%% data Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = dir(folderTASK_FLUO);
        
        len = length(d);
        NumImages = len-2;
        
        % - starting point
        RealStart_Fluo     = 1;
        
        %resampling fsample Robot Fluo
        ResampleParRobotFluo = round(Fs/FsFluo);
        %duration (it should be already ok because of the initial check done in 02_ScriptTrovaFluoRoiDaForcePeak)
        RealDurTask_Fluo = floor(length(t)/ResampleParRobotFluo);
        
        % - ending point
        LenFLUO_Im = RealStart_Fluo+RealDurTask_Fluo-1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        % %storage of each f0 for each ROI
        % MatrixImage_f0_ROI = zeros(sz_Im(1),sz_Im(2),ROI_num_Centroid);
        
        index=0;
        downsamplingfactor=1;
        
        wb = waitbar(0,'Images Loading, Please wait...');
        
        %%% load data Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LenFLUO_Im = RealStart_Fluo+RealDurTask_Fluo-1;
        ROI_Centroid_MeanSignal = [];
        
        for i=RealStart_Fluo :downsamplingfactor:   LenFLUO_Im
            
            index = index+1;
            indexImage = i+2;
            waitbar(indexImage/(LenFLUO_Im),wb)
            
            nameImage = d(indexImage,1).name;
            
            Im = imread([folderTASK_FLUO,'\',nameImage]);
            Im_Original = Im;
            sz_Im = size(Im_Original);
            
            %filtering
            Im_Original = medfilt2(Im_Original);
            
            for nR = 1:ROI_num_Centroid
                
                if index==1
                    
                    %%%%%%%%%%%% region of interest based on Centroid and Radius %%%%%%%%%%%%
                    fROI = figure('Name','Show the ROIs','NumberTitle','off');
                    hIm = imshow(Im_Original);
                    close(fROI)                    
                    
                    % Make a logical image with the selected circular region set to 1, the rest % to zero
                    [xgrid, ygrid] = meshgrid(1:sz_Im (2), 1:sz_Im (1));
                    
                    Cx = xgrid - Centroid_List(nR,1);    % offset the origin
                    Cy = ygrid - Centroid_List(nR,2);
                    circlemask = Cx.^2 + Cy.^2 <= Radius.^2;
                    
                    
                    %%% Trova f0 di questa area circolare %%%
                    MatrixImage_f0 = zeros(sz_Im(1),sz_Im(2),IntervalToFluoMean(2)-IntervalToFluoMean(1)+1);
                    indexStore_f0 = 0;
                    for fi=IntervalToFluoMean(1):IntervalToFluoMean(2)
                        
                        %index to report to images list
                        indexImage_f0 = fi+2;
                        %name
                        nameImage_f0 = d(indexImage_f0,1).name;
                        
                        %load image
                        Im_f0 = imread([folderTASK_FLUO,'\',nameImage_f0]);
                        Im_Original_f0 = Im_f0;
                        
                        %filtering
                        Im_Original_f0 = medfilt2(Im_Original_f0);
                        
                        % Use the mask to select part of the image
                        circle_image_f0 = double(Im_Original_f0) .* circlemask;
                        circle_image_f0(circle_image_f0==0) = NaN;
                        
                        %
                        indexStore_f0 = indexStore_f0+1;
                        MatrixImage_f0(:,:,indexStore_f0) = circle_image_f0;
                        
                    end
                    display('F0 extracted')
                    
                    
                    MeanMatrix_f0 = nanmean( double(MatrixImage_f0),3); %in gray tones
                    MeanMatrix_f0(MeanMatrix_f0==0) = 1;
                    
                    H_Im_Round_Fig = figure('Name','circle_ROI_around_Main_Centroid');
                    imshow(MeanMatrix_f0, []);
                    H_Im_Round_Fig_filename = ['Fig_ROI_around_Main_Centroid_',Animal_name,'_',currDay_str];
                    saveas(H_Im_Round_Fig ,[DirGCampMatrix,'\',H_Im_Round_Fig_filename],'fig');
                    close(H_Im_Round_Fig)
                    
                    %storage delle f0 delle ROI
                    MatrixImage_f0_ROI(:,:,nR) = MeanMatrix_f0;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                end
                
                % Use the mask to select part of the image
                circle_image = double(Im_Original) .* circlemask;
                %Metto tutti gli 0 pari a NaN
                circle_image(circle_image==0) = NaN;
                
                %fai delta_f/f0
                Im = (circle_image-squeeze(MatrixImage_f0_ROI(:,:,nR)))./squeeze(MatrixImage_f0_ROI(:,:,nR))*100;
                
                %compute average of the image (just the circle bc the outside the circle is all NaN
                ROI_Centroid_MeanSignal(i,nR)     = squeeze(nanmean(nanmean(Im,2),1));
                
            end
        end
        
        imFluoResMean = [];
        for nR = 1:ROI_num_Centroid
            
            imFluoResMean_Buf = resample(ROI_Centroid_MeanSignal(:,nR),ResampleParRobotFluo,1);
            imFluoResMean_Buf = imFluoResMean_Buf - median(imFluoResMean_Buf);
            
            imFluoResMean     = [imFluoResMean imFluoResMean_Buf];
        end
        
        
        %%%%%%%%%%%%%%%
        %% PLOT DATA %%
        %%%%%%%%%%%%%%%
        H_MCA_fig = figure('Name','Fig_Sig_Main_Centroid_Area');
        plot(t,dataGCamp.status/10,'r')
        hold on
        plot(t,dataGCamp.fx,'b')
        plot(t,dataGCamp.speed/100,'m')
        plot(t,imFluoResMean/100,'Color',[0 0.8 0])
        legend([{'Status','Force','Speed','Main Fluo ROI'}])
        %saving figure
        H_MCA_fig_filename = ['Fig_Sig_Main_Centroid_Area_',Animal_name,'_',currDay_str];
        saveas(H_MCA_fig,[DirGCampMatrix,'\',H_MCA_fig_filename],'fig');
        close(H_MCA_fig)
        %%%%
        
        
        %%%%%%%%%%%%%%%
        %% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%
        
        % image f0 ROI
        dataGCamp.MainCentroid.ImageROI_f0   = MatrixImage_f0;
        % centroid List
        dataGCamp.MainCentroid.Centroid_List = Centroid_List;
        % Radius
        dataGCamp.MainCentroid.Radius        = Radius;
        % signals
        dataGCamp.MainCentroid.ROI_Signal    = imFluoResMean;
        
        
        %%%%%%
        filename_GCamp = ['dataMouseGCamp_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,'_Par','_MainCentroid' ];
        save([DirGCampMatrix,'\',filename_GCamp],'dataGCamp')
        %%%%%%
        
        display(['END PROCESS for: ',filename_GCamp]);
        clear dataGCamp MatrixImage_f0 MatrixImage_f0_ROI MeanMatrix_f0 Im Cx Cy circle_image circlemask circle_image_f0 d imFluoResMean imFluoResMean_Buf t xgrid ygrid Im_Original Im_Original_f0 Im_f0
    end
end






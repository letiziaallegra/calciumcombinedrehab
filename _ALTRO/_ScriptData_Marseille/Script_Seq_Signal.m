% Script to make a matrix of the data

clear all
close all
clc

%%%%%%%%%%%
CurrDir = cd;

%%Animal list
ListAnimalTogether = {'GCaMPChR2_26_stroke',...
                      'GCaMP16_stroke_BoNT','GCaMP18_stroke_BoNT',...
                      'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
                      'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT'};

%%Folders
UserName = 'CNR-SSSUP';
UsbPort = 'I';
% UserName = 'Stefano';
% UsbPort = 'F';

%Dim new matrix
Resol    = 512;
Resol_MM = 128; %number of region for each region

DATAFOLDER = [UsbPort,':\LENS\Animals Data'];
SAVEFOLDER = ['C:\Users\CNR-SSSUP\Desktop\LENS_Dati_Marseille'];


for lat=1:length(ListAnimalTogether) %for ANIMAL
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Animal_name = ListAnimalTogether{lat};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MainDir = [DATAFOLDER,'\',Animal_name,'\'];
    NumDaysFolder = dir(MainDir);
    
    %%%%%%%% load Reference File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RefDir      = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\Script_Flip_Find_References\MAT_Rot_Trans'];
    FileRefName = [Animal_name,'_Rot_Trans_Par.mat'];
    if exist([RefDir,'\',FileRefName])
        % rot_transl
        load([RefDir,'\',FileRefName]);
    else
        error([RefDir,'\',FileRefName,' is not present in the folder']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MIP_Th_Dil_filt_All_Day = [];
    
    degree_all = [];
    trans_all  = [];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nd_i=3:length(NumDaysFolder) %for days
        
        index_day = nd_i-2;
        nday_date = NumDaysFolder(nd_i,1).name; 
        nday      = nday_date(1:2);
                
        CurrAnDayFolder = [MainDir,nday_date];
        CurrAnDayFolder_List = dir(CurrAnDayFolder);
              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%find image folder and force file        
        for cadf=3:length(CurrAnDayFolder_List)
            
            if strcmp(CurrAnDayFolder_List(cadf,1).name(1:3),'MAT')
                %data images folder
                folderTASK_FLUO = [CurrAnDayFolder,'\',CurrAnDayFolder_List(cadf,1).name];
            elseif length(CurrAnDayFolder_List(cadf,1).name)>7
                if strcmp(CurrAnDayFolder_List(cadf,1).name(end-7:end-4),'_Par')
                    %dataGCmap filename
                    folderTASK_dataGCampFilename = [CurrAnDayFolder,'\',CurrAnDayFolder_List(cadf,1).name];
                end
            end
            
        end
    
        %%load GCamp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(folderTASK_dataGCampFilename);
        
        %%% data Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = dir(folderTASK_FLUO);
        len = length(d);
        NumImages = len-2;
        
        %default
        rw = Resol;
        cl = Resol;
        
        %%baseline for deltaf/f
        if 1
            %interval of frames to compute f0
            IntervalToFluoMean = dataGCamp.Info.IntervalToFluoForceMean_PointFluoFreq;
            
            MatrixImage_f0 = zeros(rw,cl,IntervalToFluoMean(2)-IntervalToFluoMean(1)+1);
            indexStore_f0  = 0;
            
            for fi=IntervalToFluoMean(1):IntervalToFluoMean(2)
                %index to report to images list
                indexImage_f0 = fi+2;
                %name
                nameImage_f0 = d(indexImage_f0,1).name;
                
                %load image
                load([folderTASK_FLUO,'\',nameImage_f0]);
                Im_f0 = Im8_fv;
                Im_Original_f0 = Im_f0;
                
                %filtering
                Im_Original_f0 = medfilt2(Im_Original_f0);
                indexStore_f0 = indexStore_f0+1;
                MatrixImage_f0(:,:,indexStore_f0) = Im_Original_f0;
            end
            
            %baseline Matrix
            MeanMatrix_f0 = nanmean( double(MatrixImage_f0),3); %in gray tones
            MeanMatrix_f0(MeanMatrix_f0==0) = 1;
        end
        %%%
        
        
%         wb = waitbar(0,['Images Processing, Please wait...']);      
        AnimalMatrix_Day = zeros(Resol_MM,Resol_MM,NumImages);
        
        %%scroll images
        for indexImage=1:NumImages % for image stack
            
%             waitbar(indexImage/(NumImages),wb)
            display( num2str( indexImage/NumImages*100 ) )

            nameImage = d(indexImage+2,1).name;
            %load .mat
            load([folderTASK_FLUO,'\',nameImage]);
            Im8_fv = medfilt2(Im8_fv);
            
            Im_Original = double(Im8_fv);
            
            rw        = size(Im_Original,1);
            cl        = size(Im_Original,2);  
            
            %delta_f/f0
            Im_Original_Norm = (Im_Original-MeanMatrix_f0) ./ MeanMatrix_f0 *100;
            
            %%%% registration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            i_day_actual = find(rot_transl(:,1) ==   str2num(nday) );
            
            %%%%%%%%% ROTATION 
            degree = rot_transl(i_day_actual,2);
            if  degree~= 0
                Im_R = imrotate(Im_Original_Norm,degree,'crop');
            else
                Im_R = Im_Original_Norm;
            end
            
            
            Im_OR          = zeros(rw,cl);
            Im_OR_RotTrasl = zeros(rw,cl);
            
            %%%%%%%%% TRANSLATION 
            %left/right transl
            trans_lr  = rot_transl(i_day_actual,3);
            %translation Y along rows
            trans_Y   = rot_transl(i_day_actual,4);
            %translation X along columns
            trans_X   = rot_transl(i_day_actual,5);
            
            %translation along rows
            if trans_lr<0
                Im_OR(1:rw-trans_Y+1,:) = Im_R(trans_Y:end,:);
                trans_lr = -1;
            elseif trans_lr>0
                Im_OR(trans_Y:end,:) = Im_R(1:rw-trans_Y+1,:);
                trans_lr = +1;
            end
            
            %translation along columns
            Im_OR_RotTrasl(:,1:cl-trans_X+1) = Im_OR(:,trans_X:end);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            %%final 512x 512 registered and normalized image
            Frame_ith = Im_OR_RotTrasl;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %downsampling of the data matrix
            
            %%initialization of the image with a total number of region = N_reg.
            %%The third dimension of ROI_VALUE are the frames related to
            %%that day-stack
            MM = zeros(Resol_MM,Resol_MM);
            
            for n_reg_col=1:round(Resol_MM) %for # num regions-col in the Frame_ith
                
                index_col = n_reg_col-1;
                
                C_st = round(Resol/Resol_MM)*index_col+1;
                C_en = C_st+round(Resol/Resol_MM)-1;                
                
                for n_reg_row=1:round(Resol_MM) %for # num regions-row in the Frame_ith
                    
                    index_row = n_reg_row-1;
                    
                    R_st = round(Resol/Resol_MM)*index_row+1;
                    R_en = R_st+round(Resol/Resol_MM)-1;
                    
                    ROI_SQUARED = Frame_ith(  R_st:R_en, C_st:C_en);
                    
                    %%%
                    MM(n_reg_row, n_reg_col) = nanmedian(nanmedian(ROI_SQUARED));
                    
                end
            end
            
            %%store in the big matrix
            AnimalMatrix_Day(:,:,indexImage) = MM;            
            
            
        end % end for image stack
        
        filename = [Animal_name,'_',nday_date];
        save([SAVEFOLDER,'\',filename], 'AnimalMatrix_Day')
        
        display(['End Process for ', filename])
        clear AnimalMatrix_Day Im
        
    end
    display(['End Process for ', Animal_name])
    
end     
    
display('End Whole Process')
    
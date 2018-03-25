% Script to make Isochrones from Sequence_LONG
clear all
close all
clc

%%%%%%%%%%%
CurrDir = cd;
ListAnimalTogether = {  'GCaMPChR2_7_control'};
% ListAnimalTogether = {  'GMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};

% ListAnimalTogether = {      'GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                             'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke',...
%                             'GCaMP16_stroke_BoNT','GCaMP17_stroke_BoNT','GCaMP18_stroke_BoNT',...
%                             'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                             'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};

for lat=1:length(ListAnimalTogether)
%     
%     UserName = 'CNR-SSSUP';
%     UsbPort = 'I';
    UserName = 'Stefano';
    UsbPort = 'H';
    
    %%%%REGION OF INTEREST
    Resol     = 512;
    Resol_SQ  = Resol^2;
%     N_reg     = 1024;                 %total number of region
    N_reg     =128;                 %total number of region
    N_reg_Per_Border = sqrt(N_reg); %number of region for each region
    Lat       = Resol/N_reg_Per_Border; %Size of the border of each region
%     Lat       = sqrt(Resol_SQ/N_reg); %Size of the border of each region
    %%%%
    
    %%% Frequency
    Fs = 25;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Animal_name = ListAnimalTogether{lat};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%
%     MIP_tech_List = {'max','sum'};
    MIP_tech_List = {'sum'};
    %%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MainDir = [UsbPort,':\LENS\Animals Data\',Animal_name,'\'];
    NumDaysFolder = dir(MainDir);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%% load Reference File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RefDir      = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\Script_Flip_Find_References\MAT_Rot_Trans'];
    FileRefName = [Animal_name,'_Rot_Trans_Par.mat'];
    if exist([RefDir,'\',FileRefName])
        % rot_transl
        load([RefDir,'\',FileRefName]);
    else
        error([RefDir,'\',FileRefName,' is not present in the folder']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MIP_Th_Dil_filt_All_Day = [];
    
    degree_all = [];
    trans_all  = [];
    
    for iMT=1:length(MIP_tech_List)
        
        MIP_tech = MIP_tech_List{iMT};
        
        
        for nd_i=3:length(NumDaysFolder)
            
            index_day = nd_i-2;
            nday_date = NumDaysFolder(nd_i,1).name;
            nday      = nday_date(1:2);
            
            
            DirDay      = [MainDir,nday_date,'\','SequenceLong'];            
            SeqDir      = dir(DirDay);
            
            if  ~isempty(SeqDir)
                
                %seq ok
                SeqDirOKNO = 1;
                
                %load Image Sequence
                load([DirDay,'\',SeqDir(3,1).name])
                
                %%% Initial Check %%%
                if ~exist('ImageSequence')
                    error('load ImageSequence')
                end
                %%%
                
%                 %%% Dir where saving
%                 DirDayToSaveMIP = [MainDir,nday_date,'\',''];
%                 if ~isdir(DirDayToSaveMIP)
%                     %make Folder where saving
%                     mkdir(DirDayToSaveMIP)
%                 end
%                 %%%
                
                %%%%%%%%% TASK  %%%%%%%%%%%%%
                MatrixImageSequence = ImageSequence.MatrixImageSequence;
                NumTrials = size(MatrixImageSequence,1);
                TrialsUsed    = [1:NumTrials];
                % TrialsUsed    = [1];
                NumTrialsUsed = length(TrialsUsed);
                rw        = size(MatrixImageSequence{1,1},1);
                cl        = size(MatrixImageSequence{1,1},2);
                NumFrames = size(MatrixImageSequence{1,1},3);
%                 
%                 %to store median/mean
%                 Med_Frames_Matrix = zeros(rw,cl,NumFrames);
%                 %to store median/mean -> threshold
%                 Med_Frames_Matrix_Th = zeros(rw,cl,NumFrames);
                %Fixed Threshold
                % Threshold = 0.05;
                
                %                 MIP_all = zeros(rw, cl,  NumFrames*NumTrialsUsed);
                
                ind_i = 0;
                ImDelay = [];
                for int=1:NumTrialsUsed %for trials
                    
                    
                    ROI_VALUE = zeros(round(Resol/Lat),round(Resol/Lat),NumFrames);
                    
                    
                    for inf=1:NumFrames %for frames in the trial i-th
                        
                        Frame_ith = MatrixImageSequence{int,:}(:,:,inf);
                        
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        %%%%%% same rotation and translation %%%%%%%
                        i_day_actual = find(rot_transl(:,1) == str2num(nday));
                        
                        %%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        degree = rot_transl(i_day_actual,2);
                        if  degree~= 0
                            Frame_R = imrotate(Frame_ith,degree,'crop');
                        else
                            Frame_R = Frame_ith;
                        end
                        
                        Frame_OR = zeros(rw,cl);
                        Frame_OR_2 = zeros(rw,cl);
                        
                        %%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %left/right transl
                        trans_lr  = rot_transl(i_day_actual,3);
                        %translation Y along rows
                        trans_Y   = rot_transl(i_day_actual,4);
                        %translation X along columns
                        trans_X   = rot_transl(i_day_actual,5);
                        
                        %translation along rows
                        if trans_lr<0
                            Frame_OR(1:rw-trans_Y+1,:) = Frame_R(trans_Y:end,:);
                            trans_lr = -1;
                        elseif trans_lr>0
                            Frame_OR(trans_Y:end,:) = Frame_R(1:rw-trans_Y+1,:);
                            trans_lr = +1;
                        end
                        
                        %translation along columns
                        Frame_OR_2(:,1:cl-trans_X+1) = Frame_OR(:,trans_X:end);
                                                
                        %new coordinates of Bregma
                        XbN = 1;
                        YbN = 29; %==0.250 mm
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        Frame_ith = Frame_OR_2;
                        
%                         ALLFig = figure;
                        
                        for n_reg_col=1:round(Resol/Lat) %for # num regions-col in the Frame_ith // direction sx->dx
                            
                            index_col = n_reg_col-1;
                            
                            C_st = Lat*index_col+1;
                            C_en = C_st+Lat-1;
                            
                            for n_reg_row=1:round(Resol/Lat) %for # num regions-row in the Frame_ith // direction alto->basso
                                
                                index_row = n_reg_row-1;
                                
                                R_st = Lat*index_row+1;
                                R_en = R_st+Lat-1;
                                
                                ROI_SQUARED = Frame_ith( C_st:C_en, R_st:R_en);
                                
                                %%%
                                ROI_VALUE(n_reg_row, n_reg_col, inf) = median(median(ROI_SQUARED));
                                %%%
                                
                                %                                 figure
                                %                                 imagesc(ROI_SQUARED)
                                %                                 figure(ALLFig)
                                %                                 Frame_ith( R_st:R_en , C_st:C_en) = 1
                                %                                 imagesc(Frame_ith)
                                %
                            end
                            
                            
                        end
                        
                        
                    end
                    
                    
                    FF=figure
                    Seq_Delay=[];
                    inc = 0;
                    CentFrame = ImageSequence.CentFrame;
                    
                    for Roi_rw=1:size(ROI_VALUE,1)
                        
                        for Roi_cl=1:size(ROI_VALUE,2)
                            
                            figure(FF)
                            hold on
                            Seq = squeeze(ROI_VALUE(Roi_rw,Roi_cl,:));
                            SeqF = cheb2LPfilt(Seq,9,2,Fs);
                            
                            CheckSeq = SeqF-SeqF(1);
                            if all(~CheckSeq)
                                StartFluoPeak = NaN;
                                inc = inc+1;
                            else
                                
                                der_SeqF = derivative(SeqF,0.04);
                                [MAX_der_SeqF i_MAX_der_SeqF] = max(der_SeqF(1:length(der_SeqF)/1.5));
                                intPrePeak = 1:i_MAX_der_SeqF(1);
                                %%%% FLUO PEAK %%%%%%%%%
                                SeqBinary_der_SeqF      = der_SeqF(intPrePeak)>0;
                                SeqBinary_der_SeqF_Zero   = find(SeqBinary_der_SeqF==0);
                                
                                %nel caso il valore non sia proprio zero considero la derivata prima e vado a trovare il primo minimo dopo il valore massimo (== picco)) (che
                                %corrisponde all'ultimo valore del vettore der_IntervalFluoROI_Sig_Long in cui ho cambiato segno e tolto l'offset del valore di picco (in modo tale da
                                %usare la fun findpeaks che cerca i massimi)
                                if isempty(SeqBinary_der_SeqF_Zero) && length(intPrePeak)>2
                                    [pksValue pksLoc] = findpeaks( -(der_SeqF(intPrePeak) - der_SeqF(intPrePeak(end))));
                                    SeqBinary_der_SeqF_Zero = pksLoc;
                                elseif length(intPrePeak)<2
%                                     SeqBinary_der_SeqF_Zero = CentFrame;
                                    SeqBinary_der_SeqF_Zero = intPrePeak(end);
                                end
                                if isempty(SeqBinary_der_SeqF_Zero)
                                    SeqBinary_der_SeqF_Zero=1;
                                end
                                
                                StartFluoPeak                = SeqBinary_der_SeqF_Zero(end);       %[points]
%                                 StartFluoPeak = StartFluoPeak+1;                              
                                
                                SeqF = SeqF - SeqF(1);
                                
                                
%                                 plot( SeqF -5*inc);
%                                 hold on
%                                 plot( StartFluoPeak, SeqF(StartFluoPeak) -5*inc ,'-o')
                                
                                
                                inc = inc+1;
                                
                            end
                            
                            Seq_Delay = [Seq_Delay; (StartFluoPeak-CentFrame)/Fs];
                            
                        end
                    end
%                     hold on
%                     plot([CentFrame CentFrame], [5 -320],'r')
%                     ylim([-320 5])
%                     xlim([0 400])
                    
%                     figure
%                     plot(Seq_Delay,'-o');
                    
%                     figure
                    imagesc(reshape(Seq_Delay,size(ROI_VALUE,1),size(ROI_VALUE,2)))
                    colorbar
                    ImDelay(:,:,int) = reshape(Seq_Delay,size(ROI_VALUE,1),size(ROI_VALUE,2));
                    
                    
                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                

                figure
%                 ImDelay_Day = nanmedian(ImDelay,3);
                ImDelay_Day = nanmean(ImDelay,3);
                ImDelay_Day_M = medfilt2(ImDelay_Day);
                imagesc(ImDelay_Day_M)
                colorbar
                caxis([-0.3 0.3])
                saveas(gca,strcat(Animal_name,'_',nday,'_','IsoChronesFig_M',num2str(N_reg)),'fig')
                saveas(gca,strcat(Animal_name,'_',nday,'_','IsoChronesFig_M',num2str(N_reg)),'jpeg')
                close all                
               
                ImDelay_Day_M_All(:,:,index_day) = ImDelay_Day_M;
                
            end
        end
    end
end

             
                
                
                
% % %                 
% % %                 
% % %                 
% % %                 
% % %                 
% % %                 
% % %                 
% % %                 
% % %                 
% % %                 
% % %                 MIP      = max(MIP_all,[],3);
% % %                                 
% % %                 %assess threshold based on the most probable value
% % %                 Threshold = nanmedian(nanmedian(nanmedian((MIP)))) + 1*std2(MIP);
% % %                 
% % %                 
% % %                 %set origin (based on black dot)
% % %                 H_MIP_Fig = figure('Name',['MIP_Diff_TASK_Long_REST_ABS_',nday]);
% % %                 subplot(221)
% % %                 imagesc(MIP)
% % %                 colormap gray
% % %                 
% % %                 if 1
% % %                     
% % %                     %%%%%% same rotation and translation of the previous MIP%%%%%%%
% % %                     i_day_actual = find(rot_transl(:,1) == str2num(nday));
% % %                     
% % %                     %%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %                     degree = rot_transl(i_day_actual,2);
% % %                     if  degree~= 0
% % %                         Frame_R = imrotate(MIP,degree,'crop');
% % %                     else
% % %                         Frame_R = MIP;
% % %                     end
% % %                     
% % %                     Frame_OR = zeros(rw,cl);
% % %                     Frame_OR_2 = zeros(rw,cl);
% % %                     
% % %                     %%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %                     %left/right transl
% % %                     trans_lr  = rot_transl(i_day_actual,3);
% % %                     %translation Y along rows
% % %                     trans_Y   = rot_transl(i_day_actual,4);
% % %                     %translation X along columns
% % %                     trans_X   = rot_transl(i_day_actual,5);
% % %                     
% % %                     %translation along rows
% % %                     if trans_lr<0
% % %                         Frame_OR(1:rw-trans_Y+1,:) = Frame_R(trans_Y:end,:);
% % %                         trans_lr = -1;
% % %                     elseif trans_lr>0
% % %                         Frame_OR(trans_Y:end,:) = Frame_R(1:rw-trans_Y+1,:);
% % %                         trans_lr = +1;
% % %                     end
% % %                     
% % %                     %translation along columns
% % %                     Frame_OR_2(:,1:cl-trans_X+1) = Frame_OR(:,trans_X:end);
% % %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %                     
% % %                     %new coordinates of Bregma
% % %                     XbN = 1;
% % %                     YbN = 29; %==0.250 mm
% % %                     
% % %                     
% % %                 else
% % %                     %no alignement
% % %                     MIP_OR_2 = MIP;
% % %                     XbN = 1;
% % %                     YbN = 1;
% % %                 end
% % %                 
% % %                 
% % %                 figure(H_MIP_Fig)
% % %                 subplot(222)
% % %                 imagesc(Frame_OR_2)
% % %                 hold on
% % %                 text(XbN,YbN, 'X','Color','red','FontSize',10);
% % %                 changeLabel_MIP(rw,cl)
% % %                 colormap gray
% % %                 %%%%%%%%
% % %                 
% % %                 
% % %                 
% % %                 MIP = Frame_OR_2;
% % %                 
% % %                 
% % %                 %binary
% % %                 m_th = (1)/(max(max(MIP)-min(min(MIP))));
% % %                 q_th =  1 - m_th*max(max(MIP));
% % %                 MIP_preTh       = MIP.*m_th + q_th;
% % %                 Threshold_preTh = Threshold*m_th + q_th;
% % %                 MIP_Th = im2bw(MIP_preTh,Threshold_preTh);
% % %                 
% % %                 %         figure
% % %                 %         subplot(121)
% % %                 %         imagesc(MIP_max_filt_Th)
% % %                 
% % %                 %imerode
% % %                 se = strel('disk',2);
% % %                 MIP_Th_Dil = imerode(imdilate(MIP_Th,se),se);
% % %                 
% % %                 %filtering
% % %                 MIP_Th_Dil_filt = medfilt2(MIP_Th_Dil,[15 15]);
% % %                 
% % %                 
% % %                 figure(H_MIP_Fig)
% % %                 subplot(223)
% % %                 imagesc(MIP_Th_Dil_filt)
% % %                 hold on
% % %                 text(XbN,YbN, 'X','Color','red','FontSize',10);
% % %                 changeLabel_MIP(rw,cl)
% % %                 
% % %                 
% % %                 %find white pixels
% % %                 [wp_index]  = find(MIP_Th_Dil_filt == 1);
% % %                 [wp_X wp_Y] = ind2sub(size(MIP_Th_Dil_filt),wp_index);
% % %                 %find gray values corresponding to white pixels
% % %                 numTone     = MIP(wp_index);
% % %                 %max tones
% % %                 maxNumTone_i     = find(numTone == max(numTone));
% % %                 maxNumTone_index = wp_index(maxNumTone_i);
% % %                 [wp_X_max wp_Y_max] = ind2sub(size(MIP_Th_Dil_filt),maxNumTone_index);
% % %                 
% % %                 if ~isempty(wp_X_max)
% % %                     %remap intervals into 0-100 range
% % %                     max_numTone = max(numTone);
% % %                     min_numTone = min(numTone);
% % %                     m_t = (100-0)/(max_numTone(1)-min_numTone(1));
% % %                     q_t = 100-m_t*(max_numTone);
% % %                     numTone_R = m_t.*numTone + q_t;
% % %                     
% % %                     %             %multiply for the numTone
% % %                     %             wp_X_Tone = wp_X.*numTone_R;
% % %                     %             wp_Y_Tone = wp_Y.*numTone_R;
% % %                     %
% % %                     %             %find centroid coordinates (weighted on the gray tone)
% % %                     %             x = round(sum(wp_X_Tone)/sum(numTone_R));
% % %                     %             y = round(sum(wp_Y_Tone)/sum(numTone_R));
% % %                 else
% % %                     x = NaN;
% % %                     y = NaN;
% % %                     wp_index  = 1;
% % %                     numTone_R = 0;
% % %                 end
% % %                 CM           = zeros(rw,cl);
% % %                 CM(wp_index) = numTone_R;
% % %                 
% % %                 figure(H_MIP_Fig)
% % %                 subplot(224)
% % %                 imagesc(CM)
% % %                 colormap
% % %                 hold on
% % %                 text(XbN,YbN, 'X','Color','red','FontSize',10);
% % %                 changeLabel_MIP(rw,cl)
% % %                 %         hold on
% % %                 %         text(x,y, 'X','Color','red','FontSize',10);
% % %                 
% % %                 
% % %                 
% % %                 %
% % %                 H_MIP_Fig_filename = [Animal_name,'_',MIP_tech,'MIP_Diff_TASK_Long_REST_ABS_Fig_',nday];
% % %                 saveas(H_MIP_Fig,[DirDayToSaveMIP,'\',H_MIP_Fig_filename],'fig');
% % %                 saveas(H_MIP_Fig,[DirDayToSaveMIP,'\',H_MIP_Fig_filename],'tif');
% % %                 
% % %                 
% % %                 %all days together
% % %                 MIP_Th_Dil_filt_All_Day(:,:,index_day) = MIP_Th_Dil_filt;
% % %                 close(H_MIP_Fig)
% % %                 
% % %                 
% % %             end
% % %         end
% % %         
% % %         
% % %         MIP_All_Day_Fig_filename_week =  [Animal_name,'_',MIP_tech,'MIP_Diff_TASK_Long_REST_ABS_Overlapped_ROI_Fig_'];
% % %         MIP_All_Day_Fig = figure('Name',MIP_All_Day_Fig_filename_week);
% % %         
% % %         Im_med = sum(MIP_Th_Dil_filt_All_Day,3);
% % %         imagesc(Im_med)
% % %         hold on
% % %         text(XbN,YbN, 'X','Color','red','FontSize',10);
% % %         changeLabel_MIP(rw,cl)
% % %         saveas(MIP_All_Day_Fig,[MainDir,'\',MIP_All_Day_Fig_filename_week],'fig');
% % %         saveas(MIP_All_Day_Fig,[MainDir,'\',MIP_All_Day_Fig_filename_week],'tif');
% % %         
% % %         %save Im_med mat
% % %         save([MainDir,'\Im_AlongDays_',Animal_name,'_MIP_Diff_TASK_Long_REST_ABS_',MIP_tech],'Im_med')
% % %         
% % %         
% % %         %rehab group
% % %         if size(rot_transl,1)>5
% % %             
% % %             display([Animal_name,' is a rehab animal']);
% % %             
% % %             ActDay = rot_transl(:,1);
% % %             
% % %             %split days based on the week
% % %             wee = [(1:5);(6:10);(11:15);(16:20)];
% % %             for wee_i=1:size(wee,1)
% % %                 
% % %                 for wee_j=1:size(wee,2)
% % %                     
% % %                     i_wk = find(ActDay == wee(wee_i,wee_j));
% % %                     if ~isempty(i_wk)
% % %                         Im_buf(:,:,wee_j) = MIP_Th_Dil_filt_All_Day(:,:,i_wk);
% % %                     end
% % %                     
% % %                 end
% % %                 
% % %                 Im_med = sum(Im_buf,3);
% % %                 
% % %                 %save Im_med mat for each week
% % %                 MIP_All_Day_Fig_filename_week =  [Animal_name,'_',MIP_tech,'MIP_Diff_TASK_Long_REST_ABS_Week_',num2str(wee_i),'_Overlapped_ROI_Fig_'];
% % %                 MIP_All_Day_Fig = figure('Name',MIP_All_Day_Fig_filename_week);
% % %                 imagesc(Im_med)
% % %                 hold on
% % %                 text(XbN,YbN, 'X','Color','red','FontSize',10);
% % %                 changeLabel_MIP(rw,cl)
% % %                 saveas(MIP_All_Day_Fig,[MainDir,'\',MIP_All_Day_Fig_filename_week],'fig');
% % %                 saveas(MIP_All_Day_Fig,[MainDir,'\',MIP_All_Day_Fig_filename_week],'tif');
% % %                 
% % %                 save([MainDir,'\Im_AlongDays_Week_',num2str(wee_i),'_',Animal_name,'_MIP_Diff_TASK_Long_REST_ABS_',MIP_tech],'Im_med')
% % %                 
% % %                 clear Im_buf Im_med
% % %                 
% % %             end
% % %         end
% % %            
% % %         
% % %     end
% % %     
% % %     close all
% % %     clearvars -except lat ListAnimalTogether CurrDir
% % %     
display('End Process')
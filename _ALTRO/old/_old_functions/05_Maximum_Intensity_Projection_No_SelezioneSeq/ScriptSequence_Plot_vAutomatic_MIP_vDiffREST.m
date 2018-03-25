% Sequence of frames with a central one
% TrOK = [1:26];
% TrOK = [13];
% posPlot = [5 10 25 9];


clear all
close all
clc


%%%%%%%%%%%
% Animal_name = 'GCaMPChR2_7_control';
% Animal_name = 'GCaMPChR2_17_control';
% Animal_name = 'GCaMPChR2_18_control';
% Animal_name = 'GCaMPChR2_20_control';
% Animal_name = 'GCaMPChR2_21_control';
% Animal_name = 'GCaMPChR2_23_control';
% Animal_name = 'GCaMPChR2_24_control';

% Animal_name = 'GCaMPChR2_8_stroke';
% Animal_name = 'GCaMPChR2_9_stroke';
% Animal_name = 'GCaMPChR2_19_stroke';
% Animal_name = 'GCaMPChR2_22_stroke';
% Animal_name = 'GCaMPChR2_25_stroke';

% Animal_name = 'GCaMP16_stroke_BoNT';
% Animal_name = 'GCaMP17_stroke_BoNT';
% Animal_name = 'GCaMP17_stroke_BoNT';
% Animal_name = 'GCaMPChR2_3_stroke_BoNT';
% Animal_name = 'GCaMPChR2_11_stroke_BoNT';
% Animal_name = 'GCaMPChR2_12_stroke_BoNT';
% Animal_name = 'GCaMPChR2_13_stroke_BoNT';
% Animal_name = 'GCaMPChR2_14_stroke_BoNT';
% Animal_name = 'GCaMPChR2_15_stroke_BoNT';
% Animal_name = 'GCaMPChR2_16_stroke_BoNT';
%%%%%%%%%%%
CurrDir = cd;

ListAnimalTogether = {      'GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
                            'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke',...
                            'GCaMP16_stroke_BoNT','GCaMP17_stroke_BoNT','GCaMP18_stroke_BoNT',...
                            'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
                            'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};


for lat=1:length(ListAnimalTogether)
    
    UserName = 'CNR-SSSUP';
    UsbPort = 'I';
    % UserName = 'Stefano';
    % UsbPort = 'H';
    
    %%%%%%%%%%%
    Animal_name = ListAnimalTogether{lat};
    %%%%%%%%%%%
    
    %%%%%%%%%%%
    MIP_tech_List = {'max','sum'};
    %%%%%%%%%%%
    
    MainDir = [UsbPort,':\LENS\Animals Data\',Animal_name,'\'];
    NumDaysFolder = dir(MainDir);
    
    
    %%%%%%%% load Reference File %%%%%%%%%%%%%%%%%%%%%%%%
    RefDir      = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\Script_Flip_Find_References\MAT_Rot_Trans'];
    FileRefName = [Animal_name,'_Rot_Trans_Par.mat'];
    if exist([RefDir,'\',FileRefName])
        % rot_transl
        load([RefDir,'\',FileRefName]);
    else
        error([RefDir,'\',FileRefName,' is not present in the folder']);
    end
    
    
    MIP_Th_Dil_filt_All_Day = [];
    
    degree_all = [];
    trans_all  = [];
    
    for iMT=1:length(MIP_tech_List)
        
        MIP_tech = MIP_tech_List{iMT};
        
        
        for nd_i=3:length(NumDaysFolder)
            
            index_day = nd_i-2;
            nday_date = NumDaysFolder(nd_i,1).name;
            nday      = nday_date(1:2);
            
            
            DirDay      = [MainDir,nday_date,'\','Sequence'];
            DirDay_REST = [MainDir,nday_date,'\','Sequence_REST'];
            
            SeqDir      = dir(DirDay);
            SeqDir_REST = dir(DirDay_REST);
            
            if  ~isempty(SeqDir) && ~isempty(SeqDir)
                
                %seq ok
                SeqDirOKNO = 1;
                
                %load Image Sequence
                load([DirDay_REST,'\',SeqDir_REST(3,1).name])
                ImageSequence_REST = ImageSequence;
                clear ImageSequence
                load([DirDay,'\',SeqDir(3,1).name])
                
                %%% Initial Check %%%
                if ~exist('ImageSequence') | ~exist('ImageSequence_REST')
                    error('load ImageSequence')
                end
                %%%
                
                %%% Dir where saving
                DirDayToSaveMIP = [MainDir,nday_date,'\','MIP_Diff_TASK_REST_ABS'];
                if ~isdir(DirDayToSaveMIP)
                    %make Folder where saving
                    mkdir(DirDayToSaveMIP)
                end
                %%%
                
                %%%%%%%%% TASK  %%%%%%%%%%%%%
                MatrixImageSequence = ImageSequence.MatrixImageSequence;
                NumTrials = size(MatrixImageSequence,1);
                TrialsUsed    = [1:NumTrials];
                % TrialsUsed    = [1];
                NumTrialsUsed = length(TrialsUsed);
                rw        = size(MatrixImageSequence{1,1},1);
                cl        = size(MatrixImageSequence{1,1},2);
                NumFrames = size(MatrixImageSequence{1,1},3);
                
                %to store median/mean
                Med_Frames_Matrix = zeros(rw,cl,NumFrames);
                %to store median/mean -> threshold
                Med_Frames_Matrix_Th = zeros(rw,cl,NumFrames);
                %Fixed Threshold
                % Threshold = 0.05;
                
                MIP_all = zeros(rw, cl,  NumFrames*NumTrialsUsed);
                ind_i = 0;
                for int=1:NumTrialsUsed
                    for inf=1:NumFrames
                        ind_i = ind_i+1;
                        MIP_all(:,:,ind_i) = MatrixImageSequence{int,:}(:,:,inf);
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%%%%% REST %%%%%%%%%%%%%
                ImageSequence_REST = ImageSequence_REST.MatrixImageSequence;
                NumTrials = size(ImageSequence_REST,1);
                TrialsUsed    = [1:NumTrials];
                % TrialsUsed    = [1];
                NumTrialsUsed = length(TrialsUsed);
                rw        = size(ImageSequence_REST{1,1},1);
                cl        = size(ImageSequence_REST{1,1},2);
                NumFrames = ImageSequence.MinSize;
                
                %to store median/mean
                Med_Frames_Matrix = zeros(rw,cl,NumFrames);
                %to store median/mean -> threshold
                Med_Frames_Matrix_Th = zeros(rw,cl,NumFrames);
                %Fixed Threshold
                % Threshold = 0.05;
                
                MIP_all_REST = zeros(rw, cl,  NumFrames*NumTrialsUsed);
                ind_i = 0;
                for int=1:NumTrialsUsed
                    for inf=1:NumFrames
                        ind_i = ind_i+1;
                        MIP_all_REST(:,:,ind_i) = ImageSequence_REST{int,:}(:,:,inf);
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%%%%
                if strcmp(MIP_tech,'max')
                    %%%%%%%%
                    %%MAX
                    MIP_max      = max(MIP_all,[],3);
                    MIP_max_REST = max(MIP_all_REST,[],3);
                    
                    %                 Diff_MIP     = MIP_max-MIP_max_REST;
                    Diff_MIP     = abs(MIP_max)-abs(MIP_max_REST);
                    
                    minEdge = min(min(abs(MIP_max)));
                    
                    rstart = 1;
                    rstop  = 2;
                    
                    %                 %cut edges
                    %                 MIP_max(:,[rstart:rstop])    = minEdge;
                    %                 MIP_max(:,end-rstop:end)     = minEdge;
                    %
                    %                 MIP_max([rstart:rstop],:)    = minEdge;
                    %                 MIP_max(end-rstop:end,:)     = minEdge;
                    %
                    MIP = Diff_MIP;
                    
                    
                elseif strcmp(MIP_tech,'sum')
                    %%%%%%%%
                    %%SUM
                    MIP_sum      = sum(MIP_all,3);
                    MIP_sum_REST = sum(MIP_all_REST,3);
                    
                    %                 Diff_SIP     = MIP_sum-MIP_sum_REST;
                    Diff_SIP     = abs(MIP_sum)-abs(MIP_sum_REST);
                    
                    
                    minEdge = min(min(abs(MIP_sum)));
                    
                    rstart = 1;
                    rstop  = 2;
                    
                    %                 %cut edges
                    %                 MIP_sum(:,[rstart:rstop])    = minEdge;
                    %                 MIP_sum(:,end-rstop:end)     = minEdge;
                    %
                    %                 MIP_sum([rstart:rstop],:)    = minEdge;
                    %                 MIP_sum(end-rstop:end,:)     = minEdge;
                    %
                    MIP = Diff_SIP;
                    %
                    %%%%
                end
                %%%%%%%%
                
                
                %%% remove areas from MIP %%%
                %             cle = 1;
                %             frem = figure;
                %             while cle == 1
                %
                %                 imagesc(MIP)
                %                 colormap gray
                %                 ch = menu('Do you want to remove any area?','Remove it','Not, Go ahead');
                %                 switch ch
                %                     case 1
                %                         h_rect = imrect(gca, [10 10 30 30]); %-> finestra
                %                         pause
                %
                %                         rect_Remov = round(getPosition(h_rect));
                %                         MIP(rect_Remov(2):rect_Remov(2)+rect_Remov(4),rect_Remov(1):rect_Remov(1)+rect_Remov(3)) = 0;
                %
                %                     case 2
                %                         cle = 0;
                %                         close(frem)
                %                 end
                %             end
                
                %%%% ROIs %%%%%
                
                %             %filtering
                %             MIP = medfilt2(MIP);
                
                %assess threshold based on the most probable value
                Threshold = nanmedian(nanmedian(nanmedian((MIP)))) + 1*std2(MIP);
                
                
                %set origin (based on black dot)
                H_MIP_Fig = figure('Name',['MIP_Diff_TASK_REST_ABS_',nday]);
                subplot(221)
                imagesc(MIP)
                colormap gray
                
                if 1
                    
                    %%%%%% same rotation and translation of the previous MIP%%%%%%%
                    i_day_actual = find(rot_transl(:,1) == str2num(nday));
                    
                    %%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    degree = rot_transl(i_day_actual,2);
                    if  degree~= 0
                        MIP_R = imrotate(MIP,degree,'crop');
                    else
                        MIP_R = MIP;
                    end
                    
                    MIP_OR = zeros(rw,cl);
                    MIP_OR_2 = zeros(rw,cl);
                    
                    %%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %left/right transl
                    trans_lr  = rot_transl(i_day_actual,3);
                    %translation Y along rows
                    trans_Y   = rot_transl(i_day_actual,4);
                    %translation X along columns
                    trans_X   = rot_transl(i_day_actual,5);
                    
                    %translation along rows
                    if trans_lr<0
                        MIP_OR(1:rw-trans_Y+1,:) = MIP_R(trans_Y:end,:);
                        trans_lr = -1;
                    elseif trans_lr>0
                        MIP_OR(trans_Y:end,:) = MIP_R(1:rw-trans_Y+1,:);
                        trans_lr = +1;
                    end
                    
                    %translation along columns
                    MIP_OR_2(:,1:cl-trans_X+1) = MIP_OR(:,trans_X:end);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %new coordinates of Bregma
                    XbN = 1;
                    YbN = 29; %==0.250 mm
                    
                    
                else
                    %no alignement
                    MIP_OR_2 = MIP;
                    XbN = 1;
                    YbN = 1;
                end
                
                
                figure(H_MIP_Fig)
                subplot(222)
                imagesc(MIP_OR_2)
                hold on
                text(XbN,YbN, 'X','Color','red','FontSize',10);
                changeLabel_MIP(rw,cl)
                colormap gray
                %%%%%%%%
                
                
                
                MIP = MIP_OR_2;
                
                
                %binary
                m_th = (1)/(max(max(MIP)-min(min(MIP))));
                q_th =  1 - m_th*max(max(MIP));
                MIP_preTh       = MIP.*m_th + q_th;
                Threshold_preTh = Threshold*m_th + q_th;
                MIP_Th = im2bw(MIP_preTh,Threshold_preTh);
                
                %         figure
                %         subplot(121)
                %         imagesc(MIP_max_filt_Th)
                
                %imerode
                se = strel('disk',2);
                MIP_Th_Dil = imerode(imdilate(MIP_Th,se),se);
                
                %filtering
                MIP_Th_Dil_filt = medfilt2(MIP_Th_Dil,[15 15]);
                
                
                figure(H_MIP_Fig)
                subplot(223)
                imagesc(MIP_Th_Dil_filt)
                hold on
                text(XbN,YbN, 'X','Color','red','FontSize',10);
                changeLabel_MIP(rw,cl)
                
                
                %find white pixels
                [wp_index]  = find(MIP_Th_Dil_filt == 1);
                [wp_X wp_Y] = ind2sub(size(MIP_Th_Dil_filt),wp_index);
                %find gray values corresponding to white pixels
                numTone     = MIP(wp_index);
                %max tones
                maxNumTone_i     = find(numTone == max(numTone));
                maxNumTone_index = wp_index(maxNumTone_i);
                [wp_X_max wp_Y_max] = ind2sub(size(MIP_Th_Dil_filt),maxNumTone_index);
                
                if ~isempty(wp_X_max)
                    %remap intervals into 0-100 range
                    max_numTone = max(numTone);
                    min_numTone = min(numTone);
                    m_t = (100-0)/(max_numTone(1)-min_numTone(1));
                    q_t = 100-m_t*(max_numTone);
                    numTone_R = m_t.*numTone + q_t;
                    
                    %             %multiply for the numTone
                    %             wp_X_Tone = wp_X.*numTone_R;
                    %             wp_Y_Tone = wp_Y.*numTone_R;
                    %
                    %             %find centroid coordinates (weighted on the gray tone)
                    %             x = round(sum(wp_X_Tone)/sum(numTone_R));
                    %             y = round(sum(wp_Y_Tone)/sum(numTone_R));
                else
                    x = NaN;
                    y = NaN;
                    wp_index  = 1;
                    numTone_R = 0;
                end
                CM           = zeros(rw,cl);
                CM(wp_index) = numTone_R;
                
                figure(H_MIP_Fig)
                subplot(224)
                imagesc(CM)
                colormap
                hold on
                text(XbN,YbN, 'X','Color','red','FontSize',10);
                changeLabel_MIP(rw,cl)
                %         hold on
                %         text(x,y, 'X','Color','red','FontSize',10);
                
                
                
                %
                H_MIP_Fig_filename = [Animal_name,'_',MIP_tech,'MIP_Diff_TASK_REST_ABS_Fig_',nday];
                saveas(H_MIP_Fig,[DirDayToSaveMIP,'\',H_MIP_Fig_filename],'fig');
                saveas(H_MIP_Fig,[DirDayToSaveMIP,'\',H_MIP_Fig_filename],'tif');
                
                
                %all days together
                MIP_Th_Dil_filt_All_Day(:,:,index_day) = MIP_Th_Dil_filt;
                close(H_MIP_Fig)
                
                
            end
        end
        
        
        MIP_All_Day_Fig_filename_week =  [Animal_name,'_',MIP_tech,'MIP_Diff_TASK_REST_ABS_Overlapped_ROI_Fig_'];
        MIP_All_Day_Fig = figure('Name',MIP_All_Day_Fig_filename_week);
        
        Im_med = sum(MIP_Th_Dil_filt_All_Day,3);
        imagesc(Im_med)
        hold on
        text(XbN,YbN, 'X','Color','red','FontSize',10);
        changeLabel_MIP(rw,cl)
        saveas(MIP_All_Day_Fig,[MainDir,'\',MIP_All_Day_Fig_filename_week],'fig');
        saveas(MIP_All_Day_Fig,[MainDir,'\',MIP_All_Day_Fig_filename_week],'tif');
        
        %save Im_med mat
        save([MainDir,'\Im_AlongDays_',Animal_name,'_MIP_Diff_TASK_REST_ABS_',MIP_tech],'Im_med')
        
        
        %rehab group
        if size(rot_transl,1)>5
            
            display([Animal_name,' is a rehab animal']);
            
            ActDay = rot_transl(:,1);
            
            %split days based on the week
            wee = [(1:5);(6:10);(11:15);(16:20)];
            for wee_i=1:size(wee,1)
                
                for wee_j=1:size(wee,2)
                    
                    i_wk = find(ActDay == wee(wee_i,wee_j));
                    if ~isempty(i_wk)
                        Im_buf(:,:,wee_j) = MIP_Th_Dil_filt_All_Day(:,:,i_wk);
                    end
                    
                end
                
                Im_med = sum(Im_buf,3);
                
                %save Im_med mat for each week
                MIP_All_Day_Fig_filename_week =  [Animal_name,'_',MIP_tech,'MIP_Diff_TASK_REST_ABS_Week_',num2str(wee_i),'_Overlapped_ROI_Fig_'];
                MIP_All_Day_Fig = figure('Name',MIP_All_Day_Fig_filename_week);
                imagesc(Im_med)
                hold on
                text(XbN,YbN, 'X','Color','red','FontSize',10);
                changeLabel_MIP(rw,cl)
                saveas(MIP_All_Day_Fig,[MainDir,'\',MIP_All_Day_Fig_filename_week],'fig');
                saveas(MIP_All_Day_Fig,[MainDir,'\',MIP_All_Day_Fig_filename_week],'tif');
                
                save([MainDir,'\Im_AlongDays_Week_',num2str(wee_i),'_',Animal_name,'_MIP_Diff_TASK_REST_ABS_',MIP_tech],'Im_med')
                
                clear Im_buf Im_med
                
            end
        end
           
        
    end
    
    close all
    clearvars -except lat ListAnimalTogether CurrDir
    
end
display('End Process')
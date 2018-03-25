% Script to make Isochrones from Sequence_LONG
clear all
close all
clc

%%%%%%%%%%%
CurrDir = cd;

ListAnimalTogether = {      'GCaMPChR2_7_control'};

% ListAnimalTogether = {      'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                             'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke',...
%                             'GCaMP16_stroke_BoNT','GCaMP17_stroke_BoNT','GCaMP18_stroke_BoNT',...
%                             'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                             'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};

for lat=1:length(ListAnimalTogether)
    
        UserName = 'CNR-SSSUP';
        UsbPort = 'I';
%     UserName = 'Stefano';
%     UsbPort = 'H';
    
    %%%%REGION OF INTEREST
    Resol     = 512;
    Resol_SQ  = Resol^2;
    N_reg     = 1024;                 %total number of region
%         N_reg     = 128;                 %total number of region
    N_reg_Per_Border = sqrt(N_reg); %number of region for each region
    Lat       = Resol/N_reg_Per_Border; %Size of the border of each region
    %     Lat       = sqrt(Resol_SQ/N_reg); %Size of the border of each region
    %%%%
    
    %%% Frequency
    Fs = 25;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Animal_name = ListAnimalTogether{lat};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    for nd_i=length(NumDaysFolder):length(NumDaysFolder) %for days -> da togliere solo prova
%     for nd_i=3:length(NumDaysFolder) %for days
        
        index_day = nd_i-2;
        nday_date = NumDaysFolder(nd_i,1).name;
        nday      = nday_date(1:2);
        
        %Sequence of frames (all of frames are taken in SequenceLong)
        DirDay      = [MainDir,nday_date,'\','SequenceLong'];
        SeqDir      = dir(DirDay);
        
        
        if  ~isempty(SeqDir) %if SeqDir
            
            %seq ok
            SeqDirOKNO = 1;
            
            %load Image Sequence
            load([DirDay,'\',SeqDir(3,1).name])
            %%%%%%%%%%%%%%%%%%%%
            
            %%% Check %%%
            if ~exist('ImageSequence')
                error('load ImageSequence')
            end
            
            %%%%%%%%%%%%
            
            %                 %%% Dir where saving
            %                 DirDayToSaveMIP = [MainDir,nday_date,'\',''];
            %                 if ~isdir(DirDayToSaveMIP)
            %                     %make Folder where saving
            %                     mkdir(DirDayToSaveMIP)
            %                 end
            %                 %%%
            
            
            %%%%%%%%% DAILY TASK SESSION %%%%%%%%%%%%%
            MatrixImageSequence = ImageSequence.MatrixImageSequence;
            
            MatrixImageSequence = MatrixImageSequence([1 2 3 7 9 18],1);     %for days -> da togliere solo prova 
            
            NumTrials           = size(MatrixImageSequence,1);
            TrialsUsed          = [1:NumTrials];

            NumTrialsUsed       = length(TrialsUsed);
            rw                  = size(MatrixImageSequence{1,1},1);
            cl                  = size(MatrixImageSequence{1,1},2);
            NumFrames           = size(MatrixImageSequence{1,1},3);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
            ind_i = 0;
            ImDelay = [];
            
            for int=1:NumTrialsUsed %for trials
                
                %%initialization of the image with a total number of region = N_reg. 
                %%The third dimension of ROI_VALUE are the frames related to that trial
                ROI_VALUE = zeros(round(Resol/Lat),round(Resol/Lat),NumFrames);
                
                for inf=1:NumFrames %for frames in the trial i-th
                   
                    
                    Frame_ith = MatrixImageSequence{int,:}(:,:,inf);
                     
                    %%roto-translation of the frame i-th
                    Frame_ith = RotoTrans_Image( Frame_ith, rot_transl, nday, rw, cl);
                     
                    %%fill ROI_VALUE for that trials (3rd dim -> frames)                   
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
                            
% %                             figure
% %                             imagesc(ROI_SQUARED)
%                             Frame_ith_Buf = Frame_ith*0;
%                             Frame_ith_Buf( R_st:R_en , C_st:C_en) = 1;
%                             figure
%                             imagesc(Frame_ith_Buf)
%                             pause(0.05);
%                             close

                            
                        end %end for # num regions-row in the Frame_ith // direction alto->basso                       
                        
                    end %end for # num regions-col in the Frame_ith // direction sx->dx                    
                    
                end % end for frames in the trial i-th
                
                
%                 FF=figure

                %%%
                Seq_Delay = [];
                inc         = 0;
                CentFrame = ImageSequence.CentFrame;
                %%%
                
                for Roi_rw=1:size(ROI_VALUE,1) %for ROI row
                    
                    for Roi_cl=1:size(ROI_VALUE,2) %for ROI col
                        
%                         figure(FF)
%                         hold on
                        Seq  = squeeze(ROI_VALUE(Roi_rw,Roi_cl,:));
                        SeqF = cheb2LPfilt(Seq,9,2,Fs);
                        
                        CheckSeq = SeqF-SeqF(1);
                        if all(~CheckSeq)
                            StartFluoPeak = NaN;
                            inc = inc+1;
                        else
                            
                            der_SeqF = derivative(SeqF,0.04);
                            %%find max of the der_Seq_F in the first part of the interval 
                            [MAX_der_SeqF i_MAX_der_SeqF] = max(der_SeqF(1:length(der_SeqF)/1.5));
                            intPrePeak                    = 1:i_MAX_der_SeqF(1);
                            
                            %%%% FLUO PEAK %%%%%%%%%
                            SeqBinary_der_SeqF        = der_SeqF(intPrePeak)>0;
                            SeqBinary_der_SeqF_Zero   = find(SeqBinary_der_SeqF==0);
                            
                            %%nel caso il valore non sia proprio zero considero la derivata prima e vado a trovare il primo minimo dopo il valore massimo (== picco)) (che
                            %%corrisponde all'ultimo valore del vettore der_IntervalFluoROI_Sig_Long in cui ho cambiato segno e tolto l'offset del valore di picco (in modo tale da
                            %%usare la fun findpeaks che cerca i massimi)
                            if ~isempty(SeqBinary_der_SeqF_Zero)
                                
                                StartFluoPeak                = SeqBinary_der_SeqF_Zero(end);       %[points]
                                
                            elseif isempty(SeqBinary_der_SeqF_Zero)&& length(intPrePeak)>2
                                
                                [pksValue pksLoc] = findpeaks( -(der_SeqF(intPrePeak) - der_SeqF(intPrePeak(end))));
                                SeqBinary_der_SeqF_Zero = pksLoc;
                                
                                if ~isempty(SeqBinary_der_SeqF_Zero)
                                    StartFluoPeak                = SeqBinary_der_SeqF_Zero(end);       %[points]
%                                     StartFluoPeak = StartFluoPeak+1;
                                else
                                    StartFluoPeak = NaN;
                                end
                                   
                            else
                                StartFluoPeak = NaN;
                            end
                            
                            inc = inc+1;
                            
%                             %%%%%plot
%                             SeqF = SeqF - SeqF(1);
%                             plot( SeqF -5*inc);
%                             hold on
%                             if ~isnan(StartFluoPeak)
%                                 plot( StartFluoPeak, SeqF(StartFluoPeak) -5*inc ,'-o')
%                             end
%                             %%%%%%
                            
                        end
                        
                        %save the delay (onset force-onset fluo) for each ROI
                        Seq_Delay = [Seq_Delay; (StartFluoPeak-CentFrame)/Fs];
                        
                    end %end for ROI col
                    
                end %end for ROI row
                
%                 %%%%%%plot
%                 hold on
%                 plot([CentFrame CentFrame], [5 -320],'r')
%                 ylim([-320 5])
%                 xlim([0 400])
% %                 
% %                 figure
% %                 plot(Seq_Delay,'-o');
% %                 
%                 figure
%                 imagesc(reshape(Seq_Delay,size(ROI_VALUE,1),size(ROI_VALUE,2)))
%                 colorbar
%                 %%%%%%
                
                
                %save the image -> in each ROI(square) is saved the value of the Delay (for the j-th Trial)
                ImDelay(:,:,int) = reshape(Seq_Delay,size(ROI_VALUE,1),size(ROI_VALUE,2));
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            figure
%             ImDelay_Day = nanmedian(ImDelay,3);
            ImDelay_Day = nanmean(ImDelay,3);
            ImDelay_Day_M = medfilt2(ImDelay_Day);
            imagesc(ImDelay_Day_M)
            colorbar
            caxis([-0.3 0.3])
            saveas(gca,strcat(Animal_name,'_',nday,'_','IsoChronesFig_M',num2str(N_reg)),'fig')
            saveas(gca,strcat(Animal_name,'_',nday,'_','IsoChronesFig_M',num2str(N_reg)),'jpeg')
            close all
            
            
            %store the average value of the Delay (for each squared ROI) in the z-th day
            ImDelay_Day_M_All(:,:,index_day) = ImDelay_Day_M;
            
            clear ImDelay ImDelay_Day ImDelay_Day_M
            
        end %end if SeqDir
        
    end %end for days
    
    %save all info
    save(['DataIm_Mean_Delay_All_Days_',Animal_name],'ImDelay_Day_M_All');
    clear ImDelay_Day_M_All;
    
end

display('End Process')
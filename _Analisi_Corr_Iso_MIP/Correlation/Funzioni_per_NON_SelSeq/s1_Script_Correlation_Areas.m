% Script to make Correlation among functional areas (the sequences of
% images used to compute the fluo signal area the SequenceLong

tic

clear all
close all
clc

SEL_PART_LIST = {'TASK','REST'};
% SEL_PART_LIST = {'TASK'};

for i_selPART = 1:2 %for SEL_PART_LIST 


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SELECT SEPARATION IN:
%%% DAY:  to compute corr among areas for every single day (segnali messi uno dopo l'altro per avere una unica sequenza)
%%% WEEK: to compute corr among areas pooling the days of the same week (segnali di tuti i giorni messi uno dopo l'altro per avere una unica sequenza)
% SEL_SEP = 'DAY';
% SEL_SEP = 'WEEK';
%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% SELECT PART of the signal:
% %%% TASK: to compute corr among areas in the part of the signal corresponding to the task
% %%% REST:to compute corr among areas in the part of the signal corresponding to resting phase
% % SEL_PART = 'TASK';
% SEL_PART = 'REST';
SEL_PART = SEL_PART_LIST{i_selPART};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
CurrDir = cd;
% ListAnimalTogether = {  'GCaMPChR2_7_control'};
% ListAnimalTogether = {  'GCaMPChR2_11_stroke_BoNT'};

% ListAnimalTogether = {      'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                             'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke','GCaMPChR2_26_stroke',...
%                             'GCaMP16_stroke_BoNT','GCaMP17_stroke_BoNT','GCaMP18_stroke_BoNT',...
%                             'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                             'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};
                        
ListAnimalTogether = {     'GCaMPChR2_12_stroke_BoNT', 'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};


%%%%%%%% User %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UserName = 'CNR-SSSUP';
% UsbPort = 'K';
UserName = 'Stefano';
UsbPort = 'H';

% %%%%%%%% Where SAVE Folder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main_SAVE_Dir = [UsbPort,':\LENS\Correlazioni\'];
% SAVE_Dir      = [Main_SAVE_Dir,'Corr_Areas_','SEP',SEL_SEP,'_PART',SEL_PART];
% SAVE_Dir_Method1 = [SAVE_Dir,'\Method1\'];
% SAVE_Dir_Method2 = [SAVE_Dir,'\Method2\'];
%
% if ~isdir(SAVE_Dir_Method1) && ~isdir(SAVE_Dir_Method2)
%     %make Folder where saving Method 1
%     mkdir(SAVE_Dir_Method1)
%     %make Folder where saving Method 2
%     mkdir(SAVE_Dir_Method2)
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% Where SAVE Folder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PER FARLI INSIEME %%%%
Main_SAVE_Dir = [UsbPort,':\LENS\Correlazioni\'];
SEL_SEP = 'DAY';
SAVE_Dir      = [Main_SAVE_Dir,'Corr_Areas_','SEP',SEL_SEP,'_PART',SEL_PART];
SAVE_Dir_Method1 = [SAVE_Dir,'\Method1\'];
SAVE_Dir_Method2 = [SAVE_Dir,'\Method2\'];
SAVE_Dir_Method3 = [SAVE_Dir,'\Method3\'];

if ~isdir(SAVE_Dir_Method1) && ~isdir(SAVE_Dir_Method2) && ~isdir(SAVE_Dir_Method3)
    %make Folder where saving Method 1
    mkdir(SAVE_Dir_Method1)
    %make Folder where saving Method 2
    mkdir(SAVE_Dir_Method2)
    %make Folder where saving Method 3
    mkdir(SAVE_Dir_Method3)
end

SEL_SEP = 'WEEK';
SAVE_Dir      = [Main_SAVE_Dir,'Corr_Areas_','SEP',SEL_SEP,'_PART',SEL_PART];
SAVE_Dir_Method1_W = [SAVE_Dir,'\Method1\'];
SAVE_Dir_Method2_W = [SAVE_Dir,'\Method2\'];
SAVE_Dir_Method3_W = [SAVE_Dir,'\Method3\'];

if ~isdir(SAVE_Dir_Method1_W) && ~isdir(SAVE_Dir_Method2_W) && ~isdir(SAVE_Dir_Method3_W)
    %make Folder where saving Method 1
    mkdir(SAVE_Dir_Method1_W)
    %make Folder where saving Method 2
    mkdir(SAVE_Dir_Method2_W)
    %make Folder where saving Method 3
    mkdir(SAVE_Dir_Method3_W)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for lat=1:length(ListAnimalTogether)
    
    %%% Frequency
    Fs = 25;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Animal_name = ListAnimalTogether{lat};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% Data Folder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MainDir = [UsbPort,':\LENS\Animals Data\',Animal_name,'\'];
    NumDaysFolder = dir(MainDir);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% Data SAVE ANIMAL Folder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SAVE_Dir_Method1_An = [SAVE_Dir_Method1,  Animal_name,'\'];
    mkdir(SAVE_Dir_Method1_An)
    SAVE_Dir_Method2_An = [SAVE_Dir_Method2,  Animal_name,'\'];
    mkdir(SAVE_Dir_Method2_An)
    SAVE_Dir_Method3_An = [SAVE_Dir_Method3,  Animal_name,'\'];
    mkdir(SAVE_Dir_Method3_An)   
    
    SAVE_Dir_Method1_W_An = [SAVE_Dir_Method1_W,  Animal_name,'\'];
    mkdir(SAVE_Dir_Method1_W_An)
    SAVE_Dir_Method2_W_An = [SAVE_Dir_Method2_W,  Animal_name,'\'];
    mkdir(SAVE_Dir_Method2_W_An)
    SAVE_Dir_Method3_W_An = [SAVE_Dir_Method3_W,  Animal_name,'\'];
    mkdir(SAVE_Dir_Method3_W_An)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% load Reference File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RefDir      = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\Script_Flip_Find_References\MAT_Rot_Trans'];
    FileRefName = [Animal_name,'_Rot_Trans_Par.mat'];
    if exist([RefDir,'\',FileRefName])
        % rot_transl
        load([RefDir,'\',FileRefName]);
        list_real_days = rot_transl(:,1);
    else
        error([RefDir,'\',FileRefName,' is not present in the folder']);
    end
    degree_all = [];
    trans_all  = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% load the list of indeces corresponding to the functional areas%
    load('RegionArea_FunctionalMask_Index'); %-> RegionArea_FunctionalMask_Index
    FuncAres_Name = {'M2','CFA','M1Sh','M1HL','M1Tk','S1','S1Bc','RSD','LPTa','V1','AuD'};
    Num_FunAreas  = size(RegionArea_FunctionalMask_Index,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% Info Animal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [weekly_day_stop num_week] = Fun_Info_Animal(Animal_name,list_real_days);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% for SEL_SEP = 'DAY'
    Mat_corr_1_Tot = zeros(Num_FunAreas,Num_FunAreas,length(NumDaysFolder)-2);
    Mat_corr_2_Tot = zeros(Num_FunAreas,Num_FunAreas,length(NumDaysFolder)-2);
    Mat_corr_3_Tot = zeros(Num_FunAreas,Num_FunAreas,length(NumDaysFolder)-2);
    %%% for SEL_SEP = 'WEEK'
    Mat_corr_1_W_Tot = zeros(Num_FunAreas,Num_FunAreas,num_week);
    Mat_corr_2_W_Tot = zeros(Num_FunAreas,Num_FunAreas,num_week);
    Mat_corr_3_W_Tot = zeros(Num_FunAreas,Num_FunAreas,num_week);
    
    nw_ind = 0;
    
    NumFramesPrevDay = [];
    
    for nd_i=3:length(NumDaysFolder) %for days
        
        index_day = nd_i-2;
        nday_date = NumDaysFolder(nd_i,1).name;
        nday      = nday_date(1:2);
        
        %Sequence of frames (all of frames are taken in SequenceLong or in Sequence_REST)
        if strcmp(SEL_PART,'TASK')
            DirDay      = [MainDir,nday_date,'\','SequenceLong'];
        elseif strcmp(SEL_PART,'REST')
            DirDay      = [MainDir,nday_date,'\','Sequence_REST'];
        end
        SeqDir          = dir(DirDay);
        
        
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
            
            %%%%%%%%% DAILY TASK SESSION %%%%%%%%%%%%%
            MatrixImageSequence = ImageSequence.MatrixImageSequence;
            NumTrials           = size(MatrixImageSequence,1);
            TrialsUsed          = [1:NumTrials];
            
            NumTrialsUsed       = length(TrialsUsed);
            rw                  = size(MatrixImageSequence{1,1},1);
            cl                  = size(MatrixImageSequence{1,1},2);
            if strcmp(SEL_PART,'TASK')
                NumFrames       = size(MatrixImageSequence{1,1},3);
            elseif strcmp(SEL_PART,'REST')
                NumFrames       = ImageSequence.MinSize;
            end
            %number of frame for the Sequence in the previous day
            if isempty(NumFramesPrevDay)
                NumFramesPrevDay    = NumFrames;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ind_i = 0;
            
            %store day
            FuncAreaValue_Day_Met1 = [];
            FuncAreaValue_Day_Met2 = [];
            FuncAreaValue_Day_Met3 = [];
            %store week
            FuncAreaValue_Week_Met1 = [];
            FuncAreaValue_Week_Met2 = [];
            FuncAreaValue_Week_Met3 = [];
            
            for int=1:NumTrialsUsed %for trials
                
                FuncAreaValue_Tr = zeros(Num_FunAreas,NumFrames);
                for inf=1:NumFrames %for frames in the trial i-th
                    
                    Frame_ith = MatrixImageSequence{int,:}(:,:,inf);
                    
                    %%roto-translation of the frame i-th
                    Frame_ith = RotoTrans_Image( Frame_ith, rot_transl, nday, rw, cl);
                    [s1 s2 s3] = size(Frame_ith);
                    
                    %%suddivisione in func aree e storage del valore mediano
                    FuncAreaValue_singleFr = zeros(Num_FunAreas,1);
                    for i_funA = 1:Num_FunAreas
                        
                        FuncAreaValue_singleFr(i_funA) = nanmedian(Frame_ith(sub2ind([s1 s2],RegionArea_FunctionalMask_Index{i_funA,1}(:,2),RegionArea_FunctionalMask_Index{i_funA,1}(:,1))));
                        
                        %%%
                        figure
                        Frame_ith_PLOT = Frame_ith;
                        Frame_ith_PLOT(sub2ind([s1 s2],RegionArea_FunctionalMask_Index{i_funA,1}(:,2),RegionArea_FunctionalMask_Index{i_funA,1}(:,1))) = 100;
                        imagesc(Frame_ith_PLOT)
                        ciao = 1;
                        close
                        %%%
                        
                    end
                    
                    %%Metodo 1 -> correlare il segnale dei trials così
                    %%come è (i trials sono messi uno dopo l'altro per avere
                    %%un lungo segnale)
                    FuncAreaValue_Tr(:,inf) = FuncAreaValue_singleFr;
                    
                    
                end % end for frames in the trial i-th
                
                %%Metodo 2 -> correlare il segnale dei trials
                %%identificando onset picco fluo. Per ogni trial ho
                %%quindi un segnale tutto 0 eccetto nell'istante di
                %%onset del picco fluo. Anche in questo caso i segnali
                %%così fatti sono messi uno dopo l'altro per avere
                %%un lungo segnale)
                FluoOnset_Matrix_Tr = zeros(Num_FunAreas,NumFrames);
                %%Metodo 3 -> correlare il segnale dei trials
                %%identificando onset picco fluo, considerando il delay(numero)
                Delay_Tr            = zeros(Num_FunAreas,1);
                
                for i_funA = 1:Num_FunAreas %for 2nd method and %for 3rd method
                    
                    Seq = FuncAreaValue_Tr(i_funA,:);
                    SeqF = cheb2LPfilt(Seq,9,2,Fs);
                    
                    CheckSeq = SeqF-SeqF(1);
                    if all(~CheckSeq) %if checkSeq
                        StartFluoPeak = [];
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
                            else
                                StartFluoPeak = [];
                            end
                        else
                            StartFluoPeak = [];
                        end
                        
                    end %if checkSeq
                    
                    if ~isempty(StartFluoPeak)
                        %Method 2 -> array
                        FluoOnset_Matrix_Tr(i_funA,StartFluoPeak) = 1;
                        %Method 3 -> delay (scalar)
                        Delay_Tr(i_funA,1) = StartFluoPeak;
                    end
                    
                end %end for 2nd method %end for 3rd method
                %%%%%%%%%%%%%%%%
                
                
                %store signal di seguito (singolo DAY)
                FuncAreaValue_Day_Met1 = [FuncAreaValue_Day_Met1 FuncAreaValue_Tr];
                FuncAreaValue_Day_Met2 = [FuncAreaValue_Day_Met2 FluoOnset_Matrix_Tr];
                FuncAreaValue_Day_Met3 = [FuncAreaValue_Day_Met3 Delay_Tr];
                
                
                %                 %%% plot signal
                %                 if int==NumTrialsUsed
                %                     fig_M1 = figure('Name','Method 1');
                %                     fig_M2 = figure('Name','Method 2');
                %                     for fM=1:Num_FunAreas
                %                         figure(fig_M1)
                %                         plot(FuncAreaValue_Day_Met1(fM,:)+ fM*3)
                %                         hold on
                %                         figure(fig_M2)
                %                         plot(FuncAreaValue_Day_Met2(fM,:)+ fM*3)
                %                         hold on
                %                     end
                %                 end
                %                 %%%%%%
                
            end %end trials
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%% SEL_SEP separati
            %             if strcmp(SEL_SEP,'DAY') %if SEL_SEP save Fig
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %%% CORRELAZIONE %%%
            %                 [Mat_corr_1 Mat_corr_2] = Fun_Compute_Corr(FuncAreaValue_Day_Met1, FuncAreaValue_Day_Met2, Num_FunAreas);
            %
            %                 %store CORR matrix
            %                 Mat_corr_1_Tot(:,:,nd_i-2) =  Mat_corr_1;
            %                 Mat_corr_2_Tot(:,:,nd_i-2) =  Mat_corr_2;
            %
            %                 %save CORR image
            %                 Fun_Save_Im_Corr(Animal_name, SAVE_Dir_Method1_An, SAVE_Dir_Method2_An, SEL_SEP, SEL_PART, Mat_corr_1, Mat_corr_2, nd_i);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %             elseif strcmp(SEL_SEP,'WEEK')
            %
            %                 %%numero di punti nella sequenza del giorno corrente deveessere uguale a quella dei giorni precedenti (si prende la dimensione più piccola)
            %                 diff_Frames = NumFramesPrevDay - NumFrames;
            %                 if diff_Frames > 0
            %                     FuncAreaValue_Week_Met1 = FuncAreaValue_Week_Met1(1:NumFrames,:);
            %                     FuncAreaValue_Week_Met2 = FuncAreaValue_Week_Met2(1:NumFrames,:);
            %                     NumFramesPrevDay = NumFrames;
            %                 elseif diff_Frames < 0
            %                     FuncAreaValue_Day_Met1  = FuncAreaValue_Day_Met1(1:NumFramesPrevDay,:);
            %                     FuncAreaValue_Day_Met2  = FuncAreaValue_Day_Met2(1:NumFramesPrevDay,:);
            %                     NumFramesPrevDay = NumFramesPrevDay;
            %                 end
            %
            %                 FuncAreaValue_Week_Met1 = [FuncAreaValue_Week_Met1 FuncAreaValue_Day_Met1];
            %                 FuncAreaValue_Week_Met2 = [FuncAreaValue_Week_Met2 FuncAreaValue_Day_Met2];
            %
            %                 %last day of the week of interest
            %                 if num2str(nday) == weekly_day_stop
            %
            %                     nw_ind = nw_ind+1;
            %
            %                     %%%%%%%%%%%%%%%%%%%%%%%%%%
            %                     %%% CORRELAZIONE %%%
            %                     [Mat_corr_1_W Mat_corr_2_W] = Fun_Compute_Corr(FuncAreaValue_Week_Met1, FuncAreaValue_Week_Met2, Num_FunAreas);
            %
            %                     %store CORR matrix
            %                     Mat_corr_1_W_Tot(:,:,nw_ind) =  Mat_corr_1_W;
            %                     Mat_corr_2_W_Tot(:,:,nw_ind) =  Mat_corr_2_W;
            %
            %                     %save CORR image
            %                     Fun_Save_Im_Corr(Animal_name, SAVE_Dir_Method1_An, SAVE_Dir_Method2_An, SEL_SEP, SEL_PART, Mat_corr_1_W, Mat_corr_2_W, nw_ind);
            %                     %%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %
            %                 end
            %
            %             end %end if SEL_SEP save Fig
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% SEL_SEP insieme
%             if strcmp(SEL_SEP,'DAY') %if SEL_SEP save Fig
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% CORRELAZIONE %%%
                [Mat_corr_1 Mat_corr_2 Mat_corr_3] = Fun_Compute_Corr(FuncAreaValue_Day_Met1, FuncAreaValue_Day_Met2, FuncAreaValue_Day_Met3, Num_FunAreas);
                
                %store CORR matrix
                Mat_corr_1_Tot(:,:,nd_i-2) =  Mat_corr_1;
                Mat_corr_2_Tot(:,:,nd_i-2) =  Mat_corr_2;
                Mat_corr_3_Tot(:,:,nd_i-2) =  Mat_corr_3;
                
                %save CORR image
                Fun_Save_Im_Corr(Animal_name, SAVE_Dir_Method1_An, SAVE_Dir_Method2_An, SAVE_Dir_Method3_An,'DAY', SEL_PART, Mat_corr_1, Mat_corr_2, Mat_corr_3,str2num(nday));
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                
%             elseif strcmp(SEL_SEP,'WEEK')
                
                %%%%numero di punti nella sequenza del giorno corrente deveessere uguale a quella dei giorni precedenti (si prende la dimensione più piccola)
                %%% per Method 1 e 2
                diff_Frames = NumFramesPrevDay - NumFrames;
                if diff_Frames > 0
                    FuncAreaValue_Week_Met1 = FuncAreaValue_Week_Met1(1:NumFrames,:);
                    FuncAreaValue_Week_Met2 = FuncAreaValue_Week_Met2(1:NumFrames,:);
                    NumFramesPrevDay = NumFrames;
                elseif diff_Frames < 0
                    FuncAreaValue_Day_Met1  = FuncAreaValue_Day_Met1(1:NumFramesPrevDay,:);
                    FuncAreaValue_Day_Met2  = FuncAreaValue_Day_Met2(1:NumFramesPrevDay,:);
                    NumFramesPrevDay = NumFramesPrevDay;
                end
                %%%%
                FuncAreaValue_Week_Met1 = [FuncAreaValue_Week_Met1 FuncAreaValue_Day_Met1];
                FuncAreaValue_Week_Met2 = [FuncAreaValue_Week_Met2 FuncAreaValue_Day_Met2];
                FuncAreaValue_Week_Met3 = [FuncAreaValue_Week_Met3 FuncAreaValue_Day_Met3];
                
                %last day of the week of interest
                if sum(str2num(nday) == weekly_day_stop)>0
                    
                    nw_ind = nw_ind+1;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% CORRELAZIONE %%%
                    [Mat_corr_1_W Mat_corr_2_W Mat_corr_3_W] = Fun_Compute_Corr(FuncAreaValue_Week_Met1, FuncAreaValue_Week_Met2, FuncAreaValue_Week_Met3, Num_FunAreas);
                    
                    %store CORR matrix
                    Mat_corr_1_W_Tot(:,:,nw_ind) =  Mat_corr_1_W;
                    Mat_corr_2_W_Tot(:,:,nw_ind) =  Mat_corr_2_W;
                    Mat_corr_3_W_Tot(:,:,nw_ind) =  Mat_corr_3_W;
                    
                    %save CORR image
                    Fun_Save_Im_Corr(Animal_name, SAVE_Dir_Method1_W_An, SAVE_Dir_Method2_W_An, SAVE_Dir_Method3_W_An,'WEEK', SEL_PART, Mat_corr_1_W, Mat_corr_2_W, Mat_corr_3_W, nw_ind);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                end                
%             end %end if SEL_SEP save Fig
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end %end if SeqDir
        
    end %end for days
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% SAVE DATA MATRIX CORRELAZIONE %%%
%     if strcmp(SEL_SEP,'DAY')
%         filename_1 = [SAVE_Dir_Method1_An,'Data_Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_1'];
%         save(filename_1,'Mat_corr_1_Tot');
%         
%         filename_1 = [SAVE_Dir_Method2_An,'Data_Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_2'];
%         save(filename_2,'Mat_corr_2_Tot');
%         
%         clear Mat_corr_1_Tot Mat_corr_2_Tot Mat_corr_1 Mat_corr_2
%         
%     elseif strcmp(SEL_SEP,'WEEK')
%         
%         filename_1 = [SAVE_Dir_Method1_An,'Data_Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_1'];
%         save(filename_1,'Mat_corr_1_Tot_W');
%         
%         filename_1 = [SAVE_Dir_Method2_An,'Data_Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_2'];
%         save(filename_2,'Mat_corr_2_Tot_W');
%         
%         clear Mat_corr_1_W_Tot Mat_corr_2_W_Tot Mat_corr_1_W Mat_corr_2_W
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SEL_SEP insieme
    %%% SAVE DATA MATRIX CORRELAZIONE %%%
%     if strcmp(SEL_SEP,'DAY')
        filename_1 = [SAVE_Dir_Method1_An,'Data_Corr_Areas_','DAY','_',SEL_PART,'_Meth_1_',Animal_name];
        save(filename_1,'Mat_corr_1_Tot');
        
        filename_2 = [SAVE_Dir_Method2_An,'Data_Corr_Areas_','DAY','_',SEL_PART,'_Meth_2_',Animal_name];
        save(filename_2,'Mat_corr_2_Tot');
        
        filename_3 = [SAVE_Dir_Method3_An,'Data_Corr_Areas_','DAY','_',SEL_PART,'_Meth_3_',Animal_name];
        save(filename_3,'Mat_corr_3_Tot');
        
        clear Mat_corr_1_Tot Mat_corr_2_Tot Mat_corr_3_Tot Mat_corr_1 Mat_corr_2 Mat_corr_3
        
%     elseif strcmp(SEL_SEP,'WEEK')
        
        filename_1_W = [SAVE_Dir_Method1_W_An,'Data_Corr_Areas_','WEEK','_',SEL_PART,'_Meth_1_',Animal_name];
        save(filename_1_W,'Mat_corr_1_W_Tot');
        
        filename_2_W = [SAVE_Dir_Method2_W_An,'Data_Corr_Areas_','WEEK','_',SEL_PART,'_Meth_2_',Animal_name];
        save(filename_2_W,'Mat_corr_2_W_Tot');
        
        filename_3_W = [SAVE_Dir_Method3_W_An,'Data_Corr_Areas_','WEEK','_',SEL_PART,'_Meth_3_',Animal_name];
        save(filename_3_W,'Mat_corr_3_W_Tot');
        
        clear Mat_corr_1_W_Tot Mat_corr_2_W_Tot Mat_corr_3_W_Tot Mat_corr_1_W Mat_corr_2_W Mat_corr_3_W
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end

display('End Process')

end % end for SEL_PART_LIST 


display('End Complete Process')
toc
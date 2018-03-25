%%% Script to plot the average of the maximum value of fluo on the sequence 

ListAnimalTogether = {      'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
    'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke', 'GCaMPChR2_26_stroke'...
    'GCaMP16_stroke_BoNT','GCaMP18_stroke_BoNT',...
    'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
    'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};


%%%%%%%% User %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UserName = 'CNR-SSSUP';
UsbPort = 'K';
% UserName = 'Stefano';
% UsbPort = 'F';

%%%%%%%% Where SAVE Folder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PER FARLI INSIEME %%%%
Main_SAVE_Dir = [UsbPort,':\LENS\Max_Seq\'];

SEL_PART = {'TASK'};


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
    
    %%%%%%% Info Animal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [weekly_day_stop num_week] = Fun_Info_Animal(Animal_name,list_real_days);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     nw_ind = 0;
    
    NumFramesPrevDay = [];
    
    
    fF = figure('Name', Animal_name);    
    hold on
   
    for nd_i=3:length(NumDaysFolder) %for days
        
        subplot(4,5,nd_i-2)
        
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
%             load([DirDay,'\',SeqDir(3,1).name])
            
            %load Image Sequence -> ImageSequence
            if exist([DirDay,'\','ImageSequence_',Animal_name,'_',nday,'_CenFrame_1_3_ROI_1.mat'])
                load([DirDay,'\','ImageSequence_',Animal_name,'_',nday,'_CenFrame_1_3_ROI_1']);
            else
                error(['ImageSequence_',Animal_name,'_',nday,'_CenFrame_1_3_ROI_1 missing'])
            end
            %load file info Seq -> DataInfoSequence
            if exist([DirDay,'\','DataInfoSequence_',Animal_name,'_',nday,'.mat'])
                load([DirDay,'\','DataInfoSequence_',Animal_name,'_',nday]);
            else
                error(['DataInfoSequence_',Animal_name,'_',nday,' missing'])
            end
            
            
            %%%%%%%%% DAILY TASK SESSION %%%%%%%%%%%%%
            MatrixImageSequence = ImageSequence.MatrixImageSequence;
            NumTrials           = size(MatrixImageSequence,1);
            
            %trials used
            TrialsUsed = find(DataInfoSequence(:,3)==1);
            NumTrialsUsed = length(TrialsUsed);
                
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
                       
            Frame_X_ALL = [];
            Frame_Y_ALL = [];
            
            for int=1:NumTrialsUsed %for trials
                
                Frame_X = [];
                Frame_Y = [];
                
                for inf=3:3:45 %for frames in the trial i-th
                    %                 for inf=1:NumFrames %for frames in the trial i-th
                    
                    Frame_ith = MatrixImageSequence{TrialsUsed(int),:}(:,:,inf);
                    
                    %%roto-translation of the frame i-th
                    Frame_ith = RotoTrans_Image( Frame_ith, rot_transl, nday, rw, cl);
                    [s1 s2 s3] = size(Frame_ith);
                    
                    [maxFrame,ind] = max(Frame_ith(:));
                    [XF,YF] = ind2sub(size(Frame_ith),ind);
                    
                    Frame_X = [Frame_X; XF];
                    Frame_Y = [Frame_Y; YF];
                end
                
                Frame_X_ALL = [Frame_X_ALL Frame_X];
                Frame_Y_ALL = [Frame_Y_ALL Frame_Y];
                
            end %end %for trials
                
            MFX = mean(Frame_X_ALL,2);
            MFY = mean(Frame_Y_ALL,2);
            
%             scatter(MFX,MFY)
                            plot(MFX,MFY,'--o')
            hold on
            
            for tF=1:length(MFX)
                %num
                text(MFX(tF),MFY(tF),num2str(tF))
                
            end
            
            xlim([0 512])
            ylim([0 512])           
                            
        end
    end
    set(gca,'Ydir','reverse')
    saveas(gca,[Animal_name],'fig')
    saveas(gca,[Animal_name],'jpeg')
    close
    
    
end

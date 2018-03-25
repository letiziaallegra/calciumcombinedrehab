%script per plottare andamento della correlazione cotruita con i due METHOD
%nel caso WEEK in REST e TASK rispetto alla M1

clear
close all
clc

HD = 'F';
User = 'Stefano';

% HD = 'I';
% User = 'CNR-SSSUP';

PART        = {'TASK' , 'REST'};
TASK_FOLDER = [HD,':\LENS\Correlazioni_SelSeq\Corr_Areas_SEPWEEK_PART',PART{1}];
REST_FOLDER = [HD,':\LENS\Correlazioni_SelSeq\Corr_Areas_SEPWEEK_PART',PART{2}];
% TASK_FOLDER = [HD,':\LENS\Correlazioni_Parziali_SelSeq\Partial_Corr_Areas_SEPWEEK_PART',PART{1}];
% REST_FOLDER = [HD,':\LENS\Correlazioni_Parziali_SelSeq\Partial_Corr_Areas_SEPWEEK_PART',PART{2}];

METHOD = {'Method1', 'Method2', 'Method3'};
R2 = 0; %r^2 yes
ExtraLabel = '_soloChr2_NoChR2_16rehab';

FOLDER = {TASK_FOLDER REST_FOLDER};


%%%
%%% rows:      TASK ; REST
%%% coloumns:  contr; stroke; rehab1w; rehab4w;
%%% inside col: [col1 col2 ..] animals
MATRIX_INFO_M1 = cell(2,4);
MATRIX_INFO_M2 = cell(2,4);
MATRIX_INFO_M3 = cell(2,4);
indTreat = [1 2 3 4];

for i_f=1:length(FOLDER) %for FOLDER
    
    
    
    for i_m=1:length(METHOD) %for METHOD
        
        DataStoreBuf = cell(1,4);
        i_cnt    = 0;
        i_rehab1 = 0;
        i_rehab4 = 0;
        i_stk    = 0;
        
        F = [FOLDER{i_f},'\',METHOD{i_m}];
        file_F = dir(F);
        
        for i_fil=3:length(file_F) %for files
            
            AnimalName = file_F(i_fil,1).name;
            
%             if strfind(AnimalName,'GCaMP') &&  ~strcmp(AnimalName,'GCaMP16_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMP17_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMP18_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMPChR2_16_stroke_BoNT') %if GCAMP
%             if strfind(AnimalName,'GCaMP') &&  ~strcmp(AnimalName,'GCaMP16_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMP17_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMP18_stroke_BoNT')%if GCAMP
            if strfind(AnimalName,'GCaMP') &&  ~strcmp(AnimalName,'GCaMPChR2_16_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMPChR2_22_stroke')%if GCAMP
%             if strfind(AnimalName,'GCaMP')
                
                DataFileName = ['Data_Corr_Areas_WEEK_',PART{i_f},'_',[METHOD{i_m}([1:4]),'_', METHOD{i_m}([7])],'_',AnimalName];
                DataCurrDir  = [F,'\',AnimalName];
                
                %load file
                load([DataCurrDir,'\',DataFileName]) %-> Mat_corr_1_W_tot or Mat_corr_2_W_tot
                
                if i_m == 1
                    
                    MData = Mat_corr_1_W_Tot;
                elseif i_m == 2
                    
                    MData = Mat_corr_2_W_Tot;
                elseif i_m == 3
                    
                    MData = Mat_corr_3_W_Tot;
                end
                
                if strfind(AnimalName,'control') %if control
                                       
                    %control
                    i_cnt = i_cnt +1;
                    DataStoreBuf{1,indTreat(1)}(:,:,i_cnt) = MData(:,:);                       
                
                elseif strfind(AnimalName,'stroke') %if stroke and BoNT
                    
                    if strfind(AnimalName,'BoNT') %if BoNT
                        
                        %1w
                        i_rehab1 = i_rehab1 +1;
                        DataStoreBuf{1,indTreat(3)}(:,:,i_rehab1) = squeeze(MData(:,:,1));
                        %4w
                        i_rehab4 = i_rehab4 +1;
                        DataStoreBuf{1,indTreat(4)}(:,:,i_rehab4) = squeeze(MData(:,:,4));
                                                
                    else %else if stroke
                        
                        %stroke
                        i_stk = i_stk +1;
                        DataStoreBuf{1,indTreat(2)}(:,:,i_stk)    = MData(:,:);
                        
                    end %end if stroke and BoNT
                
                end %endif control
                
            end %end if GCAMP        
            
            
        end%end for files    
        
        if i_m == 1
            
            MATRIX_INFO_M1(i_f,:) = DataStoreBuf;
        elseif i_m == 2
            
            MATRIX_INFO_M2(i_f,:) = DataStoreBuf;
        elseif i_m == 3
            
            MATRIX_INFO_M3(i_f,:) = DataStoreBuf;
        end
        
        
        
    end %end for METHOD

    
end %end for FOLDER
    

%%%
AREA_NAME   = {'M2','CFA','M1Sh','M1HL','M1Tk','S1','S1Bc','RSD','LPTa','V1','AuD'};
TREAT_NAME = {'C_T' 'S_T' 'R1w_T' 'R4w_T' 'C_R' 'S_R' 'R1w_R' 'R4w_R'};

%%% plot %%%

COLOR_T =  colormap(hsv(4));
COLOR_T(4,:) = [0 0 0];


for i_M=1:3 %METHOD
    
    %method
    switch i_M
        case 1
            MATRIX_INFO = MATRIX_INFO_M1;
            figure('Name','Average Correlation with M1 area among groups - Method 1')
            INC = 0;
        case 2
            MATRIX_INFO = MATRIX_INFO_M2;
            figure('Name','Average Correlation with M1 area among groups - Method 2')
            INC = 0;
        case 3
            MATRIX_INFO = MATRIX_INFO_M3;
            figure('Name','Average Correlation with M1 area among groups - Method 3')
            INC = 0;
    end

    
    for i_P=1:size(MATRIX_INFO,1) %PART
        
        for i_T=1:size(MATRIX_INFO,2) %TREAT
            
            %R^2
            if R2==1
                MATRIX_INFO{i_P,i_T} = MATRIX_INFO{i_P,i_T}.^2;                
            end
            
            %
            Med1 = nanmean(MATRIX_INFO{i_P,i_T},3);
            SDE1 =  nanstd(MATRIX_INFO{i_P,i_T},[],3)/sqrt(size(MATRIX_INFO{i_P,i_T},3));
            
            subplot(2,size(MATRIX_INFO,2),i_T+INC)
            imagesc(Med1)
                    colorbar
            title([PART{i_P},'',num2str(i_T)])
            set(gca,'xTick',     1:length(AREA_NAME))
            set(gca,'xTickLabel',AREA_NAME)
            set(gca,'yTick',     1:length(AREA_NAME))
            set(gca,'yTickLabel',AREA_NAME)
                  
        end
%         subplot(2,length(METHOD),i_P)
        rotateXLabels( gca, 60)
%         subplot(2,length(METHOD),i_P+2)
%         rotateXLabels( gca, 60)
        INC = 4;
        
        
    end
end
subplot(2,length(METHOD),1)
legend(TREAT_NAME(1:4),'Location','Best')
subplot(2,length(METHOD),4)
legend(TREAT_NAME(5:8),'Location','Best')

%%% save
if R2 == 1
    R_name = '_R2';
else
    R_name = '_r';
end
filename_1 = ['Average_Corr_with_CFA_among_groups',R_name,'_FIG'];
% saveas(gca,[filename_1,ExtraLabel],'fig');
% saveas(gca,[filename_1,ExtraLabel],'jpeg');





%%% plot %%%

COLOR_T =  colormap(hsv(4));
COLOR_T(4,:) = [0 0 0];


for i_M=1:3 %METHOD
    
    %method
    switch i_M
        case 1
            MATRIX_INFO = MATRIX_INFO_M1;
            figure('Name','Median Correlation with M1 area among groups - Method 1')
            INC = 0;
        case 2
            MATRIX_INFO = MATRIX_INFO_M2;
            figure('Name','Median Correlation with M1 area among groups - Method 2')
            INC = 0;
        case 3
            MATRIX_INFO = MATRIX_INFO_M3;
            figure('Name','Median Correlation with M1 area among groups - Method 3')
            INC = 0;
    end

    
    for i_P=1:size(MATRIX_INFO,1) %PART
        
        for i_T=1:size(MATRIX_INFO,2) %TREAT
            
            %R^2
            if R2==1
                MATRIX_INFO{i_P,i_T} = MATRIX_INFO{i_P,i_T}.^2;                
            end
            
            %
            Med1 = nanmedian(MATRIX_INFO{i_P,i_T},3);
            SDE1 =  nanstd(MATRIX_INFO{i_P,i_T},[],3)/sqrt(size(MATRIX_INFO{i_P,i_T},3));
            
            subplot(2,size(MATRIX_INFO,2),i_T+INC)
            imagesc(Med1)
%             imagesc(SDE1)
            %         colorbar
            title([PART{i_P},'',num2str(i_T)])
            set(gca,'xTick',     1:length(AREA_NAME))
            set(gca,'xTickLabel',AREA_NAME)
            set(gca,'yTick',     1:length(AREA_NAME))
            set(gca,'yTickLabel',AREA_NAME)
                  
        end
%         subplot(2,length(METHOD),i_P)
        rotateXLabels( gca, 60)
%         subplot(2,length(METHOD),i_P+2)
%         rotateXLabels( gca, 60)
        INC = 4;
        
        
    end
end
subplot(2,length(METHOD),1)
legend(TREAT_NAME(1:4),'Location','Best')
subplot(2,length(METHOD),4)
legend(TREAT_NAME(5:8),'Location','Best')

%%% save
if R2 == 1
    R_name = '_R2';
else
    R_name = '_r';
end
filename_1 = ['Median_Corr_with_CFA_among_groups',R_name,'_FIG'];
% saveas(gca,[filename_1,ExtraLabel],'fig');
% saveas(gca,[filename_1,ExtraLabel],'jpeg');
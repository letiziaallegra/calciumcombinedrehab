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

METHOD = {'Method1', 'Method2', 'Method3'};
AreaOfInterest = 2; %CFA
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
        
        F = [FOLDER{i_f},'\',METHOD{i_m}];
        file_F = dir(F);
        
        for i_fil=3:length(file_F) %for files
            
            AnimalName = file_F(i_fil,1).name;
            
%             if strfind(AnimalName,'GCaMP') &&  ~strcmp(AnimalName,'GCaMP16_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMP17_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMP18_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMPChR2_16_stroke_BoNT') %if GCAMP
%             if strfind(AnimalName,'GCaMP') &&  ~strcmp(AnimalName,'GCaMP16_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMP17_stroke_BoNT') &&  ~strcmp(AnimalName,'GCaMP18_stroke_BoNT')%if GCAMP
            if strfind(AnimalName,'GCaMP') &&  ~strcmp(AnimalName,'GCaMPChR2_16_stroke_BoNT') %if GCAMP
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
                    DataStoreBuf{1,indTreat(1)} = [DataStoreBuf{1,indTreat(1)} MData(:,AreaOfInterest)];                       
                
                elseif strfind(AnimalName,'stroke') %if stroke and BoNT
                    
                    if strfind(AnimalName,'BoNT') %if BoNT
                        
                        %1w
                        DataStoreBuf{1,indTreat(3)} = [DataStoreBuf{1,indTreat(3)} squeeze(MData(:,AreaOfInterest,1))];
                        %4w
                        DataStoreBuf{1,indTreat(4)} = [DataStoreBuf{1,indTreat(4)} squeeze(MData(:,AreaOfInterest,4))];
                        
                        
                    else %else if stroke
                        
                        %stroke
                        DataStoreBuf{1,indTreat(2)} = [DataStoreBuf{1,indTreat(2)} MData(:,AreaOfInterest)];
                        
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
figure('Name','Average Correlation with M1 area among groups')
COLOR_T =  colormap(hsv(4));
COLOR_T(4,:) = [0 0 0];
INC = 0;

for i_P=1:size(MATRIX_INFO_M1,1) %PART   
    
    for i_T=1:size(MATRIX_INFO_M1,2) %TREAT
        
        
        %R^2
        if R2==1
            MATRIX_INFO_M1{i_P,i_T} = MATRIX_INFO_M1{i_P,i_T}.^2;
            MATRIX_INFO_M2{i_P,i_T} = MATRIX_INFO_M2{i_P,i_T}.^2;
            MATRIX_INFO_M3{i_P,i_T} = MATRIX_INFO_M3{i_P,i_T}.^2;
        end
        
        Med1 = nanmean(MATRIX_INFO_M1{i_P,i_T},2);
        SDE1 = nanstd(MATRIX_INFO_M1{i_P,i_T},[],2)/sqrt(size(MATRIX_INFO_M1{i_P,i_T},2));
        Med1(AreaOfInterest) = NaN;
        subplot(2,length(METHOD),i_P+INC)
        hold on
        errorbar(Med1,SDE1,'Color',COLOR_T(i_T,:))
        title(['Method 1 ',PART{i_P}])
        xlim([0.5 11.5])
        set(gca,'xTick',     1:length(AREA_NAME))
        set(gca,'xTickLabel',AREA_NAME)
        ylim([0 1.1])
        set(gca,'yTick',     0.1:0.1:1)        
        
        Med2 = nanmean(MATRIX_INFO_M2{i_P,i_T},2);
        SDE2 = nanstd(MATRIX_INFO_M2{i_P,i_T},[],2)/sqrt(size(MATRIX_INFO_M2{i_P,i_T},2));
        Med2(AreaOfInterest) = NaN;
        subplot(2,length(METHOD),i_P+1+INC)
        hold on
        errorbar(Med2,SDE2,'color',COLOR_T(i_T,:)) 
        xlim([0.5 11.5])
        title(['Method 2 ',PART{i_P}])
        set(gca,'xTick',     1:length(AREA_NAME))
        set(gca,'xTickLabel',AREA_NAME)
        ylim([-0.2 0.8])
        set(gca,'yTick',     -0.2:0.1:0.8)  
        
        Med3 = nanmean(MATRIX_INFO_M3{i_P,i_T},2);
        SDE3 = nanstd(MATRIX_INFO_M3{i_P,i_T},[],2)/sqrt(size(MATRIX_INFO_M3{i_P,i_T},2));
        Med3(AreaOfInterest) = NaN;
        subplot(2,length(METHOD),i_P+2+INC)
        hold on
        plot(1:11,Med3,'color',COLOR_T(i_T,:),'Marker','o')
%         errorbar(Med3,SDE3,'color',COLOR_T(i_T,:),'Marker','o')
        xlim([0.5 11.5])
        title(['Method 3 ',PART{i_P}])
        set(gca,'xTick',     1:length(AREA_NAME))
        set(gca,'xTickLabel',AREA_NAME)
        ylim([-0.2 1.2])
        set(gca,'yTick',     -0.2:0.1:0.8)
        
    end
    subplot(2,length(METHOD),i_P)
    rotateXLabels( gca, 60)
    subplot(2,length(METHOD),i_P+2)
    rotateXLabels( gca, 60)
    INC = 2;
    
    
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
figure('Name','Median Correlation with M1 area among groups')
COLOR_T =  colormap(hsv(8));
COLOR_T(4,:) = [0 0 0];
INC = 0;

for i_P=1:size(MATRIX_INFO_M1,1) %PART   
    
    for i_T=1:size(MATRIX_INFO_M1,2) %TREAT
        
        %R^2
        if R2==1
            MATRIX_INFO_M1{i_P,i_T} = MATRIX_INFO_M1{i_P,i_T}.^2;
            MATRIX_INFO_M2{i_P,i_T} = MATRIX_INFO_M2{i_P,i_T}.^2;
            MATRIX_INFO_M3{i_P,i_T} = MATRIX_INFO_M3{i_P,i_T}.^2;
        end
        
        
        Med1 = nanmedian(MATRIX_INFO_M1{i_P,i_T},2);
        SDE1 = nanstd(MATRIX_INFO_M1{i_P,i_T},[],2)/sqrt(size(MATRIX_INFO_M1{i_P,i_T},2));
        Med1(AreaOfInterest) = NaN;
        subplot(2,length(METHOD),i_P+INC)
        hold on
        errorbar(Med1,SDE1,'Color',COLOR_T(i_T,:))
        title(['Method 1 ',PART{i_P}])
        xlim([0.5 11.5])
        set(gca,'xTick',     1:length(AREA_NAME))
        set(gca,'xTickLabel',AREA_NAME)
        ylim([0 1.1])
        set(gca,'yTick',     0.1:0.1:1)  
        
        
        Med2 = nanmedian(MATRIX_INFO_M2{i_P,i_T},2);
        SDE2 = nanstd(MATRIX_INFO_M2{i_P,i_T},[],2)/sqrt(size(MATRIX_INFO_M2{i_P,i_T},2));
        Med2(AreaOfInterest) = NaN;
        subplot(2,length(METHOD),i_P+1+INC)
        hold on
        errorbar(Med2,SDE2,'color',COLOR_T(i_T,:)) 
        xlim([0.5 11.5])
        title(['Method 2 ',PART{i_P}])
        set(gca,'xTick',     1:length(AREA_NAME))
        set(gca,'xTickLabel',AREA_NAME)
        ylim([-0.2 0.8])
        set(gca,'yTick',     -0.2:0.1:0.8) 
        
        Med3 = nanmedian(MATRIX_INFO_M3{i_P,i_T},2);
        SDE3 = nanstd(MATRIX_INFO_M3{i_P,i_T},[],2)/sqrt(size(MATRIX_INFO_M3{i_P,i_T},2));
        Med3(AreaOfInterest) = NaN;
        subplot(2,length(METHOD),i_P+2+INC)
        hold on
        plot(1:11,Med3,'color',COLOR_T(i_T,:),'Marker','o')
%         errorbar(Med3,SDE3,'color',COLOR_T(i_T,:))
        xlim([0.5 11.5])
        title(['Method 3 ',PART{i_P}])
        set(gca,'xTick',     1:length(AREA_NAME))
        set(gca,'xTickLabel',AREA_NAME)
        ylim([-0.2 1.2])
        set(gca,'yTick',     -0.2:0.1:1.2)
    
    end
    subplot(2,length(METHOD),i_P+INC)
    rotateXLabels( gca, 60)
    subplot(2,length(METHOD),i_P+1+INC)
    rotateXLabels( gca, 60)
    subplot(2,length(METHOD),i_P+2+INC)
    rotateXLabels( gca, 60)

    INC = 2;
    
    
end
% subplot(2,length(METHOD),1)
% legend(TREAT_NAME(1:4),'Location','Best')
% subplot(2,length(METHOD),2)
% legend(TREAT_NAME(5:8),'Location','Best')
% subplot(2,length(METHOD),3)
% legend(TREAT_NAME(1:4),'Location','Best')
% subplot(2,length(METHOD),4)
% legend(TREAT_NAME(5:8),'Location','Best')

subplot(2,length(METHOD),2)
legend(TREAT_NAME(1:4),'Location','Best')
subplot(2,length(METHOD),5)
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








    
%script per plottare andamento della correlazione cotruita con i due METHOD
%nel caso WEEK in REST e TASK rispetto alla M1

clear
% close all
clc

HD = 'F';
User = 'Stefano';

% HD = 'I';
% User = 'CNR-SSSUP';

PART        = {'TASK'};
% TASK_FOLDER = [HD,':\LENS\Correlazioni_Parziali_SelSeq\Partial_Corr_Areas_SEPWEEK_PART',PART{1}];
TASK_FOLDER = [HD,':\LENS\Correlazioni_SelSeq\Corr_Areas_SEPWEEK_PART',PART{1}];
TASK_FOLDER_SAVE = [HD,':\LENS\Correlazioni_Parziali_SelSeq_AssiUguale'];

METHOD = {'Method1'};
AreaOfInterest = 2; %CFA
R2 = 0; %r^2 yes
ExtraLabel = '_soloChr2_NoChR2_16rehab';

FOLDER = {TASK_FOLDER};


%%%
%%% rows:      TASK ; REST
%%% coloumns:  contr; stroke; rehab1w; rehab4w;
%%% inside col: [col1 col2 ..] animals
MATRIX_INFO_M1 = cell(2,4);
MATRIX_INFO_M2 = cell(2,4);
MATRIX_INFO_M3 = cell(2,4);
indTreat = [1 2 3 4];


Rehab_Mat = zeros(11,11,6)*NaN;
i_Rehab = 0;
Control_Mat = zeros(11,11,4)*NaN;
i_Control = 0;
Stroke_Mat = zeros(11,11,6)*NaN;
i_Stroke = 0;

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
                M = Mat_corr_1_W_Tot;
                
                if size(M,3)>1
                    
                    for is=1:size(M,3)
                        
                        figure
                        imagesc(M(:,:,is))
                        caxis([0 1])
                        colorbar
%                         saveas(gca,[TASK_FOLDER_SAVE,'\',AnimalName,'_','WEEK_',num2str(is)],'jpeg')
                        close
                    end
                    if is==4
                        i_Rehab = i_Rehab+1;
                        Rehab_Mat(:,:,i_Rehab) = M(:,:,is).^2;
%                         Rehab_Mat(:,:,i_Rehab) = abs(Mat_corr_1_W_Tot(:,:,is));
                    end
                    
                else
                    figure
                    imagesc(M)
                    caxis([0 1])
                    colorbar
%                     saveas(gca,[TASK_FOLDER_SAVE,'\',AnimalName,'_','WEEK_4'],'jpeg')
                    close
                    if strfind(AnimalName,'control')
                        i_Control= i_Control+1;
                        Control_Mat(:,:,i_Control) = M.^2;
%                         Control_Mat(:,:,i_Control) = abs(Mat_corr_1_W_Tot);
                    else
                        i_Stroke= i_Stroke+1;
                        Stroke_Mat(:,:,i_Stroke) = M.^2;
%                         Stroke_Mat(:,:,i_Stroke) = abs(Mat_corr_1_W_Tot);
                    end
                    
                    
                end
                
            end
        end
    end
end


%
figure
maxC = 1.0;
minC = 0;
subplot(131)
imagesc(nanmean(Control_Mat,3))
caxis([minC maxC])
subplot(132)
imagesc(nanmean(Stroke_Mat,3))
caxis([minC maxC])
subplot(133)
imagesc(nanmean(Rehab_Mat,3))
caxis([minC maxC])



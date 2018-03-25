% Sequence of frames with a central one
% TrOK = [1:26];
% TrOK = [13];
% posPlot = [5 10 25 9];

%%% Initial Check %%%
if ~exist('ImageSequence')
    error('load ImageSequence')
end
%%%

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


for inf=1:NumFrames %-> num frames
    
    %to store the same frame of different trials
    Frame_Trials_Matrix = zeros(rw,cl,NumTrialsUsed);
    for int=1:NumTrialsUsed
        
        int_used = TrialsUsed(int);        
        Frame_Trials_Matrix(:,:,int) = MatrixImageSequence{int_used,1}(:,:,inf)/100;        
    end
    
    %median of all of the i-th frame
    Med = nanmedian(Frame_Trials_Matrix,3);
    Med = medfilt2(Med);
    
    Med_Frames_Matrix(:,:,inf) = Med; 
    
    Threshold(inf) = 3*std2(Med);    
    
    %to store median/mean -> threshold
    Med_Th = im2bw(Med,Threshold(inf));
    
    se = strel('disk',2);
    Med_Th_Dil = imerode(imdilate(Med_Th,se),se);
    
    %filtering
    Med_Th_Dil = medfilt2(Med_Th_Dil,[15 15]);
    
%     %figure
%     figure
% %     imagesc(Med_Th_Dil)
%     figure
%     ScriptExtractFluoSequence_Plot_update      
    
    [bound,label] = bwboundaries(Med_Th_Dil);
    temp1 = regionprops(label,'Centroid');
    
    %%% Max Area  %%% %%% %%% %%% %%% %%% %%%
    [n in] = max(cellfun('length',bound));
    
    %Store Cetroid of the Max
    if ~isempty(in)
        Cent_X_Y_Max = round(temp1(in,1).Centroid);          
    else    
        Cent_X_Y_Max = [NaN NaN];    
    end    
    Cent_X_Y_Max_Store(inf,:) = Cent_X_Y_Max; 
   
    %%% Size max areas  %%% %%% %%% %%% %%% %%%
    MaxSizeROI = 300;
    [inAll] = find((cellfun('length',bound)>MaxSizeROI)==1);
    
    %Store Cetroid of the all
    if ~isempty(inAll)
        for iAl=1:length(inAll)
            Cent_X_Y_All_Store_Buf(iAl,:) = round(temp1(inAll(iAl),1).Centroid);           
        end
    else    
        Cent_X_Y_All_Store_Buf = [NaN NaN];    
    end  
    Cent_X_Y_Max_All_Store{inf,:} = Cent_X_Y_All_Store_Buf; 
     %%% %%%
        
     
    
    Med_Frames_Matrix_Th(:,:,inf) = Med_Th_Dil;

end

save('Cent_X_Y_Max_Store','Cent_X_Y_Max_Store')
save('Cent_X_Y_Max_All_Store','Cent_X_Y_Max_All_Store')

    
%PLOT%
figure
sc = 6;
sr = ceil(NumFrames/sc);
for inf=1:NumFrames
    
    Im = squeeze(Med_Frames_Matrix(:,:,inf));
        
    subplot(sr,sc,inf)
        
    imagesc(Im);
    colormap hot
    caxis([-0.2 0.2])
    
    if inf == round(ImageSequence.CentFrame/ImageSequence.downsamplingfactor)
        title('Cent Frame')
    end
    
end
colorbar


% %PLOT after Threshold%
% figure
% sc = 6;
% sr = ceil(NumFrames/sc);
% for inf=1:NumFrames
%     
%     Im = squeeze(Med_Frames_Matrix_Th(:,:,inf));
%     
% 
%     
%     subplot(sr,sc,inf)
%     
%     imagesc(Im);
%     caxis([0 1])
%     colormap gray
%     
%     if inf == round(ImageSequence.CentFrame/ImageSequence.downsamplingfactor)
%         title('Cent Frame')
%     end
%     
% end
% colorbar


% %PLOT first 12%
% figure
% for inf=1:12
%     
%     Im = squeeze(Med_Frames_Matrix_Th(:,:,inf));
%     
% 
%     
%     subplot(2,6,inf)
%     
%     imagesc(Im);
%     caxis([0 1])
%     colormap gray
%     
%     if inf == round(ImageSequence.CentFrame/ImageSequence.downsamplingfactor)
%         title('Cent Frame')
%     end
%     
% end
% colorbar

%PLOT after Threshold MEDIAN%
figure('Name','Median Centroid')
sc = 6;
sr = ceil(NumFrames/sc);
for inf=1:NumFrames
    
    Im = squeeze(Med_Frames_Matrix_Th(:,:,inf));
        
    subplot(sr,sc,inf)
    
    imagesc(Im);
    caxis([0 1])
    colormap gray
    
    if inf == round(ImageSequence.CentFrame/ImageSequence.downsamplingfactor)
        title('Cent Frame')
    end
    
    [CC] = nanmedian(Cent_X_Y_Max_Store,1)
    hold on;
    plot(CC(1),CC(2),'r.','MarkerSize',10) 
    
end
colorbar

%PLOT after Threshold%
figure('Name','Single Centroid for each frame')
sc = 6;
sr = ceil(NumFrames/sc);
for inf=1:NumFrames
    
    Im = squeeze(Med_Frames_Matrix_Th(:,:,inf));
        
    subplot(sr,sc,inf)
    
    imagesc(Im);
    caxis([0 1])
    colormap gray
    
    if inf == round(ImageSequence.CentFrame/ImageSequence.downsamplingfactor)
        title('Cent Frame')
    end
    
    [CC] = nanmedian(Cent_X_Y_Max_Store(inf,:),1);
    hold on;
    plot(CC(1),CC(2),'r.','MarkerSize',10); 
    
end
colorbar


%PLOT max costruito con media pesata%
figure('Name','Weight Mean Centroid for each frame')
X_Y_Weighted_Centroid = [];
sc = 6;
sr = ceil(NumFrames/sc);
for inf=1:NumFrames
    
    Med_       = squeeze(Med_Frames_Matrix(:,:,inf));
    Med_Th_Dil = squeeze(Med_Frames_Matrix_Th(:,:,inf));
    
    %find white pixels
    [wp_index] = find(Med_Th_Dil==1);
    [wp_X wp_Y] = ind2sub(size(Med_Th_Dil),wp_index);
    
    %find gray values corresponding to white pixels
    numTone    = Med(wp_index);
    
    %max tones
    maxNumTone_i     = find(numTone == max(numTone));
    maxNumTone_index = wp_index(maxNumTone_i);
    [wp_X_max wp_Y_max] = ind2sub(size(Med_Th_Dil),maxNumTone_index);
    
    if ~isempty(wp_X_max)
    %remap intervals into 0-100 range
    max_numTone = max(numTone);
    min_numTone = min(numTone);
    m_t = (100-0)/(max_numTone(1)-min_numTone(1));
    q_t = 100-m_t*(max_numTone);
    numTone_R = m_t.*numTone + q_t;
        
    %multiplipy for the numTone
    wp_X_Tone = wp_X.*numTone_R;
    wp_Y_Tone = wp_Y.*numTone_R;
    
    %find centroid coordinates (weighted on the gray tone)
    x = round(sum(wp_X_Tone)/sum(numTone_R));
    y = round(sum(wp_Y_Tone)/sum(numTone_R));
    else
        x = NaN;
        y = NaN;
        wp_index = 1;
        numTone_R = 0;
    end
    
    CM = zeros(512,512);
    CM(wp_index) = numTone_R;  
        
    subplot(sr,sc,inf)
    imagesc(CM)
    hold on
    text(x,y, 'X','Color','red','FontSize',10);
    
    if inf == round(ImageSequence.CentFrame/ImageSequence.downsamplingfactor)
        title('Cent Frame')
    end
   
    X_Y_Weighted_Centroid = [X_Y_Weighted_Centroid;[x y]];
    
    
end
colorbar

    

% Sequence of frames with a central one
% TrOK = [1:26];
% TrOK = [13];
% posPlot = [5 10 25 9];
clear all
close all
clc

User = 'CNR-SSSUP';

%%%%%%%%%%%
Animal_name = 'GCaMPChR2_3_stroke_BoNT';
%%%%%%%%%%%

MainDir = ['I:\LENS\_data_MAT_GCamp\',Animal_name,'\'];
% MainDir = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_data_MAT_GCamp\',Animal_name,'\'];
NumDaysFolder = dir(MainDir);

X_Y_Weighted_Centroid_DAY     = [];
X_Y_Weighted_Centroid_std_DAY = [];
X_Y_Weighted_Centroid_Other_DAY = [];
X_Y_Weighted_Centroid_std_Other_DAY = [];

for nd_i=3:length(NumDaysFolder)
    
    nday = NumDaysFolder(nd_i,1).name;
    
    DirDay = [MainDir,nday,'\','Sequence_Audio'];
    
    SeqDir = dir(DirDay);
    
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
        
        int_n = 1;
        for inf=1:NumFrames
            %to store the same frame of different trials
            Frame_Trials_Matrix = zeros(rw,cl,NumTrialsUsed);
            for int=1:NumTrialsUsed                
                int_used = TrialsUsed(int);
                Frame_Trials_Matrix(:,:,int) = MatrixImageSequence{int_used,1}(:,:,inf)/100;
            end
            %median of all of the i-th frame
            Med = nanmedian(Frame_Trials_Matrix,3);
            Med = medfilt2(Med);
            
            ThreshMatrix(:,:,int_n) = Med; 
            int_n = int_n+1;
        end
        Threshold_tot = median(median(median((ThreshMatrix)))) + 2*std2(ThreshMatrix);
        
        
        inf_i = 1;
        inf_end = NumFrames;
        %gcamp3
%         if str2num(nday) == 5
%             inf_i = 5;
%             inf_end = 13;      
%         end
        
        
        for inf=inf_i:inf_end %-> num frames
            
            %to store the same frame of different trials
            Frame_Trials_Matrix = zeros(rw,cl,NumTrialsUsed);
            for int=1:NumTrialsUsed
                
                int_used = TrialsUsed(int);
                Frame_Trials_Matrix(:,:,int) = MatrixImageSequence{int_used,1}(:,:,inf)/100;
            end
            
            %median of all of the i-th frame
            Med = nanmedian(Frame_Trials_Matrix,3);
            Med = medfilt2(Med);
            
%             checkPercentage = 0;
%             Thresh_num = 4;
%             while checkPercentage == 0 && Thresh_num>1.5
%                 
%                 Threshold(inf) = Thresh_num*std2(Med);
%                 
%                 %to store median/mean -> threshold
%                 Med_Th = im2bw(Med,Threshold(inf));
%                 
%                 se = strel('disk',2);
%                 Med_Th_Dil = imerode(imdilate(Med_Th,se),se);
%                 
%                 %filtering
%                 Med_Th_Dil = medfilt2(Med_Th_Dil,[15 15]);
%                 
%                 checkPercentage =  (sum(sum(Med_Th_Dil))/(size(Med_Th_Dil,1)*size(Med_Th_Dil,2))*100) >= 1;
%                 Thresh_num = (Thresh_num-0.5);
%                 
%             end
%             Threshold(inf) = 3*std2(Med);
            
            
            Threshold(inf) = Threshold_tot  ;
            %to store median/mean -> threshold
            Med_Th = im2bw(Med,Threshold(inf));
            
            se = strel('disk',2);
            Med_Th_Dil = imerode(imdilate(Med_Th,se),se);
            
            %filtering
            Med_Th_Dil = medfilt2(Med_Th_Dil,[15 15]);
            
            
%             %     %figure
%             %     figure
%             % %     imagesc(Med_Th_Dil)
%             %     figure
%             %     ScriptExtractFluoSequence_Plot_update
%             
%             [bound,label] = bwboundaries(Med_Th_Dil);
%             temp1 = regionprops(label,'Centroid');
%             
%             %%% Max Area  %%% %%% %%% %%% %%% %%% %%%
%             [n in] = max(cellfun('length',bound));
%             
%             %Store Cetroid of the Max
%             if ~isempty(in)
%                 Cent_X_Y_Max = round(temp1(in,1).Centroid);
%             else
%                 Cent_X_Y_Max = [NaN NaN];
%             end
%             Cent_X_Y_Max_Store(inf,:) = Cent_X_Y_Max;
%             
%             %%% Size max areas  %%% %%% %%% %%% %%% %%%
%             MaxSizeROI = 300;
%             [inAll] = find((cellfun('length',bound)>MaxSizeROI)==1);
%             
%             %Store Cetroid of the all
%             if ~isempty(inAll)
%                 for iAl=1:length(inAll)
%                     Cent_X_Y_All_Store_Buf(iAl,:) = round(temp1(inAll(iAl),1).Centroid);
%                 end
%             else
%                 Cent_X_Y_All_Store_Buf = [NaN NaN];
%             end
%             Cent_X_Y_Max_All_Store{inf,:} = Cent_X_Y_All_Store_Buf;
%             %%% %%%
            
            
            %store median images per quel giorno
            Med_Frames_Matrix(:,:,inf)    = Med;
            Med_Frames_Matrix_Th(:,:,inf) = Med_Th_Dil;
            
        end
        
        
        
    else
        %seq non ok
        SeqDirOKNO = 0;
    end
    
%     save('Cent_X_Y_Max_Store','Cent_X_Y_Max_Store')
%     save('Cent_X_Y_Max_All_Store','Cent_X_Y_Max_All_Store')
    
    
    %PLOT%
    H_Seq_Fig = figure('Name',['Image_Sequence day_',nday]);
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
    H_Seq_Fig_filename = [Animal_name,'_Seq_Fig_',nday];
    saveas(H_Seq_Fig,H_Seq_Fig_filename,'fig');
    
    
    
    
    % NO %
    % %PLOT after Threshold% NO %
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
    
    
    % NO %
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
    
    % %PLOT after Threshold MEDIAN%
    % figure('Name','Median Centroid')
    % sc = 6;
    % sr = ceil(NumFrames/sc);
    % for inf=1:NumFrames
    %
    %     Im = squeeze(Med_Frames_Matrix_Th(:,:,inf));
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
    %     [CC] = nanmedian(Cent_X_Y_Max_Store,1)
    %     hold on;
    %     plot(CC(1),CC(2),'r.','MarkerSize',10)
    %
    % end
    % colorbar
    
%     %PLOT after Threshold%
%     figure('Name','Single Centroid for each frame')
%     sc = 6;
%     sr = ceil(NumFrames/sc);
%     for inf=1:NumFrames
%     
%         Im = squeeze(Med_Frames_Matrix_Th(:,:,inf));
%     
%         subplot(sr,sc,inf)
%     
%         imagesc(Im);
%         caxis([0 1])
%         colormap gray
%     
%         if inf == round(ImageSequence.CentFrame/ImageSequence.downsamplingfactor)
%             title('Cent Frame')
%         end
%     
%         [CC] = nanmedian(Cent_X_Y_Max_Store(inf,:),1);
%         hold on;
%         plot(CC(1),CC(2),'r.','MarkerSize',10);
%     
%     end
%     colorbar
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if 0
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %PLOT max costruito con media pesata%
%         % figure('Name','Weight Mean Centroid for each frame')
%         X_Y_Weighted_Centroid     = [];
%         X_Y_Weighted_Centroid_std = [];
%         sc = 6;
%         sr = ceil(NumFrames/sc);
%         for inf=1:NumFrames
%             
%             Med       = squeeze(Med_Frames_Matrix(:,:,inf));
%             Med_Th_Dil = squeeze(Med_Frames_Matrix_Th(:,:,inf));
%             
%             %find white pixels
%             [wp_index] = find(Med_Th_Dil==1);
%             [wp_X wp_Y] = ind2sub(size(Med_Th_Dil),wp_index);
%             
%             %find gray values corresponding to white pixels
%             numTone    = Med(wp_index);
%             
%             %max tones
%             maxNumTone_i     = find(numTone == max(numTone));
%             maxNumTone_index = wp_index(maxNumTone_i);
%             [wp_X_max wp_Y_max] = ind2sub(size(Med_Th_Dil),maxNumTone_index);
%             
%             if ~isempty(wp_X_max)
%                 %remap intervals into 0-100 range
%                 max_numTone = max(numTone);
%                 min_numTone = min(numTone);
%                 m_t = (100-0)/(max_numTone(1)-min_numTone(1));
%                 q_t = 100-m_t*(max_numTone);
%                 numTone_R = m_t.*numTone + q_t;
%                 
%                 %multiply for the numTone
%                 wp_X_Tone = wp_X.*numTone_R;
%                 wp_Y_Tone = wp_Y.*numTone_R;
%                 
%                 %find centroid coordinates (weighted on the gray tone)
%                 x = round(sum(wp_X_Tone)/sum(numTone_R));
%                 y = round(sum(wp_Y_Tone)/sum(numTone_R));
%             else
%                 x = NaN;
%                 y = NaN;
%                 wp_index = 1;
%                 numTone_R = 0;
%             end
%             
%             CM = zeros(512,512);
%             CM(wp_index) = numTone_R;
%             
%             %     subplot(sr,sc,inf)
%             %     imagesc(CM)
%             %     hold on
%             %     text(x,y, 'X','Color','red','FontSize',10);
%             
%             %     if inf == round(ImageSequence.CentFrame/ImageSequence.downsamplingfactor)
%             %         title('Cent Frame')
%             %     end
%             
%             X_Y_Weighted_Centroid     = [X_Y_Weighted_Centroid;[x y]];
%             X_Y_Weighted_Centroid_std = [X_Y_Weighted_Centroid_std; [std(  round(sum(wp_X_Tone)),[]) std(  round(sum(wp_Y_Tone)),[])]];
%             
%         end
%         close all
%         colorbar
%         
%         X_Y_Weighted_Centroid_DAY      = [X_Y_Weighted_Centroid_DAY    ; nanmedian(X_Y_Weighted_Centroid,1)];
%         X_Y_Weighted_Centroid_std_DAY  = [X_Y_Weighted_Centroid_std_DAY; [nanstd(X_Y_Weighted_Centroid(:,1),[]) nanstd(X_Y_Weighted_Centroid(:,2),[])]];
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT%
        X_Y_Weighted_Centroid       = [];
        X_Y_Weighted_Other_Centroid = [];
        sc = 6;
        sr = ceil(NumFrames/sc);
        H_Main_Cent_Seq_FigLabel = figure('Name',[Animal_name,'_Main_Centroid Thresholded_Image_Sequence day_',nday]);
        
        if SeqDirOKNO == 1
            for inf=1:NumFrames
                
                %median sequence for the current day
                Med        = squeeze(Med_Frames_Matrix(:,:,inf));
                Med_Th_Dil = squeeze(Med_Frames_Matrix_Th(:,:,inf));
                
                [bound,label] = bwboundaries(Med_Th_Dil);
                temp1 = regionprops(label,'Centroid');
                
                if ~isempty(temp1)
                    %%%
                    %%% Areas and Fluo values  %%% %%% %%% %%% %%% %%%
                    MetrixAreaLum = zeros(length(bound),1);
                    for i=1:length(bound)
                        currROI_ind_x_y = bound{i,1};
                        FluoValues = Med(sub2ind(size(Med),currROI_ind_x_y(:,1),currROI_ind_x_y(:,2)));
                        
                        SizeArea  = length(FluoValues);
                        if SizeArea<50
                            SizeArea = NaN;
                            FluoValues = NaN;
                        end
                        MaxLum    = max(FluoValues);
                        MeanLum   = mean(FluoValues);
                        MedianLum = median(FluoValues);
                        StdLum    = std(FluoValues,[]);
                        
                        histFluo = reshape(FluoValues,[],1);
                        [h_v h_i] = hist(histFluo);
                        
                        [h_v_max i_h_v_max] = max(h_v);
                        MostProbile_i_ok  = h_i(i_h_v_max(1));
                        
                        %metrics
                        MetrixAreaLum(i) = SizeArea*MostProbile_i_ok;
                    end
                    
                    [MaxMetrixAreaLum i_MaxMetrixAreaLum] = max(MetrixAreaLum);
                    
                    if ~isnan(MaxMetrixAreaLum)
                        CentroidSelected = temp1(i_MaxMetrixAreaLum,1).Centroid;
                    else
                        CentroidSelected = [NaN NaN];
                    end
                    
                    % centroid coordinates in the median image at frame "inf"
                    % -> main Centroid
                    x = round(CentroidSelected(1));
                    y = round(CentroidSelected(2));
                    
                    
                    
                    % check the other centroids
                    o_area_i = 0;
                    List_Areas = MetrixAreaLum;
                    List_Areas(isnan(List_Areas)) = 0; %la Metrica NaN è messa a 0 (per via del funzionamento di sortrows)
                    List_Areas = [List_Areas, [1:length(MetrixAreaLum)]'];
                    List_Areas = sortrows(List_Areas,[-1 -2]);
                    %scorre tutte le aree trovate (partendo dalla seconda perchè la prima è quella di riferimento) ordinate per metrica descrescente
                    for o_c=2:length(MetrixAreaLum)
                        
                        if List_Areas(o_c,1) ~= 0 %se non è 0, ovvero Metrica NaN
                            o_c_i = List_Areas(o_c,2);
                            Centroid_xy = temp1(o_c_i,1).Centroid;
                            %calcolo la dist tra il centroide di riferimento e il centroide dell'area corrente
                            dist_cent = sqrt((Centroid_xy(1)-x).^2+(Centroid_xy(2)-y).^2);
                            radius = 84; %pixels
                            %se la dist tra i centroidi è minore del raggio stimato allora lo considero come facente parte dell'area del centroide di ref
                            if dist_cent<=radius
                                %the same area
                            else
                                X_Y_Weighted_Other_Centroid = [X_Y_Weighted_Other_Centroid; [Centroid_xy(1) Centroid_xy(2)] ];
                            end
                        end
                    end
                    
                    
                else
                    x = NaN;
                    y = NaN;
                end
                % centroid coordinates in the median image for all frames of the image sequence (current day)
                X_Y_Weighted_Centroid     = [X_Y_Weighted_Centroid;[x y]];
                
                
                %%%%%%%%%
                Im = squeeze(Med_Frames_Matrix_Th(:,:,inf));
                subplot(sr,sc,inf)
                imagesc(Im);
                caxis([0 1])
                colormap gray
                if inf == round(ImageSequence.CentFrame/ImageSequence.downsamplingfactor)
                    title('Cent Frame')
                end
                hold on;
                plot(x,y,'r.','MarkerSize',10);
                H_Main_Cent_Seq_FigLabel_filename = [Animal_name,'_Main_Centroid_Seq_',nday];
                saveas(H_Main_Cent_Seq_FigLabel,H_Main_Cent_Seq_FigLabel_filename,'fig');
                %%%%%%%%%
                
            end
            
            
            %main area
            xmed = nanmedian(X_Y_Weighted_Centroid(:,1));
            ymed = nanmedian(X_Y_Weighted_Centroid(:,2));
            X_Y_Weighted_Centroid_DAY      = [X_Y_Weighted_Centroid_DAY    ; [str2num(nday), xmed, ymed]];
            
            %calcolo errore
            d=(X_Y_Weighted_Centroid(:,1)-xmed).^2+(X_Y_Weighted_Centroid(:,2)-ymed).^2;
            percd = prctile(d,[75]);
            Err = X_Y_Weighted_Centroid(d<percd,:);
            
            X_Y_Weighted_Centroid_std_DAY  = [X_Y_Weighted_Centroid_std_DAY; [str2num(nday), nanstd(Err(:,1),[]) nanstd(Err(:,2),[])]];
            
            
            %%% check other areas
            if size(X_Y_Weighted_Other_Centroid,1) > 3
                maxNumRegion = 3;
                [idx_centr] = kmeans(X_Y_Weighted_Other_Centroid,maxNumRegion);
                otherCentroids_med = [];
                otherCentroids_std = [];
                for i_ce=1:maxNumRegion
                    listCurrCentr = X_Y_Weighted_Other_Centroid((idx_centr==i_ce),:);
                    if length(listCurrCentr)>3
                        otherCentroid_curr_med = median(listCurrCentr,1);
                        otherCentroid_curr_std = std(listCurrCentr,[],1);
                        %calcolo la dist tra il centroide di riferimento e il centroide dell'area corrente trovato
                        dist_cent = sqrt((otherCentroid_curr_med (1)-xmed).^2+(otherCentroid_curr_med (2)-ymed).^2);
                        if dist_cent>radius
                            X_Y_Weighted_Centroid_Other_DAY         = [X_Y_Weighted_Centroid_Other_DAY;     [str2num(nday) otherCentroid_curr_med]];
                            X_Y_Weighted_Centroid_std_Other_DAY     = [X_Y_Weighted_Centroid_std_Other_DAY; [str2num(nday) otherCentroid_curr_std]];
                        end
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            X_Y_Weighted_Centroid_DAY      = [X_Y_Weighted_Centroid_DAY    ; [str2num(nday), NaN, NaN]];
            X_Y_Weighted_Centroid_std_DAY  = [X_Y_Weighted_Centroid_std_DAY; [str2num(nday), NaN, NaN]];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end

end



if 1
    
    %PLot scatter
    x = X_Y_Weighted_Centroid_DAY(:,2);
    y = X_Y_Weighted_Centroid_DAY(:,3);
    xe = X_Y_Weighted_Centroid_std_DAY(:,2);
    ye = X_Y_Weighted_Centroid_std_DAY(:,3);
    
    
    %%%%%%%% Figure 1 %%%%%%%%
    H_Main_Cent_Scatter = figure('Name',[Animal_name,'_Main_Centroid_Scatter_Along_Days']);
    subplot(121)
    %scatter with errors
    nD = length(x);
    %Make these defaults later:
    dotColor = [1 0.3 0.3]; % conservative pink
    yeColor = [0, 0.4, 0.8]; % bright navy blue
    xeColor = [0.35, 0.35, 0.35]; % not-too-dark grey
    dotSize = 23;
    % set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    set(gca, 'FontSize', 10);
    hold all;
    for i = 1:nD
        plot([(x(i) - xe(i)) (x(i) + xe(i))], [y(i) y(i)], 'Color', xeColor);
        plot([x(i) x(i)], [(y(i) - ye(i)) (y(i) + ye(i))], 'Color', yeColor);
    end
    scatter(x, y, dotSize, repmat(dotColor, nD, 1));
    set(gca,'Ydir','reverse')
    t = [X_Y_Weighted_Centroid_DAY(:,1)]; b = num2str(t); Na = cellstr(b);
    % t = [1:9,11:20]'; b = num2str(t); Na = cellstr(b);
    dx = 0.2; dy = 0.2;
    text(x+dx, y+dy, Na);
    xlim([0 512])
    ylim([0 512])
    
    %distance between time-sequence centroid
    d = [];
    d(1)= 0;
    for i=2:size(X_Y_Weighted_Centroid_DAY,1)
        d(i) = sqrt( (y(i)-y(i-1))^2 + (x(i)-x(i-1))^ 2);
    end
    d = d';
    subplot(122)
    bar(t,d)
    %saving figure
    H_Main_Cent_Scatter_filename = [Animal_name,'_Main_Centroid_Scatter_Along_Days'];
    saveas(H_Main_Cent_Scatter,H_Main_Cent_Scatter_filename,'fig');
    
    
    
    %%%%%%%% Figure 2 %%%%%%%%
    if ~isempty(X_Y_Weighted_Centroid_Other_DAY)
        
        H_All_Cent_Scatter = figure('Name',[Animal_name,'_All_Centroid_Scatter_Along_Days_']);
        subplot(121)
        hold all;
        for i = 1:nD
            plot([(x(i) - xe(i)) (x(i) + xe(i))], [y(i) y(i)], 'Color', xeColor);
            plot([x(i) x(i)], [(y(i) - ye(i)) (y(i) + ye(i))], 'Color', yeColor);
        end
        scatter(x, y, dotSize, repmat(dotColor, nD, 1));
        set(gca,'Ydir','reverse')
        t = [X_Y_Weighted_Centroid_DAY(:,1)]; b = num2str(t); Na = cellstr(b);
        % t = [1:9,11:20]'; b = num2str(t); Na = cellstr(b);
        dx = 0.2; dy = 0.2;
        text(x+dx, y+dy, Na);
        xlim([0 512])
        ylim([0 512])
        
        %other centroids
        xc = X_Y_Weighted_Centroid_Other_DAY(:,2);
        yc = X_Y_Weighted_Centroid_Other_DAY(:,3);
        xce = X_Y_Weighted_Centroid_std_Other_DAY(:,2);
        yce = X_Y_Weighted_Centroid_std_Other_DAY(:,3);
        nD = length(xc);
        %Make these defaults later:
        dotColor = [.5 .5 .5]; % conservative pink
        yeColor  = [0, 0.5, 0.5]; % bright navy blue
        xeColor  = [0, 0.5, 0.5]; % not-too-dark grey
        dotSize  = 23;
        % set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        set(gca, 'FontSize', 10);
        hold all;
        for i = 1:nD
            plot([(xc(i) - xce(i)) (xc(i) + xce(i))], [yc(i) yc(i)], 'Color', xeColor);
            plot([xc(i) xc(i)], [(yc(i) - yce(i)) (yc(i) + yce(i))], 'Color', yeColor);
        end
        scatter(xc, yc, dotSize, repmat(dotColor, nD, 1));
        t = [X_Y_Weighted_Centroid_Other_DAY(:,1)]; b = num2str(t); Na = cellstr(b);
        % t = [1:9,11:20]'; b = num2str(t); Na = cellstr(b);
        dx = 0.2; dy = 0.2;
        text(xc+dx, yc+dy, Na);
        
        %saving figure
        H_All_Cent_Scatter_filename = [Animal_name,'_All_Centroid_Scatter_Along_Days'];
        saveas(H_All_Cent_Scatter,H_All_Cent_Scatter_filename,'fig');
    end
    
    
    %saving
    save('X_Y_Weighted_Centroid_DAY','X_Y_Weighted_Centroid_DAY')
    save('X_Y_Weighted_Centroid_std_DAY','X_Y_Weighted_Centroid_std_DAY')
    save('X_Y_Weighted_Centroid_Other_DAY','X_Y_Weighted_Centroid_Other_DAY')
    save('X_Y_Weighted_Centroid_std_Other_DAY','X_Y_Weighted_Centroid_std_Other_DAY')
     
    
end



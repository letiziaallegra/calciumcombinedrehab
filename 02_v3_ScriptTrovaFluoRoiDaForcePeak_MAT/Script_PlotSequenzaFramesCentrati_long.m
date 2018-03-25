
%mappa media
MEAN = MatrixImageForcePeaks{end,1};
% rect_ROI = [5   113   504   392];
rect_ROI = [1   1   511   511];




%tutti i frames per tutti
if 0
    for k=1:size(MatrixImageForcePeaks,1)-1
        
        figure
        for j=1:size(MatrixImageForcePeaks{k,1},3)
            
            Im = MatrixImageForcePeaks{k,1}(:,:,j);
            
            Im_diff = (Im(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3))-...
                MEAN(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3)))./...
                MEAN(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3));
            
            %         Im_diff = (MatrixImageForcePeaks{k,1}(:,:,j) - MEAN)./MEAN;
            %         Im_diff = MatrixImageForcePeaks{k,1}(:,:,j) ;
            
            if j<6
                subplot(4,5,j)
            elseif j==6
                subplot(4,5,8)
            elseif j>6
                subplot(4,5,4+j)
            end
            Im_diff = medfilt2(Im_diff);
            imagesc(Im_diff);
                    colormap hot
        end
    end
end
 

%tutti i frames centrali
if 0
    figure
    for k=1:size(MatrixImageForcePeaks,1)-1
        j=6;
        
        
        Im = MatrixImageForcePeaks{k,1}(:,:,j);
        
        Im_diff = (Im(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3))-...
            MEAN(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3)))./...
            MEAN(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3));
        
        %     Im_diff = (MatrixImageForcePeaks{k,1}(:,:,j) - MEAN)./MEAN;
        %         Im_diff = MatrixImageForcePeaks{k,1}(:,:,j) ;
        
        subplot(5,5,k)
        imagesc(Im_diff);
    end
end


%mediane di tutti i frames
if 0
    figure
    MedMAtrix = zeros(rect_ROI(4)+1,rect_ROI(3)+1,16);
    for j=1:size(MatrixImageForcePeaks{1,1},3) %-> num frames
        
        ImBuff = zeros(rect_ROI(4)+1,rect_ROI(3)+1,size(MatrixImageForcePeaks{1,1},3));
        for k=1:size(MatrixImageForcePeaks,1)-1 %-> num picchi
            
            Im = MatrixImageForcePeaks{k,1}(:,:,j);
            
            Im_diff = (Im(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3))-...
                       MEAN(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3)))./...
                       MEAN(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3));
                 
            ImBuff(:,:,k) = Im_diff;                
                 
        end
        
        MedMAtrix(:,:,j) = median(ImBuff,3);
        if j<6
            subplot(4,5,j)
        elseif j==6
            subplot(4,5,8)
        elseif j>6
            subplot(4,5,4+j)
        end
        imagesc(MedMAtrix(:,:,j));
        colormap hot
        
        
        
        
    end
end


%mediane frames selezionati
% TrOK = [1:26];
% TrOK = [13];
% posPlot = [5 10 25 9];
numTr_REST = 5;
TrOK = [1:size(MatrixImageForcePeaks,1)-numTr_REST-1];
if 1
    figure
    MedMAtrix = zeros(rect_ROI(4)+1,rect_ROI(3)+1,size(MatrixImageForcePeaks{1,1},3));
    for j=1:size(MatrixImageForcePeaks{1,1},3) %-> num frames
        
        ImBuff = zeros(rect_ROI(4)+1,rect_ROI(3)+1,length(TrOK));
        
        for k=1:length(TrOK)%-> num picchi ok
            
            Im = MatrixImageForcePeaks{TrOK(k),1}(:,:,j);
            
            Im_diff = (Im(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3))-...
                       MEAN(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3)))./...
                       MEAN(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3));
                 
            ImBuff(:,:,k) = Im_diff;                
                 
        end
        
        MedMAtrix(:,:,j) = median(ImBuff,3);
        Med = (squeeze(MedMAtrix(:,:,j)));
        
        Med = medfilt2(Med);
        
        [value, location] = max(Med(:));
        [R,C] = ind2sub(size(Med),location);
        
%        Med(R,C) = 0;
%         [value, location] = max(Med(:));
%         [R,C] = ind2sub(size(Med),location);
% %         
%        Med(R,C) = 0;
%         [value, location] = max(Med(:));
%         [R,C] = ind2sub(size(Med),location);
        
%         if j<21
%             subplot(5,10,j)
%         elseif j==21
%             subplot(5,10,25)
%         elseif j>21
%             subplot(5,10,9+j)
%         end
        
        if j<11
            subplot(6,10,j)
        elseif j==11
            subplot(6,10,15)
        elseif j>11
            subplot(6,10,9+j)
        end
        
        
        
%         Med(R-5:R+5,C-5:C+5) = 1;
        imagesc(Med);
        colormap hot
        caxis([-0.2 0.2])
        
        
        
        
    end
end



%tutti i frames di rest
numTr_REST = 5;
TrOK = [size(MatrixImageForcePeaks,1)-numTr_REST-1:size(MatrixImageForcePeaks,1)-1];
if 0
    figure
    for k=1:length(TrOK)               
        
        Im = MatrixImageForcePeaks{TrOK(k),1};
        
        Im_diff = (Im(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3))-...
            MEAN(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3)))./...
            MEAN(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3));
        
        %     Im_diff = (MatrixImageForcePeaks{k,1}(:,:,j) - MEAN)./MEAN;
        %         Im_diff = MatrixImageForcePeaks{k,1}(:,:,j) ;
        
        subplot(1,5,k)
        imagesc((Im_diff));
        colormap hot
        caxis([-0.2 0.2])
    end
end

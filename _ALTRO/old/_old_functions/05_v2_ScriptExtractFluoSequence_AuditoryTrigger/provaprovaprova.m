  %PLOT%
        X_Y_Weighted_Centroid     = [];
        sc = 6;
        sr = ceil(NumFrames/sc);
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
                    
                    SizeArea  = length(FluoValues)
                    MaxLum    = max(FluoValues);
                    MeanLum   = mean(FluoValues);
                    MedianLum = median(FluoValues);
                    
                    %metrics
                    MetrixAreaLum(i) = SizeArea^MedianLum;
                end
                [MaxMetrixAreaLum i_MaxMetrixAreaLum] = max(MetrixAreaLum);
                CentroidSelected = temp1(i_MaxMetrixAreaLum,1).Centroid
                
                % centroid coordinates in the median image at frame "inf"
                x = round(CentroidSelected(1));
                y = round(CentroidSelected(2));
            else
                x = NaN;
                y = NaN;
            end
            % centroid coordinates in the median image for all frames of the image sequence (current day)
            X_Y_Weighted_Centroid     = [X_Y_Weighted_Centroid;[x y]];
            
        end
        X_Y_Weighted_Centroid_DAY      = [X_Y_Weighted_Centroid_DAY    ; nanmedian(X_Y_Weighted_Centroid,1)];
        X_Y_Weighted_Centroid_std_DAY  = [X_Y_Weighted_Centroid_std_DAY; [nanstd(X_Y_Weighted_Centroid(:,1),[]) nanstd(X_Y_Weighted_Centroid(:,2),[])]];
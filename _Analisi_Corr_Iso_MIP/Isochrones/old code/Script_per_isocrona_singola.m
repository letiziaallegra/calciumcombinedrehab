%carica il ImDelay_Day_M_All_NOMEANIMALE_control_MEDIA_TOT
%e poi lancia questo script
    
    
    IsoMean = nanmean(ImDelay_Day_M_All,3);
    
    IsoMean_Big = imresize(IsoMean, 40);
    
    figure
    imagesc(IsoMean_Big)
    changeLabel_MIP(size(IsoMean_Big,1),size(IsoMean_Big,2));
    caxis([-0.1 0.3])
    
    
    %
    % IsoMean_Big_F = medfilt2(IsoMean_Big, [50,50]);
    %
    % figure
    % imagesc(IsoMean_Big_F )
    % changeLabel_MIP(size(IsoMean_Big_F,1),size(IsoMean_Big_F,2));
    % caxis([-0.1 0.3])
    
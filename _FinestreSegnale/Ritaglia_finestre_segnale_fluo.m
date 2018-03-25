
folder_data = 'C:\Users\Stefano\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\_data_MAT_GCamp_Store';
figure


for an=1:4
    
    switch an
        case 1
            %animal CONTROL
            an_name = 'GCaMPChR2_17_control';
            dayL    = {'01','02','03','04','05'};
            treatStrg = 'control';
        case 2
            %animal STROKE
            an_name = 'GCaMPChR2_25_stroke';
            dayL    = {'02','03','04','05'};
            treatStrg = 'stroke';
        case 3
            %animal REHAB1W
            an_name = 'GCaMPChR2_15_stroke_BoNT';
            dayL    = {'01','02','03','04','05'};
            treatStrg = 'rehab1w';
        case 4
            %animal REHAB4W
            an_name = 'GCaMPChR2_15_stroke_BoNT';
            dayL    = {'17','18','19','20'};
            treatStrg = 'rehab4w';
        case 5
            %animal ROBOT4W
    end
    
    %loadGCamp
    sig_Wind = [];
    
    for i_d =1:length(dayL)
        
        load([folder_data,'\',an_name,'\','dataMouseGCamp_',an_name,'_',dayL{i_d},'_Par_MIPSIP_Par'])        
        
        %SIP rotrale
        i_Sig = 2;
        fs    = 100;
        int1  = 50;
        int2  = 300;
        t = [-int1:int2]/fs;
        
        %signal
        sig = dataGCamp.ROI_MIP_SIP.ROI_Signal(:,i_Sig);
        
        info_fluoPeak = dataGCamp.ROI_MIP_SIP_PeaksPar_Fx_Fluo{2,i_Sig};
        
        %start fluo peak
        st_fluoPeak  = info_fluoPeak(:,3)*fs;
        
        
        for i=2:length(st_fluoPeak)-1
            
            sig_Wind_buf = sig(st_fluoPeak(i)-int1:st_fluoPeak(i)+int2);
            
            sig_Wind = [sig_Wind, sig_Wind_buf];
            
        end
    end
    
    %plot animale
    sig_Wind_mean = (mean(sig_Wind,2))';
    sig_Wind_std  = (std(sig_Wind,[],2)/sqrt(size(sig_Wind,2)))';
    
    subplot(2,2,an)
    
    
    
    hold on
    fill(  [t fliplr(t)],  [ (sig_Wind_mean+sig_Wind_std)  (fliplr(sig_Wind_mean-sig_Wind_std))], [204/255 204/255 205/255])
    plot(t,sig_Wind_mean,'k','LineWidth',2)
    
    title(treatStrg)
    xlabel('time [s]')
    ylabel('deltaf/f')
    xlim([-0.5 3])
    ylim([-2 10])
    
    
end



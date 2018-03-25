%rapporto e differenza tMax rostral/ tMax caudal della FLUO (_Folder_Image_And_Stat_MEAN_MEAN_POOL)

C = 4; %rapporto delta_T_Rostral/delta_T_Caudal
C = 5; %delta_T_Rostral-delta_T_Caudal
C = 6; %delta_T_Caudal-delta_T_Rostral

M = [];
S = [];

for i=1:4
    
    switch i
        
        case 1            
            %control
            d = xlsread('DataFluoStat_Ctrl_Stk_Rehab_MaxT_SIP Rapporto_RostrCaud.xls',1);
            
        case 2            
            %stroke
            d  = xlsread('DataFluoStat_Ctrl_Stk_Rehab_MaxT_SIP Rapporto_RostrCaud.xls',2);
            
        case 3            
            %rehab1
            d = xlsread('DataFluoStat_Ctrl_Stk_Rehab_MaxT_SIP Rapporto_RostrCaud.xls',4);
        
        case 4            
            %rehab4
            d = xlsread('DataFluoStat_Ctrl_Stk_Rehab_MaxT_SIP Rapporto_RostrCaud.xls',3);
    end
    
    %il rapporto è nella quarta colonna
    %la differenza è nella quinta colonna
    d = d(:,C);
    md = nanmedian(d);
    sd = nanstd(d,[])/sqrt(length(d));
    
    M = [M, md];
    S = [S, sd];
end

subplot(321)
h = barwitherr(S , M);

set(gca,'FontSize',14)
title( 'Difference')
set(gca,'xTick',[1 2 3 4])
set(gca,'xTickLabel',{'ctrl','strk','reh1w','reh4w',})
xlabel('Week')

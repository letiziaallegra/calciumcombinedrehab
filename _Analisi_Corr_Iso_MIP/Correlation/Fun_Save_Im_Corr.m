%save CORR image
function Fun_Save_Im_Corr(Animal_name, SAVE_Dir_Method1_An, SAVE_Dir_Method2_An, SAVE_Dir_Method3_An, SEL_SEP, SEL_PART, Mat_corr_1, Mat_corr_2, Mat_corr_3, index)


%%%%%%%%%%%%%%%%%%%%%
if strcmp(SEL_SEP,'DAY')
    
    %%%%%%%%%%%%%%%%%%%%%
    %save CORR image
    figure
    imagesc(Mat_corr_1)
    caxis([0 1])
    colorbar
    filename_1 = [SAVE_Dir_Method1_An,'Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_1_DayN_',num2str(index),'_',Animal_name];
    saveas(gca,[filename_1,'Fig'],'fig')
    saveas(gca,[filename_1,'Fig'],'jpeg')
    close
    
    figure
    imagesc(Mat_corr_2)
    caxis([0 1])
    colorbar
    filename_2 = [SAVE_Dir_Method2_An,'Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_2_DayN_',num2str(index),'_',Animal_name];
    saveas(gca,[filename_2,'Fig'],'fig')
    saveas(gca,[filename_2,'Fig'],'jpeg')
    close
    
    figure
    imagesc(Mat_corr_3)
    caxis([0 1])
    colorbar
    filename_3 = [SAVE_Dir_Method3_An,'Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_3_DayN_',num2str(index),'_',Animal_name];
    saveas(gca,[filename_3,'Fig'],'fig')
    saveas(gca,[filename_3,'Fig'],'jpeg')
    close
        
    %%%%%%%%%%%%%%%%%%%%%
elseif strcmp(SEL_SEP,'WEEK')
    
    %%%%%%%%%%%%%%%%%%%%%
    %save CORR image
    figure
    imagesc(Mat_corr_1)
    caxis([0 1])
    colorbar
    filename_1 = [SAVE_Dir_Method1_An,'Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_1_WeekN_',num2str(index),'_',Animal_name];
    saveas(gca,[filename_1,'Fig'],'fig')
    saveas(gca,[filename_1,'Fig'],'jpeg')
    close
    
    figure
    imagesc(Mat_corr_2)
    caxis([0 1])
    colorbar
    filename_2 = [SAVE_Dir_Method2_An,'Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_2_WeekN_',num2str(index),'_',Animal_name];
    saveas(gca,[filename_2,'Fig'],'fig')
    saveas(gca,[filename_2,'Fig'],'jpeg')
    close
        
    figure
    imagesc(Mat_corr_3)
    caxis([0 1])
    colorbar
    filename_3 = [SAVE_Dir_Method3_An,'Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_3_WeekN_',num2str(index),'_',Animal_name];
    saveas(gca,[filename_3,'Fig'],'fig')
    saveas(gca,[filename_3,'Fig'],'jpeg')
    close
   
    
end


end
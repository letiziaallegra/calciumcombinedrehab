%save CORR image
function Fun_Save_Data_Corr(Animal_name, SAVE_Dir_Method1_An, SAVE_Dir_Method2_An, SEL_SEP, SEL_PART, Mat_corr_1, Mat_corr_2, index)

%functional areas names%
FuncAres_Name = {'1','2','3','4','5','6','7','8','9','10','11'};

%%%%%%%%%%%%%%%%%%%%%
if strcmp(SEL_SEP,'DAY')
    
    %%%%%%%%%%%%%%%%%%%%%
    %save CORR image
    figure
    imagesc(Mat_corr_1)
    set(gca,'yTick',[0.98:11.98])
    set(gca,'xTick',[0.98:11.98])
    set(gca,'xTickLabel',FuncAres_Name)
    set(gca,'yTickLabel',FuncAres_Name)
    colorbar
    filename_1 = [SAVE_Dir_Method1_An,'Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_1_DayN_',num2str(index-2),'_',Animal_name];
    saveas(gca,[filename_1,'Fig'],'fig')
    saveas(gca,[filename_1,'Fig'],'jpeg')
    close
    
    figure
    imagesc(Mat_corr_2)
    set(gca,'yTick',[0.98:11.98])
    set(gca,'xTick',[0.98:11.98])
    set(gca,'xTickLabel',FuncAres_Name)
    set(gca,'yTickLabel',FuncAres_Name)
    colorbar
    filename_2 = [SAVE_Dir_Method2_An,'Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_2_DayN_',num2str(index-2),'_',Animal_name];
    saveas(gca,[filename_2,'Fig'],'fig')
    saveas(gca,[filename_2,'Fig'],'jpeg')
    close
        
    %%%%%%%%%%%%%%%%%%%%%
elseif strcmp(SEL_SEP,'WEEK')
    
    %%%%%%%%%%%%%%%%%%%%%
    %save CORR image
    figure
    imagesc(Mat_corr_1)
    set(gca,'yTick',[0.98:11.98])
    set(gca,'xTick',[0.98:11.98])
    set(gca,'xTickLabel',FuncAres_Name)
    set(gca,'yTickLabel',FuncAres_Name)
    colorbar
    filename_1 = [SAVE_Dir_Method1_An,'Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_1_WeekN_',num2str(index),'_',Animal_name];
    saveas(gca,[filename_1,'Fig'],'fig')
    saveas(gca,[filename_1,'Fig'],'jpeg')
    close
    
    figure
    imagesc(Mat_corr_2)
    set(gca,'yTick',[0.98:11.98])
    set(gca,'xTick',[0.98:11.98])
    set(gca,'xTickLabel',FuncAres_Name)
    set(gca,'yTickLabel',FuncAres_Name)
    colorbar
    filename_2 = [SAVE_Dir_Method2_An,'Corr_Areas_',SEL_SEP,'_',SEL_PART,'_Meth_2_WeekN_',num2str(index),'_',Animal_name];
    saveas(gca,[filename_2,'Fig'],'fig')
    saveas(gca,[filename_2,'Fig'],'jpeg')
    close
   
    
end


end
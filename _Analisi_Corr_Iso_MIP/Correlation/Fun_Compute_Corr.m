%compute correlation coefficient between all the couples
function [Mat_corr_1 Mat_corr_2 Mat_corr_3] = Fun_Compute_Corr(FuncAreaValue_Day_Met1, FuncAreaValue_Day_Met2, FuncAreaValue_Day_Met3, Num_FunAreas)

Mat_corr_1 = zeros(Num_FunAreas,Num_FunAreas);
Mat_corr_2 = zeros(Num_FunAreas,Num_FunAreas);
Mat_corr_3 = zeros(Num_FunAreas,Num_FunAreas);

for fM_X = 1:Num_FunAreas
    
    SigX_1 = FuncAreaValue_Day_Met1(fM_X,:);
    SigX_2 = FuncAreaValue_Day_Met2(fM_X,:);
    SigX_3 = FuncAreaValue_Day_Met3(fM_X,:);
    
    for fM_Y = 1:Num_FunAreas
        
        SigY_1 = FuncAreaValue_Day_Met1(fM_Y,:);
        SigY_2 = FuncAreaValue_Day_Met2(fM_Y,:);
        SigY_3 = FuncAreaValue_Day_Met3(fM_Y,:);
        
        Corr_1 = corrcoef( SigX_1,  SigY_1);
        Corr_2 = corrcoef( SigX_2,  SigY_2);
        Corr_3 = corrcoef( SigX_3,  SigY_3);
        
        
        Mat_corr_1(fM_X,fM_Y) = Corr_1(1,2);
        Mat_corr_2(fM_X,fM_Y) = Corr_2(1,2);
        if length(Corr_3)>1
            Mat_corr_3(fM_X,fM_Y) = Corr_3(1,2);
        elseif ~isnan(Corr_3) 
            Mat_corr_3(fM_X,fM_Y) = Corr_3(1);
        end
        
    end
    
end

end
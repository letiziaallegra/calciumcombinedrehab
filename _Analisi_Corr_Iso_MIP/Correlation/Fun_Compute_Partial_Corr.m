%compute partial correlation coefficient between all the couples
function [Mat_corr_1 Mat_corr_2 Mat_corr_3] = Fun_Compute_Partial_Corr(FuncAreaValue_Day_Met1, FuncAreaValue_Day_Met2, FuncAreaValue_Day_Met3, Num_FunAreas)

%%%Partial correlation
Mat_corr_1 = partialcorr(FuncAreaValue_Day_Met1');
Mat_corr_2 = partialcorr(FuncAreaValue_Day_Met2');
%%%

%%%è uguale al caso non Partial perchè è una sequenza di numeri
Mat_corr_3 = zeros(Num_FunAreas,Num_FunAreas);

for fM_X = 1:Num_FunAreas
    
    SigX_3 = FuncAreaValue_Day_Met3(fM_X,:);
    
    for fM_Y = 1:Num_FunAreas
        
        SigY_3 = FuncAreaValue_Day_Met3(fM_Y,:);

        Corr_3 = corrcoef( SigX_3,  SigY_3);
        
        if length(Corr_3)>1
            Mat_corr_3(fM_X,fM_Y) = Corr_3(1,2);
        elseif ~isnan(Corr_3) 
            Mat_corr_3(fM_X,fM_Y) = Corr_3(1);
        end
        
    end
    
end
%%%

end
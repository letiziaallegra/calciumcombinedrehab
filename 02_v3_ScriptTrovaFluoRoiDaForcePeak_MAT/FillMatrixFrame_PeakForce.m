function Matrix = FillMatrixFrame_PeakForce(Matrix,folderTASK_FLUO,d,indexCurrentImage,FramePrePost,isMem)

vect = [indexCurrentImage-10:indexCurrentImage+FramePrePost];

% vect = [indexCurrentImage-FramePrePost:indexCurrentImage+FramePrePost];

for j=1:length(vect)
    
    nameImage = d(vect(j),1).name;
    Im = imread([folderTASK_FLUO,'\',nameImage]);
    
    %filter (salt and pepper removal)
    Im = medfilt2(Im);    
    
    if j==1
        Matrix{isMem,1} = zeros(size(Im,1),size(Im,2),length(vect));
    end      
    
    Matrix{isMem,1}(:,:,j) = Im;
    
end




end
nr_fr = 31;
% Initialize matrix using 'moviein'
frames = moviein(nr_fr);

MEAN = MatrixImageForcePeaks{end,1};

for i=1:nr_fr
    
    M = ((MatrixImageForcePeaks{1,1}(:,:,i)-MEAN)./MEAN)*100;
    
    imagesc(M)
    colormap hot
    frames(:, i) = getframe;
    
end

% Save the matrix so that this movie can be loaded later
save frames

movie(frames, 1, 2)
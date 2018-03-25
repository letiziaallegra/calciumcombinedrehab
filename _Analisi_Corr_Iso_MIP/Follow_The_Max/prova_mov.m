%%% Script to plot the average of the maximum value of fluo on the sequence 

ListAnimalTogether = {     'GCaMPChR2_26_stroke', 'GCaMPChR2_16_stroke_BoNT', 'GCaMPChR2_11_stroke_BoNT'};


%%%%%%%% User %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UserName = 'CNR-SSSUP';
UsbPort = 'I';
% UserName = 'Stefano';
% UsbPort = 'F';

%%%%%%%% Where SAVE Folder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PER FARLI INSIEME %%%%
Main_SAVE_Dir = [UsbPort,':\LENS\Max_Seq\'];

All_X = [];
All_Y = [];

for j=1:size(ImageSequence.MatrixImageSequence,1)
    
    I = ImageSequence.MatrixImageSequence{j,1};
    
    Tot_X = [];
    Tot_Y = [];
    
%     for im=1:size(I,3)
        for im=5:30
        
        A = squeeze(I(:,:,im));
        
        [maxA,ind] = max(A(:));
        [X,Y] = ind2sub(size(A),ind);
        
        Tot_X = [Tot_X; X];
        Tot_Y = [Tot_Y; Y];
        
    end
    
    All_X = [All_X Tot_X];
    All_Y = [All_Y Tot_Y];
    
end

MX = mean(All_X,2);
MY = mean(All_Y,2);

scatter(MX,MY)
hold on
for t=1:length(MX)
    
    text(MX(t),MY(t),num2str(t))
    
end
xlim([0 512])
ylim([0 512])

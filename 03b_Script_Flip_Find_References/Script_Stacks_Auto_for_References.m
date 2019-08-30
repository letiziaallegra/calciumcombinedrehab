%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Costruisci matrice coninfo per allineare (roto-traslazione)%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%le figure sono salvate come matrici .mat, 1 per Immagine in una cartella 

clear
close all
clc

User = getenv('username');
MainDir = ['C:\Users\asus\Desktop\to do'];
MainDir = '/Users/alessandro/Desktop/ELABORAZIONE DATA/Script_Flip_Find_References/Storage_Folder/to do';
% MainDir       = [UsbPort,':\LENS\Script_Flip_Find_References'];
%cartella dove mettere i dati 
%WorkingDir    = [MainDir,'\Working_Folder'];
WorkingDir    = [MainDir];
MainDirFolder = dir(WorkingDir);

%where placing BREGMA (y)
Bregma_y = 0.25; % (mm)
%size of the cranial window (mm)
Size_CW = 4.4; %(mm)
% Size_CW = 5.25; %(mm)

FLIP = true;

for mdf=3:length(MainDirFolder) %for MainDirFolder
    
    AnimalDir          = MainDirFolder(mdf,1).name;
    AnimalDirDayFolder = dir([WorkingDir,filesep,AnimalDir]);
    trans_all  = [];
    degree_all = [];
    indexDay = 0;
    
    for adf=3:length(AnimalDirDayFolder)  %for AnimalDirFolder
        
        indexDay = indexDay+1;
        DayImage     = AnimalDirDayFolder(adf,1).name;
        DayDirFolder = dir([WorkingDir,filesep,AnimalDir,filesep,DayImage]);    
        indexDay = str2num(DayImage(1:2));
        Im16 = imread([WorkingDir,filesep,AnimalDir,filesep,DayImage]);
        
        %16 bits -> 8 bits
        M16 = 2^16-1;
        M8 =  2^8-1;
        Im8 = uint8(Im16 * (M8/M16));
        
        %flip vertically
        if FLIP
            Im8_fv = flipdim(Im8 ,1);
        else
            Im8_fv = Im8;
        end
        Im = Im8_fv;    
        
        rw = size(Im,1);
        cl = size(Im,2);
        
        
        %set origin (based on black dot)
        H_MIP_Fig = figure('Name',['Animal:',AnimalDir,'Original Image',DayImage]);
        subplot(221)
        imagesc(Im)
        colormap gray
        
        
        %%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %rotation
        [X_gi Y_gi] = ginput;
        %coordinates of Bregma
        Xb = round(X_gi(1));
        Yb = round(Y_gi(1));
        %coordinates of posterior dot (along the medial line)
        Xr = round(X_gi(end));
        Yr = round(Y_gi(end));
        
        %compute the angle
        l_l = Yr-Yb;
        b_l = Xr-Xb;
        d_l = sqrt( (Yr-Yb)^2 + (Xr-Xb)^2 );
        
        degree = acosd(l_l/d_l) * ((b_l>0) * -1 + (b_l<=0) * 1);
        degree_all = [degree_all; [indexDay degree]];
        
        if degree ~= 0
            Im_R = imrotate(Im,degree,'crop');
        else
            Im_R = Im;
        end
        
        %new coordinates of Xb and Yb after rotation
        %rotation around z axis (with origin in the center of the image)
        [XYb_rot] = round([cosd(degree), -sind(degree); sind(degree), cosd(degree)] *  [Yb-round(cl/2); Xb-round(rw/2)]);
        XYb_rot2 =  XYb_rot + [round(cl/2); round(rw/2)];
        Xb_rot = XYb_rot2(2);
        Yb_rot = XYb_rot2(1);
        
        %new coordinates of Bregma
        XbN = 1;
        YbN = (rw*Bregma_y / Size_CW); %==0.250 mm
        
        
        %%%%%%%%%%%%%%
        subplot(221)
        hold on
        text(Xb,Yb, 'X','Color','red','FontSize',10);
        changeLabel_MIP(rw,cl)
        %%%%%%%%%%%%%%
        subplot(222)
        imagesc(Im_R)
        hold on
        text(Xb_rot,Yb_rot, 'X','Color','red','FontSize',10);
        changeLabel_MIP(rw,cl)
        %%%%%%%%%%%%%%
        
        Im_OR   = zeros(rw,cl);
        Im_OR_2 = zeros(rw,cl);
        
        %%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %translation along rows
        if Yb-YbN>0
            Im_OR(1:rw-abs(Yb_rot-YbN)+1,:) = Im_R(abs(Yb_rot-YbN):end,:);
            trans_lr = -1;
        elseif Yb-YbN<0
            Im_OR(abs(Yb_rot-YbN):end,:) = Im_R(1:rw-abs(Yb_rot-YbN)+1,:);
            trans_lr = +1;
        end
        
        %%%%%%%%%%%%%%
        subplot(223)
        if sum(sum(Im_OR)) == 0
            %ciao = 1
        else
            imagesc(Im_OR)
        end
        hold on        
        changeLabel_MIP(rw,cl)
        %%%%%%%%%%%%%%
        
        %translation along columns
        Im_OR_2(:,1:cl-abs(Xb_rot-XbN)+1) = Im_OR(:,abs(Xb_rot-XbN):end);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%
        subplot(224)
        if sum(sum(Im_OR_2)) == 0
            %ciao = 1
        else
            imagesc(Im_OR_2)
        end
        hold on
        changeLabel_MIP(rw,cl)
        %%%%%%%%%%%%%%
        
        pause
                
        
        %                        if left/right transl  translation Y    translation X
        trans_all = [trans_all; [trans_lr             abs(Yb_rot-YbN)   abs(Xb_rot-XbN)]];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
              
        
        
    end
    
        rot_transl = [degree_all trans_all];
        % Animal_name = AnimalDir(1:end-4)
        Animal_name = AnimalDir
        save([MainDir,filesep,Animal_name,'_Rot_Trans_Par'],'rot_transl')

end

display('END PROCESS')
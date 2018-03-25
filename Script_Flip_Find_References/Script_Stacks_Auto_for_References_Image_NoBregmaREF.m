%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Costruisci matrice con info per allineare (roto-traslazione)
%%% (questo script lavora su immagini in cui non era stato indicata la
%%%  posizione di bregma e gli altri due punti ortogonali)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%le figure sono salvate come matrici .mat, 1 per Immagine in una cartella 

clear
close all
clc
%% indicare l'immagine presa come ref %%
Im_REF = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
User = getenv('username');
MainDir = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\Script_Flip_Find_References'];
% MainDir       = [UsbPort,':\LENS\Script_Flip_Find_References'];
WorkingDir    = [MainDir,'\Working_Folder'];
MainDirFolder = dir(WorkingDir);

for mdf=3:length(MainDirFolder) %for MainDirFolder
    
    AnimalDir          = MainDirFolder(mdf,1).name;
    AnimalDirDayFolder = dir([WorkingDir,'\',AnimalDir]);
    trans_all          = [];
    degree_all         = [];
    indexDay           = 0;
    
    %immagine presa come ref        
    indexREF     = Im_REF+2;
    DayImage     = AnimalDirDayFolder(indexREF,1).name;
    DayDirFolder = dir([WorkingDir,'\',AnimalDir,'\',DayImage]);
    Im16 = imread([WorkingDir,'\',AnimalDir,'\',DayImage]);
    
    %16 bits -> 8 bits
    M16 = 2^16-1;
    M8 =  2^8-1;
    Im8 = uint8(Im16 * (M8/M16));
    
    %flip vertically
    Im8_fv = flipdim(Im8 ,1);
    Im = Im8_fv;
    
    rw = size(Im,1);
    cl = size(Im,2);
    
    %Point1
    H_MIP_Fig = figure('Name',['Original Image of the Day ',DayImage,' taken as frame of reference']);
    imagesc(Im)
    ImREF_plot = Im;
    colormap gray
    [X_gi_REF Y_gi_REF] = ginput;
    
    
    
    %%%%%%%%% Find angle between the two  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %rotation
    %coordinates of the First Point selected
    Xr = round(X_gi_REF(1));
    Yr = round(Y_gi_REF(1));
    %coordinates of the Second dot
    Xb = round(X_gi_REF(end));
    Yb = round(Y_gi_REF(end));
    
    %compute the angle
    l_l = Yr-Yb;
    b_l = Xr-Xb;
    d_l = sqrt( (Yr-Yb)^2 + (Xr-Xb)^2 );
    
    degreeREF = acosd(l_l/d_l) * ((b_l>0) * -1 + (b_l<=0) * 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    for adf=3:length(AnimalDirDayFolder)  %for AnimalDirFolder
        
        indexDay = indexDay+1;
        
                
        DayImage     = AnimalDirDayFolder(adf,1).name;
        DayDirFolder = dir([WorkingDir,'\',AnimalDir,'\',DayImage]);    
        
        Im16 = imread([WorkingDir,'\',AnimalDir,'\',DayImage]);
        
        %16 bits -> 8 bits
        M16 = 2^16-1;
        M8 =  2^8-1;
        Im8 = uint8(Im16 * (M8/M16));
        
        %flip vertically
        Im8_fv = flipdim(Im8 ,1);
        Im = Im8_fv;    
        
        rw = size(Im,1);
        cl = size(Im,2);
        
        
        %set origin (based on black dot)
        H_MIP_Fig = figure('Name',['Original Image',DayImage]);
        subplot(231)
        imagesc(Im)
        colormap gray
        
        if indexDay ~= Im_REF
            %Point
            [X_gi Y_gi] = ginput;
        else
            %Point1
            X_gi = X_gi_REF;
            Y_gi = Y_gi_REF;
        end
        
        
        
        %%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %rotation
%         [X_gi Y_gi] = ginput;
        %coordinates of the First Point selected
        Xb = round(X_gi(1));
        Yb = round(Y_gi(1));
        %coordinates of the Second dot
        Xr = round(X_gi(end));
        Yr = round(Y_gi(end));
        
        %compute the angle
        l_l = Yb-Yr;
        b_l = Xb-Xr;
        d_l = sqrt( (Yr-Yb)^2 + (Xr-Xb)^2 );
        
        degree = acosd(l_l/d_l) * ((b_l>0) * -1 + (b_l<=0) * 1);
        
        
        if degreeREF~=degree
            
            degree = (degree-degreeREF);            
            Im_R = imrotate(Im,degree,'crop');
            
        else
            degree = 0;
            Im_R = Im;
            
        end
        degree_all = [degree_all; [indexDay degree]];
        
        
        %new coordinates of Xb and Yb after rotation
        %rotation around z axis (with origin in the center of the image)
        [XYb_rot] = round([cosd(degree), -sind(degree); sind(degree), cosd(degree)] *  [Yb-round(cl/2); Xb-round(rw/2)]);
        XYb_rot2 =  XYb_rot + [round(cl/2); round(rw/2)];
        Xb_rot = XYb_rot2(2);
        Yb_rot = XYb_rot2(1);
        
        %new coordinates of Bregma
        XbN = round(X_gi_REF(1));
        YbN = round(Y_gi_REF(1)); %==0.250 mm
        
        
        %%%%%%%%%%%%%%
        subplot(232)
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
        elseif Yb==YbN
            Im_OR = Im_R;
            trans_lr = 0;
        end
        
        %%%%%%%%%%%%%%
        subplot(233)
        imagesc(Im_OR)
        hold on        
        changeLabel_MIP(rw,cl)
        %%%%%%%%%%%%%%
        
        %translation along colums
        if Xb-XbN>0
            Im_OR_2(:,1:cl-abs(Xb_rot-XbN)+1) = Im_OR(:,abs(Xb_rot-XbN):end);
            trans_ud = -1;
        elseif Xb-XbN<0
            Im_OR_2(:,abs(Xb_rot-XbN):end) = Im_OR(:,1:cl-abs(Xb_rot-XbN)+1);
            trans_ud = +1;
        elseif Xb==XbN
            Im_OR_2 = Im_OR;
            trans_ud = 0;
        end
        
        
%         %translation along columns
%         Im_OR_2(:,1:cl-abs(Xb_rot-XbN)+1) = Im_OR(:,abs(Xb_rot-XbN):end);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%
        
        subplot(234)
        imagesc(Im_OR_2)
        hold on
        changeLabel_MIP(rw,cl)
        %%%%%%%%%%%%%%
        
        
        
        %                        if left/right transl  translation Y    translation X
        trans_all = [trans_all; [trans_lr             abs(Yb_rot-YbN)   abs(Xb_rot-XbN)  trans_ud]];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %plot il frame ref
        subplot(235)
        imagesc(ImREF_plot) 
        
        pause
        
    end
    
        rot_transl = [degree_all trans_all];
        Animal_name = AnimalDir(1:end-4)
        save([MainDir,'\',Animal_name,'_Rot_Trans_Par'],'rot_transl')

end

display('END PROCESS')
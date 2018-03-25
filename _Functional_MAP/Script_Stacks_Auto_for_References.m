%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Costruisci matrice coninfo per allineare (roto-traslazione)%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load('Functional_MAP')

        Im = Functional_MAP;   
        
        rw = size(Im,1);
        cl = size(Im,2);
        
        
        %set origin (based on black dot)
        H_MIP_Fig = figure;
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
        
        degree = 0;
%         acosd(l_l/d_l) * ((b_l>0) * -1 + (b_l<=0) * 1);
        degree_all = [degree];
        
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
        YbN = 29; %==0.250 mm
        
        
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
        
        Im_OR   = Im_R; 
%         zeros(rw,cl);
        Im_OR_2 = Im_R;
%         
        %%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %translation along rows
        if Yb-YbN>0
            Im_OR(1:rw-abs(Yb_rot-YbN)+1,:) = Im_R(abs(Yb_rot-YbN):end,:);
            Im_OR(rw-abs(Yb_rot-YbN)+2:end,:) = [255];
            trans_lr = -1;
        elseif Yb-YbN<0
            Im_OR(abs(Yb_rot-YbN):end,:) = Im_R(1:rw-abs(Yb_rot-YbN)+1,:);
            trans_lr = +1;
        end
        
        %%%%%%%%%%%%%%
        subplot(223)
        imagesc(Im_OR)
        hold on        
        changeLabel_MIP(rw,cl)
        %%%%%%%%%%%%%%
        
        %translation along columns
        Im_OR_2(:,1:cl-abs(Xb_rot-XbN)+1,:) = Im_OR(:,abs(Xb_rot-XbN):end,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%
        subplot(224)
        imagesc(Im_OR_2)
        hold on
        changeLabel_MIP(rw,cl,1)
        %%%%%%%%%%%%%%
        
        pause
                
        
        %                        if left/right transl  translation Y    translation X
        trans_all = [[trans_lr             abs(Yb_rot-YbN)   abs(Xb_rot-XbN)]];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
        rot_transl = [degree_all trans_all];
        save('rot_transl_F','rot_transl')
        Im_Func_Mask = Im_OR_2;
        save('Im_Func_Mask','Im_Func_Mask')


display('END PROCESS')
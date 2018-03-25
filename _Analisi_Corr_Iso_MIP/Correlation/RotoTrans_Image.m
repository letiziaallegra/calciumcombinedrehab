%rotation and translation of the image
function Frame_OR_2 = RotoTrans_Image( Frame_ith, rot_transl, nday, rw, cl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% same rotation and translation %%%%%%%
i_day_actual = find(rot_transl(:,1) == str2num(nday));

%%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
degree = rot_transl(i_day_actual,2);
if  degree~= 0
    Frame_R = imrotate(Frame_ith,degree,'crop');
else
    Frame_R = Frame_ith;
end

Frame_OR = zeros(rw,cl);
Frame_OR_2 = zeros(rw,cl);

%%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%left/right transl
trans_lr  = rot_transl(i_day_actual,3);
%translation Y along rows
trans_Y   = rot_transl(i_day_actual,4);
%translation X along columns
trans_X   = rot_transl(i_day_actual,5);

%translation along rows
if trans_lr<0
    Frame_OR(1:rw-trans_Y+1,:) = Frame_R(trans_Y:end,:);
    trans_lr = -1;
elseif trans_lr>0
    Frame_OR(trans_Y:end,:) = Frame_R(1:rw-trans_Y+1,:);
    trans_lr = +1;
end

%translation along columns
Frame_OR_2(:,1:cl-trans_X+1) = Frame_OR(:,trans_X:end);

%new coordinates of Bregma
XbN = 1;
YbN = 29; %==0.250 mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
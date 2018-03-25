%load image
clear
close all
clc

load('Im_Func_Mask.mat')

Im_Func_Mask_Buf = Im_Func_Mask;
[s1 s2 s3] = size(Im_Func_Mask);

%
Thesh_Tot = [201, 192, 110, 200, 137, 196, 127, 128, 158, 164, 134]';

for i=1:11

figure
subplot(221)
imagesc(Im_Func_Mask_Buf)

%matrix
Im_gray = rgb2gray(Im_Func_Mask);
subplot(222)
imagesc(Im_gray)
colormap gray

%Im_th
% [x y] = ginput;
% Thesh = Im_gray(round(y),round(x));
% Thesh_Tot = [Thesh_Tot; Thesh];

Thesh = Thesh_Tot(i);
Im_Thesh_original = Im_gray==Thesh;

subplot(223)
imagesc(Im_Thesh_original)

if i==3 %shoulder
    Im_Thesh = Im_Thesh_original;
    se = strel('disk',4);
    Im_Thesh = imdilate(Im_Thesh,se);
else
    %erode
    se = strel('disk',4);
    Im_Thesh = imerode(Im_Thesh_original,se);
    %fill hole
    Im_Thesh =  imfill(Im_Thesh,'holes');
    %dilate
    %erode
    se = strel('disk',4);
    Im_Thesh = imdilate(Im_Thesh,se);
end


subplot(224)
imagesc(Im_Thesh)

[bound,label] = bwboundaries(Im_Thesh);
% temp1 = regionprops(label,'Centroid');
Areas = regionprops(label,'PixelList');

%%% Max Area
[n in] = max(cellfun('length',bound));

RegionArea_FunctionalMask_Index{i,1} = Areas(in,1).PixelList;
RegionArea_FunctionalMask_Index{i,2} = i;
RegionArea_FunctionalMask_Index{i,3} = [bound{in,1}(:,2)  bound{in,1}(:,1)];

end
save('RegionArea_FunctionalMask_Index','RegionArea_FunctionalMask_Index')
close all

%check plot
for i_Mask=1:11
    
    Im_MASK = zeros(s1,s2);
    Im_MASK(sub2ind([s1 s2],RegionArea_FunctionalMask_Index{i_Mask,1}(:,2),RegionArea_FunctionalMask_Index{i_Mask,1}(:,1))) = 1;
%     Im_MASK = imfill(Im_MASK);

    subplot(3,4,i_Mask)
    imagesc(Im_MASK)
    
end
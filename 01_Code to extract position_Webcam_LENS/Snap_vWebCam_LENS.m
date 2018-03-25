function [P_Marker, rect_Marker, rect_LED] = Snap_vWebCam_LENS(I)
close all

fM = figure('Name','First Frame','NumberTitle','off');
imshow(I);

%MARKER
%rect_ROI_default_pos =   [7   152];
rect_ROI_default_pos =   [7   200];
%ROI_size = [352 158];
ROI_size = [352 185];
% %nel caso di LED fisso non visibile 
% rect_ROI_default_pos =   [1   1];
% ROI_size = [639 479];

%area di movimento del marker
h_rect = imrect(gca, [rect_ROI_default_pos(1) rect_ROI_default_pos(2) ROI_size(1) ROI_size(2)]); %-> finestra di 150x150 o più grande px
pause

rect_ROI = round(getPosition(h_rect));
ROI_size             = rect_ROI(1,[1 2]);
rect_ROI_default_pos = rect_ROI(1,[3 4]);

rect_Marker = round(getPosition(h_rect));

I_Marker = I(rect_Marker(2):rect_Marker(2)+rect_Marker(4),rect_Marker(1):rect_Marker(1)+rect_Marker(3));
[bound,label] = bwboundaries(I_Marker);

%if it finds the region
if ~isempty(bound)    
    Cen = regionprops(label,'Centroid');
    P_Marker = round([Cen(1,1).Centroid(1) Cen(1,1).Centroid(2)]);
else
    error('marker not found')
end    

%LED
%rect_LED_default_pos =   [491   277];
%LED_size = [137   171];
rect_LED_default_pos =   [330   237];
LED_size = [290   211];

%%ok
rect_LED = [rect_LED_default_pos LED_size];
imshow(I);
h_rect = imrect(gca, [rect_LED]); 
rect_LED = round(getPosition(h_rect));


% %%per videowebcam piccoli
% imshow(I);
% h_rect = imrect(gca,  [rect_ROI_default_pos(2) rect_ROI_default_pos(2) ROI_size(1) ROI_size(2)])%-> finestra di 150x150 o più grande px
% pause
% rect_LED = round(getPosition(h_rect));

close(fM)

end
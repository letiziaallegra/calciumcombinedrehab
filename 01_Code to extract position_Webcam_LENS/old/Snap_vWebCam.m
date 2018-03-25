function [PLED, P, rect_Marker, rect_LED, BW] = Snap_vWebCam(bw,G)
close all

figure('Name','Indicate the MARKER Rectangular Area','NumberTitle','off')
imshow(G);
[I_Marker rect_Marker] = imcrop;
rect_Marker = round(rect_Marker);

close
figure('Name','Indicate the Marker position','NumberTitle','off')
% imshow(bw)
imshow(I_Marker)
[Px,Py] = ginput
P = round([Px(1,1) Py(1,1)]);
BW = double(I_Marker(P(2), P(1)))/255;
%updated del 10%
BW = BW*1.20;

close
figure('Name','Indicate the LED Rectangular Area','NumberTitle','off')
imshow(G);
[I_LED rect_LED] = imcrop;
rect_LED = round(rect_LED);


close
figure('Name','Indicate the LED position','NumberTitle','off')
imshow(I_LED);
[PxLED, PyLED] = ginput
PLED = [PxLED(1,1) PyLED(1,1)];

close







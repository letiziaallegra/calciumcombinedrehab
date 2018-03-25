function res=viewSingleFrame(filename,framenumber,strplot,mode,threshold)

filtwind=6;

 obj = mmreader(filename);
 if(strcmp(filename(end-3:end),'.avi'))
    totalframes = obj.NumberOfFrames;
 else
    totalframes = 10000;
 end
 
 if( framenumber <=  totalframes)     
 	video = read(obj,framenumber);
    
    grayImm = rgb2gray(video);  
    if(strcmp(mode,'gray'))
        filtImm = medfilt2(grayImm,[filtwind filtwind]);
        res =filtImm;
        if(strcmp(strplot,'plot'))
            figure,
            imshow(filtImm);
        end
    elseif(strcmp(mode,'bw'))
        filtImm = medfilt2(grayImm,[filtwind filtwind]);
%         threshold = 50/255;
        bw = im2bw(filtImm,threshold);
        se = strel('disk',3);
        imdil = imdilate(bw,se);
        imdiler = imerode(imdil,se);
        bw = imdiler;
        if(strcmp(strplot,'plot'))
            figure,
            imshow(bw);
        end
        res =bw;
    end
%     lenx = 40;
%     leny = 20;
%     Px = [291 276];
%     Mx = bw(round(Px(2))-leny:round(Px(2))+leny, round(Px(1))-lenx:round(Px(1))+lenx);
%     [bound,label] = bwboundaries(Mx);
%     temp1 = regionprops(label,'Centroid');
%     res1 = temp1(2,1).Centroid;
%     Px = Px+res1-[lenx leny];
%     hold on
%     scatter(Px(1),Px(2),'r')
 end
 
 return
function P = findPositionMarker(bw)

for i=1:size(bw,1)
    for j=1:size(bw,2)
        if(bw(i,j)==true)
            BW(i,j)=1;
        else
            BW(i,j)=0;
        end
    end
end
clear bw
load marker

len = floor(size(marker,1)/2);
for i=len+1:size(BW,1)-len
    for j=len+1:size(BW,2)-len
        bw = BW(i-len:i+len,j-len:j+len);
%         MSE(i-len,j-len)=sum(sum(bw-MARKER).^2);
        MSE(i-len,j-len)=sum(sum( (bw-marker).^2 ));
    end
end

m = min(min(MSE));
[i j] = find(MSE==m);
P(2)=mean(i)+len;
P(1)=mean(j)+len;

hold on
plot(P(1),P(2),'or','MarkerSize',40)
return


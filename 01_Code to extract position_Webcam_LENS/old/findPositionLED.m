function PLED = findPositionLED(GrayImm, expected_pos)
% function PLED = findPositionLED(GrayImm)

load LED

len = 40;%floor(size(LED,1)/2);

MSE = -1*ones(2*floor(len/4));

% expected_pos = [119 63];

% for i= len+1 : size(GrayImm,1)-len
%     for j= len+1 : size(GrayImm,2)-len
%         g = GrayImm(i-len:i+len,j-len:j+len);
%         MSE(i-len,j-len) = sum(sum( (LED-g).^2 ));
%     end
% end

found=false;
for i= expected_pos(1)-floor(len/4) : expected_pos(1)+floor(len/4)
    for j= expected_pos(2)-floor(len/4) : expected_pos(2)+floor(len/4)
        g = GrayImm(i-len:i+len,j-len:j+len);
        MSE(i-(expected_pos(1)-floor(len/4))+1,j-(expected_pos(2)-floor(len/4))+1) = sum(sum( (LED-g).^2 ));
        m = min(min(MSE));
        if(m<3) && (m>=0)
            found=true;
        end
    end
end
if(~found)
    for i= len+1 : size(GrayImm,1)-len
        for j= len+1 : size(GrayImm,2)-len
            g = GrayImm(i-len:i+len,j-len:j+len);
            MSE(i-len,j-len) = sum(sum( (LED-g).^2 ));
            m = min(min(MSE));
        end
    end
end

[r,c] = find(MSE==m);
PLED(2)=mean(r)+expected_pos(1)-floor(len/4);
PLED(1)=mean(c)+expected_pos(2)-floor(len/4);

imshow(GrayImm)
hold on
plot(PLED(1),PLED(2),'or','MarkerSize',40)

return


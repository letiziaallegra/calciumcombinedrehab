
%find white pixels
[wp_index] = find(Med_Th_Dil==1);
[wp_X wp_Y] = ind2sub(size(Med_Th_Dil),wp_index);

%find gray values corresponding to white pixels
numTone    = Med(wp_index);

%max tones
maxNumTone_i     = find(numTone == max(numTone));
maxNumTone_index = wp_index(maxNumTone_i);
[wp_X_max wp_Y_max] = ind2sub(size(Med_Th_Dil),maxNumTone_index);

%remap intervals into 0-100 range
max_numTone = max(numTone);
min_numTone = min(numTone);
m_t = (100-0)/(max_numTone(1)-min_numTone(1));
q_t = 100-m_t*(max_numTone);
numTone_R = m_t.*numTone + q_t;


%multiplipy for the numTone
wp_X_Tone = wp_X.*numTone_R;
wp_Y_Tone = wp_Y.*numTone_R;

%find centroid coordinates (weighted on the gray tone)
x = round(sum(wp_X_Tone)/sum(numTone_R));
y = round(sum(wp_Y_Tone)/sum(numTone_R));

figure
% subplot(121)
% imagesc(Med_Th_Dil)
% subplot(122)
CM = zeros(512,512);
CM(wp_index) = numTone_R;
hold on
imagesc(CM)
text(x,y, 'X','Color','red','FontSize',15);
ciao=1
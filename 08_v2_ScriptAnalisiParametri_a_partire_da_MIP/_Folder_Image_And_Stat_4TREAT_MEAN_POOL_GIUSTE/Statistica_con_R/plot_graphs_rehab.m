clear
figure

for j=1:3
    
    switch j
        case 1
            Par = load('max.csv','-ASCII');
            Te = 'amplitude';
        case 2
            Par = load('slope.csv','-ASCII');
            Te = 'slope';
        case 3
            Par = load('deltaT.csv','-ASCII');
            Te = 'delay';
    end

MEAN = [];
ERRSTD = [];


for i=1:4
    
    index = find(Par(:,2)==i);
    MEAN_b = mean(Par(index,1));
    STD_b    = std(Par(index,1),[]);
    ERRSTD_b = STD_b/sqrt(length(Par(index,1)));
    
    MEAN   = [MEAN MEAN_b];
    ERRSTD = [ERRSTD ERRSTD_b];
    
end

subplot(1,3,j)
bw = barwitherr(ERRSTD, MEAN);
title(Te);
xticks([1:4])
xticklabels({'C','S','R4w','R1w'})
xlabel('Week')

end